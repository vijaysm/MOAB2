#include "MBInterface.hpp"
#include "MBParallelComm.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBReadUtilIface.hpp"
#include "EntitySequenceManager.hpp"
#include "EntitySequence.hpp"
#include "TagServer.hpp"
#include "MBTagConventions.hpp"
#include "MBSkinner.hpp"
#include "MBParallelConventions.h"
#include "MBCore.hpp"
#include "MBError.hpp"

#define MAX_SHARING_PROCS 10  

extern "C" 
{
#include "gs.h"
#include "tuple_list.h"
}

#include <assert.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

#define INITIAL_BUFF_SIZE 1024

#define PACK_INT(buff, int_val) {int tmp_val = int_val; PACK_INTS(buff, &tmp_val, 1);}

#define PACK_INTS(buff, int_val, num) {memcpy(buff, int_val, num*sizeof(int)); buff += num*sizeof(int);}

#define PACK_DBL(buff, dbl_val, num) {memcpy(buff, dbl_val, num*sizeof(double)); buff += num*sizeof(double);}

#define PACK_EH(buff, eh_val, num) {memcpy(buff, eh_val, num*sizeof(MBEntityHandle)); buff += num*sizeof(MBEntityHandle);}

#define PACK_CHAR_64(buff, char_val) {strcpy((char*)buff, char_val); buff += 64;}

#define PACK_VOID(buff, val, num) {memcpy(buff, val, num); buff += num;}

#define PACK_RANGE(buff, rng) {int num_subs = num_subranges(rng); PACK_INTS(buff, &num_subs, 1); \
          for (MBRange::const_pair_iterator cit = rng.const_pair_begin(); cit != rng.const_pair_end(); cit++) { \
            MBEntityHandle eh = (*cit).first; PACK_EH(buff_ptr, &eh, 1); \
            eh = (*cit).second; PACK_EH(buff_ptr, &eh, 1);}}

#define UNPACK_INT(buff, int_val) {UNPACK_INTS(buff, &int_val, 1);}

#define UNPACK_INTS(buff, int_val, num) {memcpy(int_val, buff, num*sizeof(int)); buff += num*sizeof(int);}

#define UNPACK_DBL(buff, dbl_val, num) {memcpy(dbl_val, buff, num*sizeof(double)); buff += num*sizeof(double);}

#define UNPACK_EH(buff, eh_val, num) {memcpy(eh_val, buff, num*sizeof(MBEntityHandle)); buff += num*sizeof(MBEntityHandle);}

#define UNPACK_CHAR_64(buff, char_val) {strcpy(char_val, (char*)buff); buff += 64;}

#define UNPACK_VOID(buff, val, num) {memcpy(val, buff, num); buff += num;}

#define UNPACK_RANGE(buff, rng) {int num_subs; UNPACK_INTS(buff, &num_subs, 1); MBEntityHandle _eh[2]; \
          for (int i = 0; i < num_subs; i++) { UNPACK_EH(buff_ptr, _eh, 2); rng.insert(_eh[0], _eh[1]);}}

#define RR(a) if (MB_SUCCESS != result) {\
          dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(a);\
          return result;}

MBParallelComm::MBParallelComm(MBInterface *impl, MPI_Comm comm) 
    : mbImpl(impl), procConfig(comm)
{
  myBuffer.reserve(INITIAL_BUFF_SIZE);

  tagServer = dynamic_cast<MBCore*>(mbImpl)->tag_server();
  sequenceManager = dynamic_cast<MBCore*>(mbImpl)->sequence_manager();
}

MBParallelComm::MBParallelComm(MBInterface *impl,
                               std::vector<unsigned char> &tmp_buff, 
                               MPI_Comm comm) 
    : mbImpl(impl), procConfig(comm)
{
  myBuffer.swap(tmp_buff);
}

//! assign a global id space, for largest-dimension or all entities (and
//! in either case for vertices too)
MBErrorCode MBParallelComm::assign_global_ids(const int dimension, 
                                              const int start_id,
                                              const bool largest_dim_only) 
{
  MBRange entities[4];
  int local_num_elements[4];
  MBErrorCode result;
  for (int dim = 0; dim <= dimension; dim++) {
    if (dim == 0 || !largest_dim_only || dim == dimension) {
      result = mbImpl->get_entities_by_dimension(0, dim, entities[dim]); 
      RR("Failed to get vertices in assign_global_ids.");
    }

      // need to filter out non-locally-owned entities!!!
    MBRange dum_range;
    for (MBRange::iterator rit = entities[dim].begin(); rit != entities[dim].end(); rit++)
      if (mbImpl->handle_utils().rank_from_handle(*rit) != 
          (unsigned int) mbImpl->proc_rank()) 
        dum_range.insert(*rit);
    entities[dim] = entities[dim].subtract(dum_range);
    
    local_num_elements[dim] = entities[dim].size();
  }
  
    // communicate numbers
  std::vector<int> num_elements(procConfig.proc_size()*4);
#ifdef USE_MPI
  if (procConfig.proc_size() > 1) {
    int retval = MPI_Alltoall(local_num_elements, 4, MPI_INTEGER,
                              &num_elements[0], procConfig.proc_size()*4, 
                              MPI_INTEGER, procConfig.proc_comm());
    if (0 != retval) return MB_FAILURE;
  }
  else
#endif
    for (int dim = 0; dim < 4; dim++) num_elements[dim] = local_num_elements[dim];
  
    // my entities start at one greater than total_elems[d]
  int total_elems[4] = {start_id, start_id, start_id, start_id};
  
  for (unsigned int proc = 0; proc < procConfig.proc_rank(); proc++) {
    for (int dim = 0; dim < 4; dim++) total_elems[dim] += num_elements[4*proc + dim];
  }
  
    //.assign global ids now
  MBTag gid_tag;
  int zero = 0;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), 
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &zero, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  
  for (int dim = 0; dim < 4; dim++) {
    if (entities[dim].empty()) continue;
    num_elements.reserve(entities[dim].size());
    int i = 0;
    for (MBRange::iterator rit = entities[dim].begin(); rit != entities[dim].end(); rit++)
      num_elements[i++] = total_elems[dim]++;
    
    result = mbImpl->tag_set_data(gid_tag, entities[dim], &num_elements[0]); 
    RR("Failed to set global id tag in assign_global_ids.");
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::communicate_entities(const int from_proc, const int to_proc,
                                                 MBRange &entities,
                                                 const bool adjacencies,
                                                 const bool tags) 
{
#ifndef USE_MPI
  return MB_FAILURE;
#else
  
  MBErrorCode result = MB_SUCCESS;
  
    // if I'm the from, do the packing and sending
  if ((int)procConfig.proc_rank() == from_proc) {
    allRanges.clear();
    vertsPerEntity.clear();
    setRange.clear();
    setRanges.clear();
    allTags.clear();
    setSizes.clear();
    optionsVec.clear();
    setPcs.clear();

    MBRange whole_range;

    int buff_size;
    
    result = pack_buffer(entities, adjacencies, tags, true, 
                         whole_range, buff_size); 
    RR("Failed to compute buffer size in communicate_entities.");

      // if the message is large, send a first message to tell how large
    if (INITIAL_BUFF_SIZE < buff_size) {
      int tmp_buff_size = -buff_size;
      MPI_Request send_req;
      int success = MPI_Isend(&tmp_buff_size, sizeof(int), MPI_UNSIGNED_CHAR, to_proc, 
                              0, procConfig.proc_comm(), &send_req);
      if (!success) return MB_FAILURE;
    }
    
      // allocate space in the buffer
    myBuffer.reserve(buff_size);

      // pack the actual buffer
    int actual_buff_size;
    result = pack_buffer(entities, adjacencies, tags, false, 
                         whole_range, actual_buff_size); 
    RR("Failed to pack buffer in communicate_entities.");
    
      // send it
    MPI_Request send_req;
    int success = MPI_Isend(&myBuffer[0], actual_buff_size, MPI_UNSIGNED_CHAR, to_proc, 
                            0, procConfig.proc_comm(), &send_req);
    if (!success) return MB_FAILURE;
  }
  else if ((int)procConfig.proc_rank() == to_proc) {
    int buff_size;
    
      // get how much to allocate
    MPI_Status status;
    int success = MPI_Recv(&myBuffer[0], myBuffer.size(), MPI_UNSIGNED_CHAR, from_proc, 
                           MPI_ANY_TAG, procConfig.proc_comm(), &status);
    int num_recd;
    success = MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &num_recd);
    
    if (sizeof(int) == num_recd && 0 > *((int*)&myBuffer[0])) {
        // this was just the size of the next message; prepare buffer then receive that message
      buff_size = myBuffer[0];
      myBuffer.reserve(buff_size);
    
      // receive the real message
      success = MPI_Recv(&myBuffer[0], buff_size, MPI_UNSIGNED_CHAR, from_proc, 
                         MPI_ANY_TAG, procConfig.proc_comm(), &status);
    }
    
      // unpack the buffer
    result = unpack_buffer(entities); 
    RR("Failed to unpack buffer in communicate_entities.");
  }
  
  return result;

#endif
}

MBErrorCode MBParallelComm::broadcast_entities( const int from_proc,
                                                MBRange &entities,
                                                const bool adjacencies,
                                                const bool tags) 
{
#ifndef USE_MPI
  return MB_FAILURE;
#else
  
  MBErrorCode result = MB_SUCCESS;
  int success;
  MBRange whole_range;
  int buff_size;
  
  allRanges.clear();
  vertsPerEntity.clear();
  setRange.clear();
  setRanges.clear();
  allTags.clear();
  setSizes.clear();
  optionsVec.clear();
  setPcs.clear();

  if ((int)procConfig.proc_rank() == from_proc) {
    result = pack_buffer( entities, adjacencies, tags, true, whole_range, buff_size ); 
    RR("Failed to compute buffer size in broadcast_entities.");
  }

  success = MPI_Bcast( &buff_size, 1, MPI_INT, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success)
    return MB_FAILURE;
  
  if (!buff_size) // no data
    return MB_SUCCESS;
  
  myBuffer.reserve( buff_size );
  
  if ((int)procConfig.proc_rank() == from_proc) {
    int actual_buffer_size;
    result = pack_buffer( entities, adjacencies, tags, false, 
                          whole_range, actual_buffer_size );
    RR("Failed to pack buffer in broadcast_entities.");
  }

  success = MPI_Bcast( &myBuffer[0], buff_size, MPI_UNSIGNED_CHAR, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success)
    return MB_FAILURE;
  
  if ((int)procConfig.proc_rank() != from_proc) {
    result = unpack_buffer( entities );
    RR("Failed to unpack buffer in broadcast_entities.");
  }

  return MB_SUCCESS;
#endif
}

MBErrorCode MBParallelComm::pack_buffer(MBRange &entities, 
                                        const bool adjacencies,
                                        const bool tags,
                                        const bool just_count,
                                        MBRange &whole_range,
                                        int &buff_size) 
{
    // pack the buffer with the entity ranges, adjacencies, and tags sections
  MBErrorCode result;

  buff_size = 0;
  MBRange::const_iterator rit;
  unsigned char *buff_ptr = NULL;
  if (!just_count) buff_ptr = &myBuffer[0];
  
    // entities
  result = pack_entities(entities, rit, whole_range, buff_ptr, 
                         buff_size, just_count); 
  RR("Packing entities failed.");
  
    // sets
  int tmp_size;
  result = pack_sets(entities, rit, whole_range, buff_ptr, tmp_size, just_count); 
  RR("Packing sets failed.");
  buff_size += tmp_size;
  
    // adjacencies
  if (adjacencies) {
    result = pack_adjacencies(entities, rit, whole_range, buff_ptr, 
                              tmp_size, just_count);
    RR("Packing adjs failed.");
    buff_size += tmp_size;
  }
    
    // tags
  if (tags) {
    result = pack_tags(entities, rit, whole_range, buff_ptr, 
                       tmp_size, just_count);
    RR("Packing tags failed.");
    buff_size += tmp_size;
  }

  return result;
}
 
MBErrorCode MBParallelComm::unpack_buffer(MBRange &entities) 
{
  if (myBuffer.capacity() == 0) return MB_FAILURE;
  
  unsigned char *buff_ptr = &myBuffer[0];
  MBErrorCode result = unpack_entities(buff_ptr, entities);
  RR("Unpacking entities failed.");
  result = unpack_sets(buff_ptr, entities);
  RR("Unpacking sets failed.");
  result = unpack_tags(buff_ptr, entities);
  RR("Unpacking tags failed.");
  
  return MB_SUCCESS;
}

int MBParallelComm::num_subranges(const MBRange &this_range)
{
    // ok, have all the ranges we'll pack; count the subranges
  int num_sub_ranges = 0;
  for (MBRange::const_pair_iterator pit = this_range.const_pair_begin(); 
       pit != this_range.const_pair_end(); pit++)
    num_sub_ranges++;

  return num_sub_ranges;
}

MBErrorCode MBParallelComm::pack_entities(MBRange &entities,
                                          MBRange::const_iterator &start_rit,
                                          MBRange &whole_range,
                                          unsigned char *&buff_ptr,
                                          int &count,
                                          const bool just_count) 
{
  count = 0;
  unsigned char *orig_buff_ptr = buff_ptr;
  MBErrorCode result;
  MBWriteUtilIface *wu = NULL;
  if (!just_count) {
    result = mbImpl->query_interface(std::string("MBWriteUtilIface"), 
                                     reinterpret_cast<void**>(&wu));
    RR("Couldn't get MBWriteUtilIface.");

  }
  
    // pack vertices
  if (just_count) {
    entTypes.push_back(MBVERTEX);
    vertsPerEntity.push_back(1);
    allRanges.push_back(entities.subset_by_type(MBVERTEX));
  }
  else {
    PACK_INT(buff_ptr, MBVERTEX);
    PACK_RANGE(buff_ptr, allRanges[0]);
    int num_verts = allRanges[0].size();
    std::vector<double*> coords(3);
    for (int i = 0; i < 3; i++)
      coords[i] = reinterpret_cast<double*>(buff_ptr + i * num_verts * sizeof(double));

    assert(NULL != wu);
    
    result = wu->get_node_arrays(3, num_verts, allRanges[0], 0, 0, coords);
    RR("Couldn't allocate node space.");

    buff_ptr += 3 * num_verts * sizeof(double);

    whole_range = allRanges[0];
  }

    // place an iterator at the first non-vertex entity
  if (!allRanges[0].empty()) {
    start_rit = entities.find(*allRanges[0].rbegin());
    start_rit++;
  }
  else {
    start_rit = entities.begin();
  }
  
  MBRange::const_iterator end_rit = start_rit;
  if (allRanges[0].size() == entities.size()) return MB_SUCCESS;

  std::vector<MBRange>::iterator allr_it = allRanges.begin();
  
    // pack entities
  if (just_count) {    

      // get all ranges of entities that have different #'s of vertices or different types
    while (end_rit != entities.end() && TYPE_FROM_HANDLE(*start_rit) != MBENTITYSET) {

        // get the sequence holding this entity
      MBEntitySequence *seq;
      ElementEntitySequence *eseq;
      result = sequenceManager->find(*start_rit, seq);
      RR("Couldn't find entity sequence.");

      if (NULL == seq) return MB_FAILURE;
      eseq = dynamic_cast<ElementEntitySequence*>(seq);

        // if type and nodes per element change, start a new range
      if (eseq->get_type() != *entTypes.rbegin() || (int) eseq->nodes_per_element() != *vertsPerEntity.rbegin()) {
        entTypes.push_back(eseq->get_type());
        vertsPerEntity.push_back(eseq->nodes_per_element());
        allRanges.push_back(MBRange());
        allr_it++;
      }
    
        // get position in entities list one past end of this sequence
      end_rit = entities.lower_bound(start_rit, entities.end(), eseq->get_end_handle()+1);

        // put these entities in the last range
      eseq->get_entities(*allRanges.rbegin());
      whole_range.merge(*allRanges.rbegin());
      
        // now start where we last left off
      start_rit = end_rit;
    }

      // update vertex range and count those data, now that we know which entities get communicated
    result = mbImpl->get_adjacencies(whole_range, 0, false, allRanges[0], 
                                     MBInterface::UNION);
    RR("Failed get_adjacencies.");
    whole_range.merge(allRanges[0]);
    count += 3 * sizeof(double) * allRanges[0].size();
    
      // space for the ranges
    std::vector<MBRange>::iterator vit = allRanges.begin();
    std::vector<int>::iterator iit = vertsPerEntity.begin();
    std::vector<MBEntityType>::iterator eit = entTypes.begin();
    for (; vit != allRanges.end(); vit++, iit++, eit++) {
        // subranges of entities
      count += 2*sizeof(MBEntityHandle)*num_subranges(*vit);
        // connectivity of subrange
      if (iit != vertsPerEntity.begin()) {
        if (*eit != MBPOLYGON && *eit != MBPOLYHEDRON) 
            // for non-poly's: #verts/ent * #ents * sizeof handle
          count += *iit * (*vit).size() * sizeof(MBEntityHandle);
          // for poly's:  length of conn list * handle size + #ents * int size (for offsets)
        else count += *iit * sizeof(MBEntityHandle) + (*vit).size() * sizeof(int);
      }
    }
      //                                num_verts per subrange    ent type in subrange
    count += (vertsPerEntity.size() + 1) * (sizeof(int) + sizeof(MBEntityType));

      // extra entity type at end
    count += sizeof(int);
  }
  else {
      // for each range beyond the first
    allr_it++;
    std::vector<int>::iterator nv_it = vertsPerEntity.begin();
    std::vector<MBEntityType>::iterator et_it = entTypes.begin();
    nv_it++; et_it++;
    
    for (; allr_it != allRanges.end(); allr_it++, nv_it++, et_it++) {
        // pack the entity type
      PACK_INT(buff_ptr, *et_it);
      
        // pack the range
      PACK_RANGE(buff_ptr, (*allr_it));

        // pack the nodes per entity
      PACK_INT(buff_ptr, *nv_it);
      
        // pack the connectivity
      const MBEntityHandle *connect;
      int num_connect;
      if (*et_it == MBPOLYGON || *et_it == MBPOLYHEDRON) {
        std::vector<int> num_connects;
        for (MBRange::const_iterator rit = allr_it->begin(); rit != allr_it->end(); rit++) {
          result = mbImpl->get_connectivity(*rit, connect, num_connect);
          RR("Failed to get connectivity.");
          num_connects.push_back(num_connect);
          PACK_EH(buff_ptr, &connect[0], num_connect);
        }
        PACK_INTS(buff_ptr, &num_connects[0], num_connects.size());
      }
      else {
        for (MBRange::const_iterator rit = allr_it->begin(); rit != allr_it->end(); rit++) {
          result = mbImpl->get_connectivity(*rit, connect, num_connect);
          RR("Failed to get connectivity.");
          assert(num_connect == *nv_it);
          PACK_EH(buff_ptr, &connect[0], num_connect);
        }
      }

      whole_range.merge(*allr_it);
    }

      // pack MBMAXTYPE to indicate end of ranges
    PACK_INT(buff_ptr, MBMAXTYPE);

    count = buff_ptr - orig_buff_ptr;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_entities(unsigned char *&buff_ptr,
                                            MBRange &entities) 
{
  MBErrorCode result;
  bool done = false;
  MBReadUtilIface *ru = NULL;
  result = mbImpl->query_interface(std::string("MBReadUtilIface"), 
                                   reinterpret_cast<void**>(&ru));
  RR("Failed to get MBReadUtilIface.");

  
  while (!done) {
    MBEntityType this_type;
    UNPACK_INT(buff_ptr, this_type);
    assert(this_type >= MBVERTEX && 
           (this_type == MBMAXTYPE || this_type < MBENTITYSET));

      // MBMAXTYPE signifies end of entities data
    if (MBMAXTYPE == this_type) break;
    
      // get the range
    MBRange this_range;
    UNPACK_RANGE(buff_ptr, this_range);
    
    if (MBVERTEX == this_type) {
        // unpack coords
      int num_verts = this_range.size();
      std::vector<double*> coords(3*num_verts);
      for (MBRange::const_pair_iterator pit = this_range.const_pair_begin(); 
           pit != this_range.const_pair_end(); pit++) {
          // allocate handles
        int start_id = mbImpl->handle_utils().id_from_handle((*pit).first);
        int start_proc = mbImpl->handle_utils().rank_from_handle((*pit).first);
        MBEntityHandle actual_start;
        int tmp_num_verts = (*pit).second - (*pit).first + 1;
        result = ru->get_node_arrays(3, tmp_num_verts, start_id, start_proc, actual_start,
                                     coords);
        RR("Failed to allocate node arrays.");

        if (actual_start != (*pit).first)
          return MB_FAILURE;

        entities.insert((*pit).first, (*pit).second);
        
          // unpack the buffer data directly into coords
        for (int i = 0; i < 3; i++) 
          memcpy(coords[i], buff_ptr+i*num_verts*sizeof(double), 
                 tmp_num_verts*sizeof(double));

        buff_ptr += tmp_num_verts * sizeof(double);
      }

        // increment the buffer ptr beyond the y and z coords
      buff_ptr += 2 * num_verts * sizeof(double);
    }

    else {
      
      int verts_per_entity;
      
        // unpack the nodes per entity
      UNPACK_INT(buff_ptr, verts_per_entity);
      
        // unpack the connectivity
      for (MBRange::const_pair_iterator pit = this_range.const_pair_begin(); 
           pit != this_range.const_pair_end(); pit++) {
          // allocate handles, connect arrays
        int start_id = mbImpl->handle_utils().id_from_handle((*pit).first);
        int start_proc = mbImpl->handle_utils().rank_from_handle((*pit).first);
        MBEntityHandle actual_start;
        int num_elems = (*pit).second - (*pit).first + 1;
        MBEntityHandle *connect;
        int *connect_offsets;
        if (this_type == MBPOLYGON || this_type == MBPOLYHEDRON) {
          result = ru->get_poly_element_array(num_elems, verts_per_entity, this_type,
                                              start_id, start_proc, actual_start,
                                              connect_offsets, connect);
          RR("Failed to allocate poly element arrays.");
        }

        else {
          result = ru->get_element_array(num_elems, verts_per_entity, this_type,
                                         start_id, start_proc, actual_start,
                                         connect);
          RR("Failed to allocate element arrays.");
        }

          // copy connect arrays
        if (this_type != MBPOLYGON && this_type != MBPOLYHEDRON) {
          UNPACK_EH(buff_ptr, connect, num_elems * verts_per_entity);
        }
        else {
          UNPACK_EH(buff_ptr, connect, verts_per_entity);
          assert(NULL != connect_offsets);
            // and the offsets
          UNPACK_INTS(buff_ptr, connect_offsets, num_elems);
        }

        entities.insert((*pit).first, (*pit).second);
      }
      
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_sets(MBRange &entities,
                                      MBRange::const_iterator &start_rit,
                                      MBRange &whole_range,
                                      unsigned char *&buff_ptr,
                                      int &count,
                                      const bool just_count)
{
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  count = 0;
  unsigned char *orig_buff_ptr = buff_ptr;
  MBErrorCode result;

  if (just_count) {
    for (; start_rit != entities.end(); start_rit++) {
      setRange.insert(*start_rit);
      count += sizeof(MBEntityHandle);
    
      unsigned int options;
      result = mbImpl->get_meshset_options(*start_rit, options);
      RR("Failed to get meshset options.");
      optionsVec.push_back(options);
      count += sizeof(unsigned int);
    
      if (options & MESHSET_SET) {
          // range-based set; count the subranges
        setRanges.push_back(MBRange());
        result = mbImpl->get_entities_by_handle(*start_rit, *setRanges.rbegin());
        RR("Failed to get set entities.");
        count += 2 * sizeof(MBEntityHandle) * num_subranges(*setRanges.rbegin()) + sizeof(int);
      }
      else if (options & MESHSET_ORDERED) {
          // just get the number of entities in the set
        int num_ents;
        result = mbImpl->get_number_entities_by_handle(*start_rit, num_ents);
        RR("Failed to get number entities in ordered set.");
        
        count += sizeof(int);
        
        setSizes.push_back(num_ents);
        count += sizeof(MBEntityHandle) * num_ents + sizeof(int);
      }
      whole_range.insert(*start_rit);

        // get numbers of parents/children
      int num_par, num_ch;
      result = mbImpl->num_child_meshsets(*start_rit, &num_ch);
      RR("Failed to get num children.");

      result = mbImpl->num_parent_meshsets(*start_rit, &num_par);
      RR("Failed to get num parents.");
      count += 2*sizeof(int) + (num_par + num_ch) * sizeof(MBEntityHandle);
    
    }
  }
  else {
    
    std::vector<unsigned int>::const_iterator opt_it = optionsVec.begin();
    std::vector<MBRange>::const_iterator rit = setRanges.begin();
    std::vector<int>::const_iterator mem_it = setSizes.begin();
    static std::vector<MBEntityHandle> members;

      // set handle range
    PACK_RANGE(buff_ptr, setRange);

    for (MBRange::const_iterator set_it = setRange.begin(); set_it != setRange.end(); 
         set_it++, opt_it++) {
        // option value
      PACK_VOID(buff_ptr, &(*opt_it), sizeof(unsigned int));
      
      if ((*opt_it) & MESHSET_SET) {
          // pack entities as a range
        PACK_RANGE(buff_ptr, (*rit));
        rit++;
      }
      else if ((*opt_it) & MESHSET_ORDERED) {
          // pack entities as vector, with length
        PACK_INT(buff_ptr, *mem_it);
        members.clear();
        result = mbImpl->get_entities_by_handle(*set_it, members);
        RR("Failed to get set entities.");
        PACK_EH(buff_ptr, &members[0], *mem_it);
        mem_it++;
      }
      
        // pack parents
      members.clear();
      result = mbImpl->get_parent_meshsets(*set_it, members);
      RR("Failed to pack parents.");
      PACK_INT(buff_ptr, members.size());
      if (!members.empty()) {
        PACK_EH(buff_ptr, &members[0], members.size());
      }
      
        // pack children
      members.clear();
      result = mbImpl->get_child_meshsets(*set_it, members);
      RR("Failed to pack children.");
      PACK_INT(buff_ptr, members.size());
      if (!members.empty()) {
        PACK_EH(buff_ptr, &members[0], members.size());
      }
      
    }
    
    count = buff_ptr - orig_buff_ptr;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_sets(unsigned char *&buff_ptr,
                                        MBRange &entities)
{
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  MBErrorCode result;

  std::vector<unsigned int>::const_iterator opt_it = optionsVec.begin();
  std::vector<MBRange>::const_iterator rit = setRanges.begin();
  std::vector<int>::const_iterator mem_it = setSizes.begin();

  MBRange set_handles, new_sets;
  UNPACK_RANGE(buff_ptr, set_handles);
  std::vector<MBEntityHandle> members;
  
  for (MBRange::const_iterator rit = set_handles.begin(); 
       rit != set_handles.end(); rit++) {
    
      // option value
    unsigned int opt;
    UNPACK_VOID(buff_ptr, &opt, sizeof(unsigned int));
      
      // create the set
    MBEntityHandle set_handle = *rit;
    result = sequenceManager->allocate_mesh_set(set_handle, opt);
    RR("Failed to create set in unpack.");
    new_sets.insert(set_handle);

    int num_ents;
    if (opt & MESHSET_SET) {
        // unpack entities as a range
      MBRange set_range;
      UNPACK_RANGE(buff_ptr, set_range);
      result = mbImpl->add_entities(*rit, set_range);
      RR("Failed to add ents to set in unpack.");
    }
    else if (opt & MESHSET_ORDERED) {
        // unpack entities as vector, with length
      UNPACK_INT(buff_ptr, num_ents);
      members.reserve(num_ents);
      UNPACK_EH(buff_ptr, &members[0], num_ents);
      result = mbImpl->add_entities(*rit, &members[0], num_ents);
      RR("Failed to add ents to ordered set in unpack.");
    }
      
      // unpack parents/children
    UNPACK_INT(buff_ptr, num_ents);
    members.reserve(num_ents);
    UNPACK_EH(buff_ptr, &members[0], num_ents);
    for (int i = 0; i < num_ents; i++) {
      result = mbImpl->add_parent_meshset(*rit, members[i]);
      RR("Failed to add parent to set in unpack.");
    }
    UNPACK_INT(buff_ptr, num_ents);
    members.reserve(num_ents);
    UNPACK_EH(buff_ptr, &members[0], num_ents);
    for (int i = 0; i < num_ents; i++) {
      result = mbImpl->add_child_meshset(*rit, members[i]);
      RR("Failed to add child to set in unpack.");
    }
  }

  entities.merge(new_sets);
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_adjacencies(MBRange &entities,
                                             MBRange::const_iterator &start_rit,
                                             MBRange &whole_range,
                                             unsigned char *&buff_ptr,
                                             int &count,
                                             const bool just_count)
{
  return MB_FAILURE;
}

MBErrorCode MBParallelComm::unpack_adjacencies(unsigned char *&buff_ptr,
                                               MBRange &entities)
{
  return MB_FAILURE;
}

MBErrorCode MBParallelComm::pack_tags(MBRange &entities,
                                      MBRange::const_iterator &start_rit,
                                      MBRange &whole_range,
                                      unsigned char *&buff_ptr,
                                      int &count,
                                      const bool just_count)
{
    // tags
    // get all the tags
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  count = 0;
  unsigned char *orig_buff_ptr = buff_ptr;
  MBErrorCode result;

  if (just_count) {

    std::vector<MBTag> all_tags;
    result = tagServer->get_tags(all_tags);
    RR("Failed to get tags in pack_tags.");


    for (std::vector<MBTag>::iterator tag_it = all_tags.begin(); tag_it != all_tags.end(); tag_it++) {
      const TagInfo *tinfo = tagServer->get_tag_info(*tag_it);
      int this_count = 0;
      MBRange tmp_range;
      result = tagServer->get_entities(*tag_it, tmp_range);
      RR("Failed to get entities for tag in pack_tags.");
      tmp_range = tmp_range.intersect(whole_range);
      if (!tmp_range.empty()) this_count = tmp_range.size() * tinfo->get_size();

      if (0 == this_count) continue;

        // ok, we'll be sending this tag

        // tag handle
      allTags.push_back(*tag_it);
      count += sizeof(MBTag);
      
        // default value
      count += sizeof(int);
      if (NULL != tinfo->default_value()) count += tinfo->get_size();
      
        // size, data type
      count += sizeof(int);
      
        // data type
      count += sizeof(MBDataType);

        // name
      count += 64;

      if (!tmp_range.empty()) {
        tagRanges.push_back(tmp_range);
          // range of tag
        count += sizeof(int) + 2 * num_subranges(tmp_range) * sizeof(MBEntityHandle);
      }
      
          // tag data values for range or vector
      count += this_count;
    }

      // number of tags
    count += sizeof(int);
  }

  else {
    static std::vector<int> tag_data;
    std::vector<MBRange>::const_iterator tr_it = tagRanges.begin();

    PACK_INT(buff_ptr, allTags.size());
    
    for (std::vector<MBTag>::const_iterator tag_it = allTags.begin(); tag_it != allTags.end(); tag_it++) {

      const TagInfo *tinfo = tagServer->get_tag_info(*tag_it);

        // tag handle
      PACK_EH(buff_ptr, &(*tag_it), 1);
      
        // size, data type
      PACK_INT(buff_ptr, tinfo->get_size());
      PACK_INT(buff_ptr, tinfo->get_data_type());
      
        // default value
      if (NULL == tinfo->default_value()) {
        PACK_INT(buff_ptr, 0);
      }
      else {
        PACK_INT(buff_ptr, 1);
        PACK_VOID(buff_ptr, tinfo->default_value(), tinfo->get_size());
      }
      
        // name
      PACK_CHAR_64(buff_ptr, tinfo->get_name().c_str());
      
      tag_data.reserve((tr_it->size()+1) * tinfo->get_size() / sizeof(int));
      result = mbImpl->tag_get_data(*tag_it, *tr_it, &tag_data[0]);
      RR("Failed to get tag data in pack_tags.");
      PACK_RANGE(buff_ptr, (*tr_it));
      PACK_VOID(buff_ptr, &tag_data[0], tr_it->size()*tinfo->get_size());
      tr_it++;
    }

    count = buff_ptr - orig_buff_ptr;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_tags(unsigned char *&buff_ptr,
                                        MBRange &entities)
{
    // tags
    // get all the tags
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  MBErrorCode result;
  
  int num_tags;
  UNPACK_INT(buff_ptr, num_tags);
  std::vector<int> tag_data;

  for (int i = 0; i < num_tags; i++) {
    
        // tag handle
    MBTag tag_handle;
    UNPACK_EH(buff_ptr, &tag_handle, 1);

      // size, data type
    int tag_size, tag_data_type;
    UNPACK_INT(buff_ptr, tag_size);
    UNPACK_INT(buff_ptr, tag_data_type);
      
      // default value
    int has_def_value;
    UNPACK_INT(buff_ptr, has_def_value);
    void *def_val_ptr = NULL;
    if (1 == has_def_value) {
      def_val_ptr = buff_ptr;
      buff_ptr += tag_size;
    }
    
      // name
    char *tag_name = reinterpret_cast<char *>(buff_ptr);
    buff_ptr += 64;

      // create the tag
    MBTagType tag_type;
    result = mbImpl->tag_get_type(tag_handle, tag_type);
    RR("Failed to get tag type in unpack_tags.");
    result = mbImpl->tag_create(tag_name, tag_size, tag_type, 
                                (MBDataType) tag_data_type, tag_handle,
                                def_val_ptr);
    if (MB_ALREADY_ALLOCATED == result) {
        // already allocated tag, check to make sure it's the same size, type, etc.
      const TagInfo *tag_info = tagServer->get_tag_info(tag_name);
      if (tag_size != tag_info->get_size() ||
          tag_data_type != tag_info->get_data_type() ||
          (def_val_ptr && !tag_info->default_value() ||
           !def_val_ptr && tag_info->default_value())) {
        RR("Didn't get correct tag info when unpacking tag.");
      }
      MBTagType this_type;
      result = mbImpl->tag_get_type(tag_handle, this_type);
      if (MB_SUCCESS != result || this_type != tag_type) {
        RR("Didn't get correct tag type when unpacking tag.");
      }
    }
    else if (MB_SUCCESS != result) return result;
    
      // set the tag data; don't have to worry about dense tags with default
      // values, as those aren't sent
    MBRange tag_range;
    UNPACK_RANGE(buff_ptr, tag_range);
    result = mbImpl->tag_set_data(tag_handle, tag_range, buff_ptr);
    RR("Trouble setting range-based tag data when unpacking tag.");
    buff_ptr += tag_range.size() * tag_size;
  }
  
  return MB_SUCCESS;
}

bool MBParallelComm::buffer_size(const unsigned int new_size) 
{
  unsigned int old_size = myBuffer.size();
  myBuffer.reserve(new_size);
  return (new_size == old_size);
}

void MBParallelComm::take_buffer(std::vector<unsigned char> &new_buffer) 
{
  new_buffer.swap(myBuffer);
}

MBErrorCode MBParallelComm::resolve_shared_ents(int dim,
                                                int shared_dim) 
{
  MBErrorCode result;
  MBRange proc_ents;
  if (-1 == dim) {
    int this_dim = 3;
    while (proc_ents.empty() && this_dim >= 0) {
      result = mbImpl->get_entities_by_dimension(0, this_dim, proc_ents);
      if (MB_SUCCESS != result) return result;
      this_dim--;
    }
  }
  else {
    result = mbImpl->get_entities_by_dimension(0, dim, proc_ents);
    if (MB_SUCCESS != result) return result;
  }

  if (proc_ents.empty()) return MB_SUCCESS;
  
  return resolve_shared_ents(proc_ents, shared_dim);
}
  
MBErrorCode MBParallelComm::resolve_shared_ents(MBRange &proc_ents,
                                                int shared_dim) 
{
  MBRange::iterator rit;
  MBSkinner skinner(mbImpl);
  
    // get the skin entities by dimension
  MBRange skin_ents[4];
  MBErrorCode result;
  int upper_dim = MBCN::Dimension(TYPE_FROM_HANDLE(*proc_ents.begin()));

  int skin_dim;
  if (shared_dim < upper_dim) {
      // if shared entity dimension is less than maximal dimension,
      // start with skin entities
    skin_dim = upper_dim-1;
    result = skinner.find_skin(proc_ents, skin_ents[skin_dim],
                               skin_ents[skin_dim]);
    if (MB_SUCCESS != result) return result;
  }
  else {
      // otherwise start with original entities
    skin_ents[upper_dim] = proc_ents;
    skin_dim = upper_dim;
  }

    // get entities adjacent to skin ents from shared_dim down to
    // zero; don't create them if they don't exist already
  for (int this_dim = shared_dim; this_dim >= 0; this_dim--) {

    if (this_dim == skin_dim) continue;
      
    result = mbImpl->get_adjacencies(skin_ents[skin_dim], this_dim,
                                     false, skin_ents[this_dim],
                                     MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
  }
  
    // global id tag
  MBTag gid_tag; int def_val = -1;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &def_val, true);
  if (MB_FAILURE == result) return result;

  else if (MB_ALREADY_ALLOCATED != result) {
      // just created it, so we need global ids
    result = assign_global_ids(upper_dim);
    if (MB_SUCCESS != result) return result;
  }

    // store index in temp tag; reuse gid_data 
  std::vector<int> gid_data(skin_ents[0].size());
  int idx = 0;
  for (MBRange::iterator rit = skin_ents[0].begin(); 
       rit != skin_ents[0].end(); rit++) 
    gid_data[idx] = idx, idx++;
  MBTag idx_tag;
  result = mbImpl->tag_create("__idx_tag", sizeof(int), MB_TAG_DENSE,
                              MB_TYPE_INTEGER, idx_tag, &def_val, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  result = mbImpl->tag_set_data(idx_tag, skin_ents[0], &gid_data[0]);
  if (MB_SUCCESS != result) return result;

    // get gids for skin verts in a vector, to pass to gs
  result = mbImpl->tag_get_data(gid_tag, skin_ents[0], &gid_data[0]);
  if (MB_SUCCESS != result) return result;

    // get a crystal router
  crystal_data *cd = procConfig.crystal_router();
  
    // get total number of verts; will overshoot highest global id, but
    // that's ok
  int nverts_total, nverts_local;
  result = mbImpl->get_number_entities_by_dimension(0, 0, nverts_local);
  if (MB_SUCCESS != result) return result;
  int failure = MPI_Allreduce(&nverts_local, &nverts_total, 1,
                              MPI_INTEGER, MPI_SUM, procConfig.proc_comm());
  if (failure) return MB_FAILURE;
  
    // call gather-scatter to get shared ids & procs
  gs_data *gsd = gs_data_setup(skin_ents[0].size(), (const ulong_*)&gid_data[0], 1, cd);
  if (NULL == gsd) return MB_FAILURE;
  
    // get shared proc tags
  int def_vals[2] = {-10*procConfig.proc_size(), -10*procConfig.proc_size()};
  MBTag sharedp_tag, sharedps_tag;
  result = mbImpl->tag_create(PARALLEL_SHARED_PROC_TAG_NAME, 2*sizeof(int), 
                              MB_TAG_DENSE,
                              MB_TYPE_INTEGER, sharedp_tag, &def_vals, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  result = mbImpl->tag_create(PARALLEL_SHARED_PROCS_TAG_NAME, 
                              MAX_SHARING_PROCS*sizeof(int), 
                              MB_TAG_SPARSE,
                              MB_TYPE_INTEGER, sharedps_tag, NULL, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;

    // load shared vertices into a tuple, then sort by index
  tuple_list shared_verts;
  tuple_list_init_max(&shared_verts, 0, 2, 0, 
                      skin_ents[0].size()*MAX_SHARING_PROCS);
  int i = 0, j = 0;
  for (unsigned int p = 0; p < gsd->nlinfo->np; p++) 
    for (unsigned int np = 0; np < gsd->nlinfo->nshared[p]; np++) 
      shared_verts.vl[i++] = gsd->nlinfo->sh_ind[j++],
        shared_verts.vl[i++] = gsd->nlinfo->target[p],
        shared_verts.n++;
  std::vector<int> sort_buffer(skin_ents[0].size()*MAX_SHARING_PROCS);
  tuple_list_sort(&shared_verts, 0,(buffer*)&sort_buffer[0]);

    // set sharing procs tags on skin vertices
  int maxp = -10*procConfig.proc_size();
  int sharing_procs[MAX_SHARING_PROCS] = {maxp};
  j = 0;
  while (j < 2*shared_verts.n) {
      // count & accumulate sharing procs
    int nump = 0, this_idx = shared_verts.vl[j];
    while (shared_verts.vl[j] == this_idx)
      j++, sharing_procs[nump++] = shared_verts.vl[j++];

    sharing_procs[nump++] = procConfig.proc_rank();
    MBEntityHandle this_ent = skin_ents[0][this_idx];
    if (2 == nump)
      result = mbImpl->tag_set_data(sharedp_tag, &this_ent, 1,
                                    sharing_procs);
    else
      result = mbImpl->tag_set_data(sharedps_tag, &this_ent, 1,
                                    sharing_procs);
    if (MB_SUCCESS != result) return result;

      // reset sharing proc(s) tags
    std::fill(sharing_procs, sharing_procs+nump, maxp);
  }
  
    // set sharing procs tags on other skin ents
  const MBEntityHandle *connect; int num_connect;
  for (int d = shared_dim; d > 0; d--) {
    for (MBRange::iterator rit = skin_ents[d].begin();
         rit != skin_ents[d].end(); rit++) {
        // get connectivity
      result = mbImpl->get_connectivity(*rit, connect, num_connect);
      if (MB_SUCCESS != result) return result;
      MBRange sp_range, vp_range;
      for (int nc = 0; nc < num_connect; nc++) {
          // get sharing procs
        result = mbImpl->tag_get_data(sharedp_tag, &(*rit), 1, sharing_procs);
        if (MB_SUCCESS != result) return result;
        if (sharing_procs[0] == maxp) {
          result = mbImpl->tag_get_data(sharedps_tag, &(*rit), 1, sharing_procs);
          if (MB_SUCCESS != result) return result;
        }
          // build range of sharing procs for this vertex
        unsigned int p = 0; vp_range.clear();
        while (sharing_procs[p] != maxp && p < MAX_SHARING_PROCS)
          vp_range.insert(sharing_procs[p]), p++;
        assert(p < MAX_SHARING_PROCS);
          // intersect with range for this skin ent
        if (0 != nc) sp_range = sp_range.intersect(vp_range);
        else sp_range = vp_range;
      }
        // intersection is the owning proc(s) for this skin ent; should
        // not be empty
      assert(!sp_range.empty());
      MBRange::iterator rit2;
        // set tag for this ent
      for (j = 0, rit2 = sp_range.begin(); 
           rit2 != sp_range.end(); rit2++, j++)
        sharing_procs[j] = *rit;
      if (2 >= j)
        result = mbImpl->tag_set_data(sharedp_tag, &(*rit), 1,
                                      sharing_procs);
      else
        result = mbImpl->tag_set_data(sharedps_tag, &(*rit), 1,
                                      sharing_procs);

      if (MB_SUCCESS != result) return result;
      
        // reset sharing proc(s) tags
      std::fill(sharing_procs, sharing_procs+j, maxp);
    }
  }

    // done
  return result;
}

MBErrorCode MBParallelComm::get_shared_entities(int dim,
                                                MBRange &shared_ents) 
{
    // check shared entities
  MBTag sharedproc_tag = 0, sharedprocs_tag = 0;
  MBErrorCode result = mbImpl->tag_get_handle(PARALLEL_SHARED_PROC_TAG_NAME, 
                                              sharedproc_tag);

  result = mbImpl->tag_get_handle(PARALLEL_SHARED_PROCS_TAG_NAME, 
                                  sharedprocs_tag);

  if (0 == sharedproc_tag && 0 == sharedprocs_tag) 
    return MB_SUCCESS;

    // get the tag values
  MBEntityType start_type = MBCN::TypeDimensionMap[dim].first,
    end_type = MBCN::TypeDimensionMap[dim].second;
  std::vector<int> proc_tags;
  for (MBEntityType this_type = start_type; this_type <= end_type;
       this_type++) {
    MBRange tmp_ents;

      // PARALLEL_SHARED_PROC is a dense tag, so all ents will have a
      // value (the default value)
    if (0 != sharedproc_tag) {
      result = mbImpl->get_entities_by_type(0, this_type, tmp_ents);
      RR("Trouble getting entities for shared entities.");
      proc_tags.resize(2*tmp_ents.size());
      if (!tmp_ents.empty()) {
        result = mbImpl->tag_get_data(sharedproc_tag, 
                                      tmp_ents, &proc_tags[0]);
        RR("Trouble getting tag data for shared entities.");
      }
      int i;
      MBRange::iterator rit;
      for (i = 0, rit = tmp_ents.begin(); rit != tmp_ents.end(); i+=2, rit++) 
        if (proc_tags[i] > -1) shared_ents.insert(*rit);
    }
    if (0 != sharedprocs_tag) {
      // PARALLEL_SHARED_PROCS is a sparse tag, so only entities with this
      // tag set will have one
      result = mbImpl->get_entities_by_type_and_tag(0, this_type, 
                                                    &sharedprocs_tag,
                                                    NULL, 1, tmp_ents,
                                                    MBInterface::UNION);
      RR("Trouble getting sharedprocs_tag for shared entities.");
    }
  }

  return MB_SUCCESS;
}
  
#ifdef TEST_PARALLELCOMM

#include <iostream>

#include "MBCore.hpp"
#include "MBParallelComm.hpp"
#include "MBRange.hpp"

#define PM {std::cerr << "Test failed; error message:" << std::endl;\
          std::string errmsg; \
          dynamic_cast<MBCore*>(my_impl)->get_last_error(errmsg); \
          std::cerr << errmsg << std::endl;\
          return 1;}

int main(int argc, char* argv[])
{

    // Check command line arg
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " <mesh_file_name>" << std::endl;
    exit(1);
  }

  const char* file = argv[1];
  MBCore *my_impl = new MBCore(0, 2);
  MBInterface* mbImpl = my_impl;

    // create a communicator class, which will start mpi too
  MBParallelComm pcomm(mbImpl, my_impl->tag_server(), my_impl->sequence_manager());
  MBErrorCode result;

    // load the mesh
  result = mbImpl->load_mesh(file, 0, 0);
  if (MB_SUCCESS != result) return result;

    // get the mesh
  MBRange all_mesh, whole_range;
  result = mbImpl->get_entities_by_dimension(0, 3, all_mesh);
  if (MB_SUCCESS != result) return result;
    
  int buff_size;
  result = pcomm.pack_buffer(all_mesh, false, true, true, whole_range, buff_size);
  PM;


    // allocate space in the buffer
  pcomm.buffer_size(buff_size);

    // pack the actual buffer
  int actual_buff_size;
  result = pcomm.pack_buffer(whole_range, false, true, false, all_mesh, 
                             actual_buff_size);
  PM;

    // list the entities that got packed
  std::cout << "ENTITIES PACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);

    // get the buffer
  std::vector<unsigned char> tmp_buffer;
  pcomm.take_buffer(tmp_buffer);
    
    // stop and restart MOAB
  delete mbImpl;
  my_impl = new MBCore(1, 2);
  mbImpl = my_impl;
    
    // create a new communicator class, using our old buffer
  MBParallelComm pcomm2(mbImpl, my_impl->tag_server(), my_impl->sequence_manager(),
                        tmp_buffer);

    // unpack the results
  all_mesh.clear();
  result = pcomm2.unpack_buffer(all_mesh);
  PM;
  
  std::cout << "ENTITIES UNPACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);
  
  std::cout << "Success, processor " << mbImpl->proc_rank() << "." << std::endl;
  
  return 1;
}
#endif
