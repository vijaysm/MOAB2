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

#define RR if (MB_SUCCESS != result) return result

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
      result = mbImpl->get_entities_by_dimension(0, dim, entities[dim]); RR;
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
    
    result = mbImpl->tag_set_data(gid_tag, entities[dim], &num_elements[0]); RR;
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
    
    result = pack_buffer(entities, adjacencies, tags, true, whole_range, buff_size); RR;

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
    result = pack_buffer(entities, adjacencies, tags, false, whole_range, actual_buff_size); RR;
    
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
    result = unpack_buffer(entities); RR;
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
    result = pack_buffer( entities, adjacencies, tags, true, whole_range, buff_size ); RR;
  }

  success = MPI_Bcast( &buff_size, 1, MPI_INT, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success)
    return MB_FAILURE;
  
  if (!buff_size) // no data
    return MB_SUCCESS;
  
  myBuffer.reserve( buff_size );
  
  if ((int)procConfig.proc_rank() == from_proc) {
    int actual_buffer_size;
    result = pack_buffer( entities, adjacencies, tags, false, whole_range, actual_buffer_size ); RR;
  }

  success = MPI_Bcast( &myBuffer[0], buff_size, MPI_UNSIGNED_CHAR, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success)
    return MB_FAILURE;
  
  if ((int)procConfig.proc_rank() != from_proc) {
    result = unpack_buffer( entities ); RR;
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
  result = pack_entities(entities, rit, whole_range, buff_ptr, buff_size, just_count); RR;
  
    // sets
  int tmp_size;
  result = pack_sets(entities, rit, whole_range, buff_ptr, tmp_size, just_count); RR;
  buff_size += tmp_size;
  
    // adjacencies
  if (adjacencies) {
    result = pack_adjacencies(entities, rit, whole_range, buff_ptr, tmp_size, just_count); RR;
    buff_size += tmp_size;
  }
    
    // tags
  if (tags) {
    result = pack_tags(entities, rit, whole_range, buff_ptr, tmp_size, just_count); RR;
    buff_size += tmp_size;
  }

  return result;
}
 
MBErrorCode MBParallelComm::unpack_buffer(MBRange &entities) 
{
  if (myBuffer.capacity() == 0) return MB_FAILURE;
  
  unsigned char *buff_ptr = &myBuffer[0];
  MBErrorCode result = unpack_entities(buff_ptr, entities); RR;
  result = unpack_sets(buff_ptr, entities); RR;
  result = unpack_tags(buff_ptr, entities); RR;
  
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
    result = mbImpl->query_interface(std::string("MBWriteUtilIface"), reinterpret_cast<void**>(&wu)); RR;
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
    
    result = wu->get_node_arrays(3, num_verts, allRanges[0], 0, 0, coords); RR;

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
      result = sequenceManager->find(*start_rit, seq); RR;
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
    result = mbImpl->get_adjacencies(whole_range, 0, false, allRanges[0], MBInterface::UNION); RR;
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
          result = mbImpl->get_connectivity(*rit, connect, num_connect); RR;
          num_connects.push_back(num_connect);
          PACK_EH(buff_ptr, &connect[0], num_connect);
        }
        PACK_INTS(buff_ptr, &num_connects[0], num_connects.size());
      }
      else {
        for (MBRange::const_iterator rit = allr_it->begin(); rit != allr_it->end(); rit++) {
          result = mbImpl->get_connectivity(*rit, connect, num_connect); RR;
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
  result = mbImpl->query_interface(std::string("MBReadUtilIface"), reinterpret_cast<void**>(&ru)); RR;
  
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
                                     coords); RR;
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
        if (this_type == MBPOLYGON || this_type == MBPOLYHEDRON)
          result = ru->get_poly_element_array(num_elems, verts_per_entity, this_type,
                                              start_id, start_proc, actual_start,
                                              connect_offsets, connect); RR;
        else
          result = ru->get_element_array(num_elems, verts_per_entity, this_type,
                                         start_id, start_proc, actual_start,
                                         connect); RR;

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
      result = mbImpl->get_meshset_options(*start_rit, options); RR;
      optionsVec.push_back(options);
      count += sizeof(unsigned int);
    
      if (options & MESHSET_SET) {
          // range-based set; count the subranges
        setRanges.push_back(MBRange());
        result = mbImpl->get_entities_by_handle(*start_rit, *setRanges.rbegin()); RR;
        count += 2 * sizeof(MBEntityHandle) * num_subranges(*setRanges.rbegin()) + sizeof(int);
      }
      else if (options & MESHSET_ORDERED) {
          // just get the number of entities in the set
        int num_ents;
        result = mbImpl->get_number_entities_by_handle(*start_rit, num_ents); RR;
        count += sizeof(int);
        
        setSizes.push_back(num_ents);
        count += sizeof(MBEntityHandle) * num_ents + sizeof(int);
      }
      whole_range.insert(*start_rit);

        // get numbers of parents/children
      int num_par, num_ch;
      result = mbImpl->num_child_meshsets(*start_rit, &num_ch); RR;
      result = mbImpl->num_parent_meshsets(*start_rit, &num_par); RR;
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
        result = mbImpl->get_entities_by_handle(*set_it, members); RR;
        PACK_EH(buff_ptr, &members[0], *mem_it);
        mem_it++;
      }
      
        // pack parents
      members.clear();
      result = mbImpl->get_parent_meshsets(*set_it, members); RR;
      PACK_INT(buff_ptr, members.size());
      if (!members.empty()) {
        PACK_EH(buff_ptr, &members[0], members.size());
      }
      
        // pack children
      members.clear();
      result = mbImpl->get_child_meshsets(*set_it, members); RR;
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

  MBRange set_handles;
  UNPACK_RANGE(buff_ptr, set_handles);
  std::vector<MBEntityHandle> members;
  
  for (MBRange::const_iterator rit = set_handles.begin(); rit != set_handles.end(); rit++) {
    
      // option value
    unsigned int opt;
    UNPACK_VOID(buff_ptr, &opt, sizeof(unsigned int));
      
      // create the set
    MBEntityHandle set_handle;
    result = mbImpl->create_meshset(opt, set_handle, 
                                    mbImpl->handle_utils().id_from_handle(*rit), 
                                    mbImpl->handle_utils().rank_from_handle(*rit)); RR;
    if (set_handle != *rit)
      return MB_FAILURE;

    int num_ents;
    if (opt & MESHSET_SET) {
        // unpack entities as a range
      MBRange set_range;
      UNPACK_RANGE(buff_ptr, set_range);
      result = mbImpl->add_entities(*rit, set_range); RR;
    }
    else if (opt & MESHSET_ORDERED) {
        // unpack entities as vector, with length
      UNPACK_INT(buff_ptr, num_ents);
      members.reserve(num_ents);
      UNPACK_EH(buff_ptr, &members[0], num_ents);
      result = mbImpl->add_entities(*rit, &members[0], num_ents); RR;
    }
      
      // unpack parents/children
    UNPACK_INT(buff_ptr, num_ents);
    members.reserve(num_ents);
    UNPACK_EH(buff_ptr, &members[0], num_ents);
    for (int i = 0; i < num_ents; i++) {
      result = mbImpl->add_parent_meshset(*rit, members[i]); RR;
    }
    UNPACK_INT(buff_ptr, num_ents);
    members.reserve(num_ents);
    UNPACK_EH(buff_ptr, &members[0], num_ents);
    for (int i = 0; i < num_ents; i++) {
      result = mbImpl->add_child_meshset(*rit, members[i]); RR;
    }
  }
  
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
  int whole_size = whole_range.size();

  if (just_count) {

    std::vector<MBTag> all_tags;
    result = tagServer->get_tags(all_tags); RR;

    for (std::vector<MBTag>::iterator tag_it = all_tags.begin(); tag_it != all_tags.end(); tag_it++) {
      const TagInfo *tinfo = tagServer->get_tag_info(*tag_it);
      int this_count = 0;
      MBRange tmp_range;
      if (PROP_FROM_TAG_HANDLE(*tag_it) == MB_TAG_DENSE) {
        this_count += whole_size * tinfo->get_size();
      }
      else {
        result = tagServer->get_entities(*tag_it, MBMAXTYPE, tmp_range); RR;
        tmp_range = tmp_range.intersect(whole_range);
        if (!tmp_range.empty()) this_count = tmp_range.size() * tinfo->get_size();
      }

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
      
      if (PROP_FROM_TAG_HANDLE(*tag_it) == MB_TAG_DENSE) {
        tag_data.reserve((whole_size+1) * tinfo->get_size() / sizeof(int));
        result = mbImpl->tag_get_data(*tag_it, whole_range, &tag_data[0]);
        PACK_VOID(buff_ptr, &tag_data[0], whole_size*tinfo->get_size());
      }
      else {
        tag_data.reserve((tr_it->size()+1) * tinfo->get_size() / sizeof(int));
        result = mbImpl->tag_get_data(*tag_it, *tr_it, &tag_data[0]); RR;
        PACK_RANGE(buff_ptr, (*tr_it));
        PACK_VOID(buff_ptr, &tag_data[0], tr_it->size()*tinfo->get_size());
        tr_it++;
      }
      
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
    result = mbImpl->tag_get_type(tag_handle, tag_type); RR;

    result = mbImpl->tag_create(tag_name, tag_size, tag_type, (MBDataType) tag_data_type, tag_handle,
                                def_val_ptr);
    if (MB_ALREADY_ALLOCATED == result) {
        // already allocated tag, check to make sure it's the same size, type, etc.
      const TagInfo *tag_info = tagServer->get_tag_info(tag_name);
      if (tag_size != tag_info->get_size() ||
          tag_data_type != tag_info->get_data_type() ||
          (def_val_ptr && !tag_info->default_value() ||
           !def_val_ptr && tag_info->default_value()))
        return MB_FAILURE;
      MBTagType this_type;
      result = mbImpl->tag_get_type(tag_handle, this_type);
      if (MB_SUCCESS != result || this_type != tag_type) return MB_FAILURE;
    }
    else if (MB_SUCCESS != result) return result;
    
      // set the tag data
    if (PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_DENSE) {
      if (NULL != def_val_ptr && tag_data_type != MB_TYPE_OPAQUE) {
          // only set the tags whose values aren't the default value; only works
          // if it's a known type
        MBRange::iterator start_rit = entities.begin(), end_rit = start_rit;
        MBRange set_ents;
        while (end_rit != entities.end()) {
          while (start_rit != entities.end() &&
                 ((tag_data_type == MB_TYPE_INTEGER && *((int*)def_val_ptr) == *((int*)buff_ptr)) ||
                  (tag_data_type == MB_TYPE_DOUBLE && *((double*)def_val_ptr) == *((double*)buff_ptr)) ||
                  (tag_data_type == MB_TYPE_HANDLE && *((MBEntityHandle*)def_val_ptr) == *((MBEntityHandle*)buff_ptr)))) {
            start_rit++;
            buff_ptr += tag_size;
          }
          end_rit = start_rit;
          void *end_ptr = buff_ptr;
          while (start_rit != entities.end() && end_rit != entities.end() &&
                 ((tag_data_type == MB_TYPE_INTEGER && *((int*)def_val_ptr) == *((int*)end_ptr)) ||
                  (tag_data_type == MB_TYPE_DOUBLE && *((double*)def_val_ptr) == *((double*)end_ptr)) ||
                  (tag_data_type == MB_TYPE_HANDLE && *((MBEntityHandle*)def_val_ptr) == *((MBEntityHandle*)end_ptr)))) {
            set_ents.insert(*end_rit);
            end_rit++;
            buff_ptr += tag_size;
          }
          
          if (!set_ents.empty()) {
            result = mbImpl->tag_set_data(tag_handle, set_ents, buff_ptr); RR;
          }
          if (start_rit != entities.end()) {
            end_rit++;
            start_rit = end_rit;
            buff_ptr += tag_size;
          }
        }
      }
      else {
        result = mbImpl->tag_set_data(tag_handle, entities, buff_ptr); RR;
        buff_ptr += entities.size() * tag_size;
      }
    }
    else {
      MBRange tag_range;
      UNPACK_RANGE(buff_ptr, tag_range);
      result = mbImpl->tag_set_data(tag_handle, tag_range, buff_ptr); RR;
      buff_ptr += tag_range.size() * tag_size;
    }
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

MBErrorCode MBParallelComm::resolve_shared_ents(MBRange &proc_ents,
                                                const int dim) 
{
  MBRange::iterator rit;
  MBSkinner skinner(mbImpl);
  
    // get the skin entities by dimension
  MBRange skin_ents[3];
  MBErrorCode result;
  int upper_dim = MBCN::Dimension(TYPE_FROM_HANDLE(*proc_ents.begin()));

  if (upper_dim > 0) {
      // first get d-1 skin ents
    result = skinner.find_skin(proc_ents, skin_ents[upper_dim-1],
                               skin_ents[upper_dim-1]);
    if (MB_SUCCESS != result) return result;
      // then get d-2, d-3, etc. entities adjacent to skin ents 
    for (int this_dim = upper_dim-1; this_dim >= 0; this_dim--) {
      result = mbImpl->get_adjacencies(skin_ents[upper_dim-1], this_dim,
                                       true, skin_ents[this_dim]);
      if (MB_SUCCESS != result) return result;
    }
  }
  else skin_ents[0] = proc_ents;
  
    // global id tag
  MBTag gid_tag; int def_val = -1;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &def_val, true);
  if (MB_FAILURE == result) return result;
  else if (MB_ALREADY_ALLOCATED != result) {
      // just created it, so we need global ids
    result = assign_global_ids(dim);
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
#define MAX_SHARING_PROCS 10  
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
  for (int d = dim-1; d > 0; d--) {
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


  
  
#ifdef TEST_PARALLELCOMM

#include <iostream>

#include "MBCore.hpp"
#include "MBParallelComm.hpp"
#include "MBRange.hpp"

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
  result = pcomm.pack_buffer(all_mesh, false, true, true, whole_range, buff_size); RR;

    // allocate space in the buffer
  pcomm.buffer_size(buff_size);

    // pack the actual buffer
  int actual_buff_size;
  result = pcomm.pack_buffer(whole_range, false, true, false, all_mesh, actual_buff_size); RR;

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
  result = pcomm2.unpack_buffer(all_mesh); RR;
  std::cout << "ENTITIES UNPACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);
  
  std::cout << "Success, processor " << mbImpl->proc_rank() << "." << std::endl;
  
  return 1;
}
#endif
