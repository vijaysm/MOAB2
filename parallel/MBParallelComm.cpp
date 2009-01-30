#include "MBInterface.hpp"
#include "MBParallelComm.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBReadUtilIface.hpp"
#include "SequenceManager.hpp"
#include "EntitySequence.hpp"
#include "TagServer.hpp"
#include "MBTagConventions.hpp"
#include "MBSkinner.hpp"
#include "MBParallelConventions.h"
#include "MBCore.hpp"
#include "MBError.hpp"
#include "ElementSequence.hpp"
#include "MBCN.hpp"
#include "RangeMap.hpp"
#include "MeshTopoUtil.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>

#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
const bool debug = false;

#include <math.h>
#include <assert.h>


extern "C" 
{
#include "minmax.h"
#include "gs.h"
#include "errmem.h"
#include "types.h"
#include "sort.h"
#include "tuple_list.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

#define INITIAL_BUFF_SIZE 1024

#undef DEBUG_PACKING
#ifdef DEBUG_PACKING
unsigned int __PACK_num = 0, __UNPACK_num = 0, __PACK_count = 0, __UNPACK_count = 0;
std::string __PACK_string, __UNPACK_string;

#define PC(n, m) {\
          if (__PACK_num == (unsigned int)n && __PACK_string == m) __PACK_count++;\
          else {\
            if (__PACK_count > 1) std::cerr << " (" << __PACK_count << "x)";\
            __PACK_count = 1; __PACK_string = m; __PACK_num = n;\
            std::cerr << std::endl << "PACK: " << n << m;\
          }}
#define UPC(n, m) {\
          if (__UNPACK_num == (unsigned int)n && __UNPACK_string == m) __UNPACK_count++;\
          else {\
            if (__UNPACK_count > 1) std::cerr << "(" << __UNPACK_count << "x)";\
            __UNPACK_count = 1; __UNPACK_string = m; __UNPACK_num = n;\
            std::cerr << std::endl << "UNPACK: " << n << m;\
          }}
#else
#define PC(n, m)
#define UPC(n, m)
#endif

#define PACK_INT(buff, int_val) {int tmp_val = int_val; PACK_INTS(buff, &tmp_val, 1);}

#define PACK_INTS(buff, int_val, num) {memcpy(buff, int_val, (num)*sizeof(int)); buff += (num)*sizeof(int); PC(num, " ints");}

#define PACK_DBL(buff, dbl_val, num) {memcpy(buff, dbl_val, (num)*sizeof(double)); buff += (num)*sizeof(double); PC(num, " doubles");}

#define PACK_EH(buff, eh_val, num) {memcpy(buff, eh_val, (num)*sizeof(MBEntityHandle)); buff += (num)*sizeof(MBEntityHandle); PC(num, " handles");}

#define PACK_CHAR_64(buff, char_val) {strcpy((char*)buff, char_val); buff += 64; PC(64, " chars");}

#define PACK_VOID(buff, val, num) {memcpy(buff, val, num); buff += num; PC(num, " void");}

#define PACK_BYTES(buff, val, num) PACK_INT(buff, num) PACK_VOID(buff, val, num)

#define PACK_RANGE(buff, rng) {int num_subs = num_subranges(rng); PACK_INTS(buff, &num_subs, 1); PC(num_subs, "-subranged range"); \
          for (MBRange::const_pair_iterator cit = rng.const_pair_begin(); cit != rng.const_pair_end(); cit++) { \
            MBEntityHandle eh = (*cit).first; PACK_EH(buff, &eh, 1); \
            eh = (*cit).second; PACK_EH(buff, &eh, 1);}; }

#define UNPACK_INT(buff, int_val) {UNPACK_INTS(buff, &int_val, 1);}

#define UNPACK_INTS(buff, int_val, num) {memcpy(int_val, buff, (num)*sizeof(int)); buff += (num)*sizeof(int); UPC(num, " ints");}

#define UNPACK_DBL(buff, dbl_val, num) {memcpy(dbl_val, buff, (num)*sizeof(double)); buff += (num)*sizeof(double); UPC(num, " doubles");}

#define UNPACK_EH(buff, eh_val, num) {memcpy(eh_val, buff, (num)*sizeof(MBEntityHandle)); buff += (num)*sizeof(MBEntityHandle); UPC(num, " handles");}

#define UNPACK_CHAR_64(buff, char_val) {strcpy(char_val, (char*)buff); buff += 64; UPC(64, " chars");}

#define UNPACK_VOID(buff, val, num) {memcpy(val, buff, num); buff += num; UPC(num, " void");}

#define UNPACK_RANGE(buff, rng) {int num_subs; UNPACK_INTS(buff, &num_subs, 1); UPC(num_subs, "-subranged range"); MBEntityHandle _eh[2]; \
          for (int i = 0; i < num_subs; i++) { UNPACK_EH(buff, _eh, 2); rng.insert(_eh[0], _eh[1]);}}
#define CHECK_BUFF_SPACE(buff_vec, buff_ptr, addl_space) { \
      unsigned int _new_size = buff_ptr - &buff_vec[0] + (addl_space);  \
    if (_new_size > buff_vec.capacity()) {  \
      buff_vec.reserve(1.5*_new_size); \
      buff_vec.resize(_new_size - (addl_space));            \
      buff_ptr = &buff_vec[_new_size-(addl_space)];}}

#define RANGE_SIZE(rng) (2*sizeof(MBEntityHandle)*num_subranges(rng)+sizeof(int))
#define RR(a) if (MB_SUCCESS != result) {\
          dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(a);\
          return result;}

#define RRA(a) if (MB_SUCCESS != result) {\
      std::string tmp_str; mbImpl->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(tmp_str.c_str()); \
      return result;}

/** Name of tag used to store MBParallelComm Index on mesh paritioning sets */
const char* PARTITIONING_PCOMM_TAG_NAME = "__PRTN_PCOMM";
 
/** \brief Tag storing parallel communication objects
 *
 * This tag stores pointers to MBParallelComm communication
 * objects; one of these is allocated for each different
 * communicator used to read mesh.  MBParallelComm stores
 * partition and interface sets corresponding to its parallel mesh.
 * By default, a parallel read uses the first MBParallelComm object
 * on the interface instance; if instantiated with one, ReadParallel
 * adds this object to the interface instance too.
 *
 * Tag type: opaque
 * Tag size: MAX_SHARING_PROCS*sizeof(MBParallelComm*)
 */
#define PARALLEL_COMM_TAG_NAME "__PARALLEL_COMM"


enum MBMessageTag {MB_MESG_ANY=MPI_ANY_TAG, 
                   MB_MESG_SIZE,
                   MB_MESG_ENTS,
                   MB_MESG_REMOTE_HANDLES,
                   MB_MESG_SHAREDHPS,
                   MB_MESG_TAGS };
    
MBParallelComm::MBParallelComm(MBInterface *impl, MPI_Comm comm, int* id ) 
        : mbImpl(impl), procConfig(comm),
          sharedpTag(0), sharedpsTag(0),
          sharedhTag(0), sharedhsTag(0), pstatusTag(0), ifaceSetsTag(0),
          partitionTag(0), globalPartCount(-1), partitioningSet(0)
{
  myBuffer.resize(INITIAL_BUFF_SIZE);

  tagServer = dynamic_cast<MBCore*>(mbImpl)->tag_server();
  sequenceManager = dynamic_cast<MBCore*>(mbImpl)->sequence_manager();

  int flag = 1;
  int retval = MPI_Initialized(&flag);
  if (MPI_SUCCESS != retval || !flag) {
    int argc = 0;
    char **argv = NULL;
    
      // mpi not initialized yet - initialize here
    retval = MPI_Init(&argc, &argv);
  }

  pcommID = add_pcomm(this);
  if (id)
    *id = pcommID;
}

MBParallelComm::MBParallelComm(MBInterface *impl,
                               std::vector<unsigned char> &tmp_buff, 
                               MPI_Comm comm,
                               int* id) 
    : mbImpl(impl), procConfig(comm),
      sharedpTag(0), sharedpsTag(0),
      sharedhTag(0), sharedhsTag(0), pstatusTag(0), ifaceSetsTag(0),
      partitionTag(0), globalPartCount(-1), partitioningSet(0)
{
  myBuffer.swap(tmp_buff);
  int flag = 1;
  int retval = MPI_Initialized(&flag);
  if (MPI_SUCCESS != retval || !flag) {
    int argc = 0;
    char **argv = NULL;
    
      // mpi not initialized yet - initialize here
    retval = MPI_Init(&argc, &argv);
  }

  pcommID = add_pcomm(this);
  if (id)
    *id = pcommID;
}

MBParallelComm::~MBParallelComm() 
{
  remove_pcomm(this);
}

int MBParallelComm::add_pcomm(MBParallelComm *pc) 
{
    // add this pcomm to instance tag
  std::vector<MBParallelComm *> pc_array(MAX_SHARING_PROCS, 
                                         (MBParallelComm*)NULL);
  MBTag pc_tag = pcomm_tag(mbImpl, true);
  assert(0 != pc_tag);
  
  MBErrorCode result = mbImpl->tag_get_data(pc_tag, 0, 0, (void*)&pc_array[0]);
  if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) 
    return -1;
  int index = 0;
  while (index < MAX_SHARING_PROCS && pc_array[index]) index++;
  if (index == MAX_SHARING_PROCS) {
    index = -1;
    assert(false);
  }
  else {
    pc_array[index] = pc;
    mbImpl->tag_set_data(pc_tag, 0, 0, (void*)&pc_array[0]);
  }
  return index;
}

void MBParallelComm::remove_pcomm(MBParallelComm *pc) 
{
    // remove this pcomm from instance tag
  std::vector<MBParallelComm *> pc_array(MAX_SHARING_PROCS);
  MBTag pc_tag = pcomm_tag(mbImpl, true);
  
  MBErrorCode result = mbImpl->tag_get_data(pc_tag, 0, 0, (void*)&pc_array[0]);
  std::vector<MBParallelComm*>::iterator pc_it = 
    std::find(pc_array.begin(), pc_array.end(), pc);
  assert(MB_SUCCESS == result && 
         pc_it != pc_array.end());
  *pc_it = NULL;
  mbImpl->tag_set_data(pc_tag, 0, 0, (void*)&pc_array[0]);
}

//! assign a global id space, for largest-dimension or all entities (and
//! in either case for vertices too)
MBErrorCode MBParallelComm::assign_global_ids(MBEntityHandle this_set,
                                              const int dimension, 
                                              const int start_id,
                                              const bool largest_dim_only,
                                              const bool parallel) 
{
  MBRange entities[4];
  int local_num_elements[4];
  MBErrorCode result;
  std::vector<unsigned char> pstatus;
  for (int dim = 0; dim <= dimension; dim++) {
    if (dim == 0 || !largest_dim_only || dim == dimension) {
      result = mbImpl->get_entities_by_dimension(this_set, dim, entities[dim]); 
      RRA("Failed to get vertices in assign_global_ids.");
    }

      // need to filter out non-locally-owned entities!!!
    pstatus.resize(entities[dim].size());
    result = mbImpl->tag_get_data(pstatus_tag(), entities[dim], &pstatus[0]);
    RRA("Failed to get pstatus in assign_global_ids.");
    
    MBRange dum_range;
    MBRange::iterator rit;
    unsigned int i;
    for (rit = entities[dim].begin(), i = 0; rit != entities[dim].end(); rit++, i++)
      if (pstatus[i] & PSTATUS_NOT_OWNED)
        dum_range.insert(*rit);
    entities[dim] = entities[dim].subtract(dum_range);
    
    local_num_elements[dim] = entities[dim].size();
  }
  
    // communicate numbers
  std::vector<int> num_elements(procConfig.proc_size()*4);
#ifdef USE_MPI
  if (procConfig.proc_size() > 1 && parallel) {
    int retval = MPI_Alltoall(local_num_elements, 4, MPI_INT,
                              &num_elements[0], procConfig.proc_size()*4, 
                              MPI_INT, procConfig.proc_comm());
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
  
    //assign global ids now
  MBTag gid_tag;
  int zero = 0;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), 
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &zero, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  
  for (int dim = 0; dim < 4; dim++) {
    if (entities[dim].empty()) continue;
    num_elements.resize(entities[dim].size());
    int i = 0;
    for (MBRange::iterator rit = entities[dim].begin(); rit != entities[dim].end(); rit++)
      num_elements[i++] = total_elems[dim]++;
    
    result = mbImpl->tag_set_data(gid_tag, entities[dim], &num_elements[0]); 
    RRA("Failed to set global id tag in assign_global_ids.");
  }
  
  return MB_SUCCESS;
}

/*
MBErrorCode MBParallelComm::send_entities(const int to_proc,
                                          MBRange &orig_ents,
                                          const bool adjs,
                                          const bool tags,
                                          const bool store_remote_handles,
                                          MBRange &final_ents,
                                          bool wait_all) 
{

    // find a spot in the proc buffers vector
  int index = get_buffers(to_proc);
  
    // pack/send the entities to the destination
  MPI_Request req_vec[2] = {MPI_REQUEST_NULL,  // send req (ownerSBuff)
                             MPI_REQUEST_NULL}; // recv req (remote handles)
  MBErrorCode result = pack_send_entities(to_proc, orig_ents, 
                                          adjs, tags, store_remote_handles, false,
                                          ownerSBuffs[index], ownerRBuffs[index],
                                          req_vec[0], req_vec[1], final_ents);
  RRA("Failed to pack-send entities.");

  MPI_Status statuses[2];

    // if we're storing remote handles, process the return message
  if (store_remote_handles) {
    assert(MPI_REQUEST_NULL != req_vec[1]);
    MPI_Wait(&req_vec[1], &statuses[1]);
    MBRange remote_range;
    unsigned char *buff_ptr = &ownerRBuffs[index][0];
    UNPACK_RANGE(buff_ptr, remote_range);
    result = set_remote_data(final_ents, remote_range, to_proc);
    RRA("Trouble setting remote data range on sent entities.");
  }

  if (wait_all) {
      // wait for all to finish
    int success = MPI_Waitall(2, req_vec, statuses);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed in waitall in send_entities.");
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::recv_entities(const int from_proc,
                                          const bool store_remote_handles,
                                          MBRange &final_ents,
                                          bool wait_all) 
{
  MBErrorCode result;
  
    // find a spot in the proc buffers vector
  int index = get_buffers(from_proc);
  
    // recv/unpack entities from the sender
  MPI_Request req = MPI_REQUEST_NULL;
  MPI_Status status;
  int success = MPI_Recv(&ghostRBuffs[index][0], ghostRBuffs[index].size(), 
                         MPI_UNSIGNED_CHAR, from_proc, 
                         MB_MESG_ANY, procConfig.proc_comm(), 
                         &status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Recv failed in recv_ents.");
  }

    // check type, and if it's just the size, resize buffer and re-call recv
  if (MB_MESG_SIZE == status.MPI_TAG) {
    result = recv_size_buff(from_proc,
                            ghostRBuffs[index],
                            req, MB_MESG_ANY);
    RRA("Failed to resize recv buffer.");
    MPI_Wait(&req, &status);
  }
  
    // ok, got the actual entities, now unpack
  result = recv_unpack_entities(from_proc, store_remote_handles, false,
                                ghostRBuffs[index], ghostSBuffs[index], 
                                req, final_ents);
  RRA("Failed to recv-unpack message.");

  if (wait_all) {
      // wait for last message to finish
    int success = MPI_Wait(&req, &status);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed in wait in recv_entities.");
    }
  }

  return result;
}
*/

int MBParallelComm::get_buffers(int to_proc, bool *is_new) 
{
  int ind = -1;
  std::vector<unsigned int>::iterator vit = 
    std::find(buffProcs.begin(), buffProcs.end(), to_proc);
  if (vit == buffProcs.end()) {
    ind = buffProcs.size();
    buffProcs.push_back((unsigned int)to_proc);
    ownerSBuffs.push_back(std::vector<unsigned char>(INITIAL_BUFF_SIZE));
    ghostRBuffs.push_back(std::vector<unsigned char>(INITIAL_BUFF_SIZE));
      // allocate these other buffs in case we're storing remote handles
    ownerRBuffs.push_back(std::vector<unsigned char>(INITIAL_BUFF_SIZE));
    ghostSBuffs.push_back(std::vector<unsigned char>(INITIAL_BUFF_SIZE));
    sharedEnts.push_back(GhostStruct());
    if (is_new) *is_new = true;
  }
  else {
    ind = vit - buffProcs.begin();
    if (is_new) *is_new = false;
  }
  assert(ind < MAX_SHARING_PROCS);
  return ind;
}

MBErrorCode MBParallelComm::pack_send_entities(const int to_proc,
                                               MBRange &orig_ents,
                                               std::set<int> &addl_procs,
                                               const bool store_remote_handles,
                                               const bool iface_layer,
                                               std::vector<unsigned char> &send_buff,
                                               std::vector<unsigned char> &recv_buff,
                                               MPI_Request &send_req,
                                               MPI_Request &recv_req) 
{
#ifndef USE_MPI
  return MB_FAILURE;
#else

  MBErrorCode result = MB_SUCCESS;
  MBRange whole_range;
  
/*
    // get an estimate of the buffer size
  int buff_size = estimate_ents_buffer_size(orig_ents, store_remote_handles);
  send_buff.clear();
  send_buff.reserve(buff_size);

  result = pack_buffer(orig_ents, addl_procs, 
                       store_remote_handles, iface_layer, to_proc,
                       send_buff, buff_size); 
  RRA("Failed to pack buffer in pack_send.");

    // now that we know how many entities, post a receive for the 
    // remote range, if requested
  int success;
  if (store_remote_handles && iface_layer) {
    recv_buff.resize(final_ents.size()*sizeof(MBEntityHandle) + sizeof(int));
    
    success = MPI_Irecv(&recv_buff[0], recv_buff.size(), MPI_UNSIGNED_CHAR, to_proc, 
                        MB_MESG_REMOTE_HANDLES_VECTOR, procConfig.proc_comm(), 
                        &recv_req);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
  else if (store_remote_handles && !iface_layer) {
    recv_buff.resize(2*num_subranges(final_ents)*sizeof(MBEntityHandle) + sizeof(int));
    
    success = MPI_Irecv(&recv_buff[0], recv_buff.size(), MPI_UNSIGNED_CHAR, to_proc, 
                        MB_MESG_REMOTE_HANDLES_RANGE, procConfig.proc_comm(), 
                        &recv_req);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }


    // if the message is large, send a first message to tell how large
  int success;
  if (INITIAL_BUFF_SIZE < buff_size) {
    int tmp_buff_size = -buff_size;
    success = MPI_Isend(&tmp_buff_size, sizeof(int), MPI_UNSIGNED_CHAR, 
                            to_proc, MB_MESG_SIZE, procConfig.proc_comm(),
                            &send_req);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
    
    // send the buffer
  success = MPI_Isend(&send_buff[0], buff_size, MPI_UNSIGNED_CHAR, to_proc, 
                      MB_MESG_ENTS, procConfig.proc_comm(), &send_req);
  if (success != MPI_SUCCESS) return MB_FAILURE;
*/

  return result;
#endif
}

MBErrorCode MBParallelComm::send_buffer(const unsigned int to_proc,
                                        const unsigned char *send_buff,
                                        const unsigned int buff_size,
                                        const int msg_type,
                                        MPI_Request &send_req) 
{
#ifndef USE_MPI
  return MB_FAILURE;
#else

  MBErrorCode result = MB_SUCCESS;
  int success;

    // if the message is large, send a first message to tell how large
  if (INITIAL_BUFF_SIZE < buff_size) {
    int tmp_buff_size = -buff_size;
    int success = MPI_Isend(&tmp_buff_size, sizeof(int), MPI_UNSIGNED_CHAR, 
                            to_proc, MB_MESG_SIZE, procConfig.proc_comm(),
                            &send_req);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
    
    // send the buffer
  success = MPI_Isend(const_cast<unsigned char*>(send_buff), buff_size, MPI_UNSIGNED_CHAR, to_proc, 
                      msg_type, procConfig.proc_comm(), &send_req);
  if (success != MPI_SUCCESS) return MB_FAILURE;

  return result;
#endif
}

MBErrorCode MBParallelComm::recv_size_buff(const int from_proc,
                                           std::vector<unsigned char> &recv_buff,
                                           MPI_Request &recv_req,
                                           int mesg_tag) 
{
    // use the received size to resize buffer, then post another irecv
  recv_buff.resize(-(*((int*)&recv_buff[0])));
  int success = MPI_Irecv(&recv_buff[0], recv_buff.size(), MPI_UNSIGNED_CHAR, from_proc, 
                          mesg_tag, procConfig.proc_comm(), &recv_req);
  if (MPI_SUCCESS != success) {
    MBErrorCode result = MB_FAILURE;
    RRA("Failed call to Irecv in recv_size_buff.");
  }
  
  return MB_SUCCESS;
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
  int buff_size;
  
  std::vector<unsigned char> buff;
  std::vector<int> addl_procs;
  if ((int)procConfig.proc_rank() == from_proc) {
    MBRange tmp_ents;
    result = add_verts(entities);
    RRA("Failed to add adj vertices.");
    entities.swap(tmp_ents);

    result = pack_buffer( entities, addl_procs, adjacencies, tags, 
                          false, -1, tmp_ents, buff, buff_size); 
    RRA("Failed to compute buffer size in broadcast_entities.");
  }

  success = MPI_Bcast( &buff_size, 1, MPI_INT, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("MPI_Bcast of buffer size failed.");
  }
  
  if (!buff_size) // no data
    return MB_SUCCESS;

  if ((int)procConfig.proc_rank() != from_proc) 
    buff.resize(buff_size);

  success = MPI_Bcast( &buff[0], buff_size, MPI_UNSIGNED_CHAR, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("MPI_Bcast of buffer failed.");
  }

  if ((int)procConfig.proc_rank() != from_proc) {
    std::set<unsigned int> my_addl_procs;
    result = unpack_buffer(&buff[0], false, from_proc, -1, entities,
                           NULL, my_addl_procs);
    RRA("Failed to unpack buffer in broadcast_entities.");
  }

  return MB_SUCCESS;
#endif
}

MBErrorCode MBParallelComm::pack_buffer(MBRange &orig_ents, 
                                        std::vector<int> &addl_procs,
                                        const bool adjacencies,
                                        const bool tags,
                                        const bool store_remote_handles,
                                        const int to_proc,
                                        MBRange &owned_shared,
                                        std::vector<unsigned char> &buff,
                                        int &buff_size) 
{
    // pack the buffer with the entity ranges, adjacencies, and tags sections
    // if do_filter is true, remove already shared entities and add adj vertices
    // 
    // Note: new entities used in subsequent connectivity lists, sets, or tags, 
    //   are referred to as (MBMAXTYPE + index), where index is into vector 
    //   of new entities, 0-based
    //
    // DATA LAYOUT IN BUFFER:
    // . w/ handles ? (0-no, 1-yes, range, 2-yes, vector)
    // OWNING PROC INFO:
    // . # owning proc tuples
    // . for each tuple:
    //   - index in entities in this message
    //   - proc owning this entity
    //   - handle on owning proc
    // ENTITIES:
    // . for all types:
    //   - ent_type -- end of data indicated by ent_type == MBMAXTYPE
    //   - range(ent type) -- #ents = range.size()
    //   - if (ent_type == MBVERTEX) xxx[#ents], yyy[#ents], zzz[#ents]
    //   - else {
    //     . nodes_per_ent
    //     . connect[#ents * nodes_per_ent]
    //   - }
    //   - if (handles) range/vector of remote handles
    // SETS:
    // . range of set handles -- (#sets = range.size())
    // . options[#sets] (unsigned int)
    // . if (unordered) set range 
    // . else if ordered
    //   - #ents in set
    //   - handles[#ents]
    // . #parents
    // . if (#parents) handles[#parents]
    // . #children
    // . if (#children) handles[#children]
    // ADJACENCIES:
    // (no adjs for now)
    // TAGS:
    // . #tags
    // . for each tag
    //   - tag size
    //   - tag type
    //   - data type
    //   - if (default value)
    //     . 1
    //     . default value
    //   - else
    //     . 0
    //   - name (null-terminated string, in char[64])
    //   - range (size = #ents)
    //   - tag_vals[#ents]

  MBErrorCode result;

  MBRange set_range;
  std::vector<MBRange> set_ranges;
  std::vector<MBTag> all_tags;
  std::vector<MBRange> tag_ranges;
  std::vector<int> set_sizes;
  std::vector<unsigned int> options_vec;

  MBRange::const_iterator rit;

    // get an estimate of the buffer size
  buff_size = estimate_ents_buffer_size(orig_ents, store_remote_handles);
  buff.clear();
  buff.resize(buff_size);

  unsigned char *buff_ptr = &buff[0];
  
    // entities
  result = pack_entities(orig_ents, buff, buff_ptr,
                         store_remote_handles, to_proc, true, owned_shared); 
  RRA("Packing entities failed.");
  
    // sets
  result = pack_sets(orig_ents, buff, buff_ptr,
                     store_remote_handles, to_proc); 
  RRA("Packing sets (count) failed.");

    // tags
  MBRange final_ents;
  if (tags) {
    result = get_tag_send_list(orig_ents, all_tags, tag_ranges );
    RRA("Failed to get tagged entities.");
    result = pack_tags(orig_ents, all_tags, tag_ranges, 
                       buff, buff_ptr, store_remote_handles, to_proc);
    RRA("Packing tags (count) failed.");
  }

  return result;
}
 
MBErrorCode MBParallelComm::unpack_buffer(unsigned char *buff_ptr,
                                          const bool store_remote_handles,
                                          const int from_proc,
                                          const int ind,
                                          MBRange &final_ents,
                                          MBRange *nonowned_ghosts,
                                          std::set<unsigned int> &my_addl_procs)
{
  if (myBuffer.capacity() == 0) return MB_FAILURE;
  
#ifdef DEBUG_PACKING
    unsigned char *tmp_buff = buff_ptr;
#endif  
    MBErrorCode result;
    result = unpack_entities(buff_ptr, store_remote_handles,
                             from_proc, ind, final_ents, nonowned_ghosts, my_addl_procs);
  RRA("Unpacking entities failed.");
#ifdef DEBUG_PACKING
    std::cerr << "unpack_entities buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  
  result = unpack_sets(buff_ptr, final_ents, store_remote_handles, 
                       from_proc);
  RRA("Unpacking sets failed.");
#ifdef DEBUG_PACKING
    std::cerr << "unpack_sets buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  
  result = unpack_tags(buff_ptr, final_ents, store_remote_handles,
                       from_proc);
  RRA("Unpacking tags failed.");
#ifdef DEBUG_PACKING
    std::cerr << "unpack_tags buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  

#ifdef DEBUG_PACKING
  std::cerr << std::endl;
#endif
  
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

int MBParallelComm::estimate_ents_buffer_size(MBRange &entities,
                                              const bool store_remote_handles) 
{
  int buff_size = 0;
  std::vector<MBEntityHandle> dum_connect_vec;
  const MBEntityHandle *connect;
  int num_connect;

  int num_verts = entities.num_of_type(MBVERTEX);
    // # verts + coords + handles
  buff_size += 2*sizeof(int) + 3*sizeof(double)*num_verts;
  if (store_remote_handles) buff_size += sizeof(MBEntityHandle)*num_verts;

    // do a rough count by looking at first entity of each type
  for (MBEntityType t = MBEDGE; t < MBENTITYSET; t++) {
    const MBRange::iterator rit = entities.lower_bound(t);
    if (TYPE_FROM_HANDLE(*rit) != t) continue;
    
    MBErrorCode result = mbImpl->get_connectivity(*rit, connect, num_connect, 
                                                  true, &dum_connect_vec);
    RRA("Failed to get connectivity to estimate buffer size.");

      // number, type, nodes per entity
    buff_size += 3*sizeof(int);
    int num_ents = entities.num_of_type(t);
      // connectivity, handle for each ent
    buff_size += (num_connect+1)*sizeof(MBEntityHandle)*num_ents;
  }

      // extra entity type at end, passed as int
  buff_size += sizeof(int);

  return buff_size;
}

int MBParallelComm::estimate_sets_buffer_size(MBRange &entities,
                                              const bool store_remote_handles) 
{
    // number of sets
  int buff_size = sizeof(int);
  
    // do a rough count by looking at first entity of each type
  MBRange::iterator rit = entities.lower_bound(MBENTITYSET);
  MBErrorCode result;
  
  for (; rit != entities.end(); rit++) {
    unsigned int options;
    result = mbImpl->get_meshset_options(*rit, options);
    RRA("Failed to get meshset options.");

    buff_size += sizeof(int);
    
    MBRange set_range;
    if (options & MESHSET_SET) {
        // range-based set; count the subranges
      result = mbImpl->get_entities_by_handle(*rit, set_range);
      RRA("Failed to get set entities.");

        // set range
      buff_size += RANGE_SIZE(set_range);
    }
    else if (options & MESHSET_ORDERED) {
        // just get the number of entities in the set
      int num_ents;
      result = mbImpl->get_number_entities_by_handle(*rit, num_ents);
      RRA("Failed to get number entities in ordered set.");

        // set vec
      buff_size += sizeof(MBEntityHandle) * num_ents + sizeof(int);
    }

      // get numbers of parents/children
    int num_par, num_ch;
    result = mbImpl->num_child_meshsets(*rit, &num_ch);
    RRA("Failed to get num children.");

    result = mbImpl->num_parent_meshsets(*rit, &num_par);
    RRA("Failed to get num parents.");

    buff_size += (num_ch + num_par) * sizeof(MBEntityHandle) + 2*sizeof(int);
  }

  return buff_size;
}

MBErrorCode MBParallelComm::pack_entities(MBRange &entities,
                                          std::vector<unsigned char> &buff,
                                          unsigned char *&buff_ptr,
                                          const bool store_remote_handles,
                                          const int to_proc,
                                          const bool is_ghost,
                                          MBRange &owned_shared) 
{
    // pack a range of entities into a message, appending to buffer;
    // buffer contains:
    // . # non-owned (by sender) entities
    // . for each NO entity: index in new ents, owning proc
    // . for each NO entity: handle on owner
    // . repeat:
    //   - ent type
    //   - # ents
    //   - if type == vertex, xxx[#ents], yyy[#ents], zzz[#ents]
    //   - else
    //     . # nodes per entity
    //     . connect[#ents * nodes_per_ent]
    // . until type == MBMAXTYPE
    // . if (store_remote_handles) range of local handles on sending proc

  MBWriteUtilIface *wu;
  MBErrorCode result = mbImpl->query_interface(std::string("MBWriteUtilIface"), 
                                               reinterpret_cast<void**>(&wu));
  RRA("Couldn't get MBWriteUtilIface.");

    // non-owned entities don't matter if remote handles aren't stored
    // or if passing interface layer
  if (store_remote_handles && is_ghost) {
      // get owned, non-owned ents
    MBRange owned_ents, nonowned_ents;
    result = filter_owned_shared(entities, 
                                 true, true, false, false,
                                 -1, &owned_ents);
    RRA("Failed to get owned entities.");
    if (owned_ents.size() != entities.size())
      nonowned_ents = entities.subtract(owned_ents);
  
      // pack index, proc, handle for nonowned ents
    CHECK_BUFF_SPACE(buff, buff_ptr,
                     sizeof(int) + (2*sizeof(int) + 
                                    sizeof(MBEntityHandle))*nonowned_ents.size());

    PACK_INT(buff_ptr, nonowned_ents.size());
    if (!nonowned_ents.empty()) {
      int proc;
      std::vector<MBEntityHandle> handles;
      handles.reserve(nonowned_ents.size());
      MBEntityHandle handle;
      for (MBRange::const_iterator rit = nonowned_ents.begin();
           rit != nonowned_ents.end(); rit++) {
        int idx = entities.index(*rit);
        assert(-1 != idx);
        result = get_owner_handle(*rit, proc, handle);
        RRA("Failed to get owner and handle for entity.");
        assert(proc != (int) procConfig.proc_rank() && 
               to_proc != proc && 0 != handle);
        PACK_INT(buff_ptr, idx);
        PACK_INT(buff_ptr, proc);
        handles.push_back(handle);
      }
      PACK_EH(buff_ptr, &handles[0], nonowned_ents.size());
    }

      // put owned shared in sharedEnts list
    owned_shared.merge(owned_ents);
  }

  if (store_remote_handles) {
      // pack the local handles
    int tmp_space = RANGE_SIZE(entities);
    CHECK_BUFF_SPACE(buff, buff_ptr, tmp_space);
    PACK_RANGE(buff_ptr, entities);
  }

    // pack vertices
  MBRange these_ents = entities.subset_by_type(MBVERTEX);
  int num_ents = these_ents.size();
  int buff_size;

    // if iface, shouldn't be any vertices here
  assert(is_ghost || num_ents == 0);
  
  if (num_ents) {
    buff_size = 2*sizeof(int) + 3*num_ents*sizeof(double);
    CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);

    // type, # ents
    PACK_INT(buff_ptr, ((int) MBVERTEX));
    PACK_INT(buff_ptr, ((int) num_ents));

    std::vector<double*> coords(3);
    for (int i = 0; i < 3; i++)
      coords[i] = reinterpret_cast<double*>(buff_ptr + i * num_ents * sizeof(double));
    assert(NULL != wu);
    result = wu->get_node_arrays(3, num_ents, these_ents, 0, 0, coords);
    RRA("Couldn't allocate node space.");
    PC(3*num_ents, " doubles");

    buff_ptr += 3 * num_ents * sizeof(double);

#ifdef DEBUG_PACKING
  std::cerr << "Packed " << these_ents.size() << " ents of type " 
            << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*these_ents.begin())) << std::endl;
#endif      
  }

    // now entities; go through range, packing by type and equal # verts per element
  MBRange::iterator start_rit = entities.find(*these_ents.rbegin());
  start_rit++;
  int last_nodes = -1;
  MBEntityType last_type = MBMAXTYPE;
  these_ents.clear();
  MBRange::iterator end_rit = start_rit;
  EntitySequence *seq;
  ElementSequence *eseq;
  
  while (start_rit != entities.end() || !these_ents.empty()) {
      // cases:
      // A: !end, last_type == MAXTYPE, seq: save contig sequence in these_ents
      // B: !end, last type & nodes same, seq: save contig sequence in these_ents
      // C: !end, last type & nodes different: pack these_ents, then save contig sequence in these_ents
      // D: end: pack these_ents

      // find the sequence holding current start entity, if we're not at end
    eseq = NULL;
    if (start_rit != entities.end()) {
      result = sequenceManager->find(*start_rit, seq);
      RRA("Couldn't find entity sequence.");
      if (NULL == seq) return MB_FAILURE;
      eseq = dynamic_cast<ElementSequence*>(seq);
    }

      // pack the last batch if at end or next one is different
    if (!these_ents.empty() &&
        (!eseq || eseq->type() != last_type ||
         last_nodes != (int) eseq->nodes_per_element())) {
      result = pack_entity_seq(last_nodes, store_remote_handles,
                               to_proc, these_ents, entities, buff, buff_ptr);
      RRA("Failed to pack entities from a sequence.");
      these_ents.clear();
    }

    if (eseq) {
        // continuation of current range, just save these entities
        // get position in entities list one past end of this sequence
      end_rit = entities.lower_bound(start_rit, entities.end(), eseq->end_handle()+1);

        // put these entities in the range
      std::copy(start_rit, end_rit, mb_range_inserter(these_ents));

      last_type = eseq->type();
      last_nodes = eseq->nodes_per_element();
    }
    else if (start_rit != entities.end() &&
             TYPE_FROM_HANDLE(*start_rit) == MBENTITYSET)
      break;

    start_rit = end_rit;
  }

    // pack MBMAXTYPE to indicate end of ranges
  CHECK_BUFF_SPACE(buff, buff_ptr, sizeof(int));
  PACK_INT(buff_ptr, ((int)MBMAXTYPE));

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_addl_procs(std::set<unsigned int> &addl_procs,
                                            std::vector<unsigned char> &buff,
                                            unsigned char *&buff_ptr) 
{
  int buff_size = (1 + addl_procs.size())*sizeof(int);
  CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);

    // pack extra procs to_proc must communicate with
  PACK_INT(buff_ptr, addl_procs.size());
  
  for (std::set<unsigned int>::iterator sit = addl_procs.begin(); 
       sit != addl_procs.end(); sit++)
    PACK_INT(buff_ptr, *sit);

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_addl_procs(unsigned char *&buff_ptr,
                                              std::set<unsigned int> &addl_procs) 
{
  int num_addl, tmp_addl;
  UNPACK_INT(buff_ptr, num_addl);
  for (int i = 0; i < num_addl; i++) {
    UNPACK_INT(buff_ptr, tmp_addl);
    addl_procs.insert(tmp_addl);
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_entity_seq(const int nodes_per_entity,
                                            const bool store_remote_handles,
                                            const int to_proc,
                                            MBRange &these_ents,
                                            MBRange &entities,
                                            std::vector<unsigned char> &buff,
                                            unsigned char *&buff_ptr) 
{
  int tmp_space = 3*sizeof(int) + nodes_per_entity*these_ents.size()*sizeof(MBEntityHandle);
  CHECK_BUFF_SPACE(buff, buff_ptr, tmp_space);
  
    // pack the entity type
  PACK_INT(buff_ptr, ((int)TYPE_FROM_HANDLE(*these_ents.begin())));

    // pack # ents
  PACK_INT(buff_ptr, these_ents.size());
      
    // pack the nodes per entity
  PACK_INT(buff_ptr, nodes_per_entity);
      
    // pack the connectivity
  const MBEntityHandle *connect;
  int num_connect;
  std::vector<MBEntityHandle> dum_connect;
  MBEntityHandle *start_vec = (MBEntityHandle*)buff_ptr;
  MBErrorCode result = MB_SUCCESS;
  for (MBRange::const_iterator rit = these_ents.begin(); rit != these_ents.end(); rit++) {
    result = mbImpl->get_connectivity(*rit, connect, num_connect, false,
                                      &dum_connect);
    RRA("Failed to get connectivity.");
    assert(num_connect == nodes_per_entity);
    PACK_EH(buff_ptr, connect, num_connect);
  }

    // substitute destination handles
  result = get_remote_handles(store_remote_handles, start_vec, start_vec,
                              nodes_per_entity*these_ents.size(), to_proc,
                              entities);
  RRA("Trouble getting remote handles when packing entities.");

#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Packed " << these_ents.size() << " ents of type " 
            << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*these_ents.begin())) << std::endl;
#endif      

  return result;
}


MBErrorCode MBParallelComm::get_remote_handles(const bool store_remote_handles,
                                               MBEntityHandle *from_vec, 
                                               MBEntityHandle *to_vec_tmp,
                                               int num_ents, int to_proc,
                                               const MBRange &new_ents) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE RANGE-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  if (0 == num_ents) return MB_SUCCESS;
  
    // use a local destination ptr in case we're doing an in-place copy
  std::vector<MBEntityHandle> tmp_vector;
  MBEntityHandle *to_vec = to_vec_tmp;
  if (to_vec == from_vec) {
    tmp_vector.resize(num_ents);
    to_vec = &tmp_vector[0];
  }

  if (!store_remote_handles) {
    int err;
      // in this case, substitute position in new_ents list
    for (int i = 0; i < num_ents; i++) {
      int ind = new_ents.index(from_vec[i]);
      to_vec[i] = CREATE_HANDLE(MBMAXTYPE, ind, err);
      assert(to_vec[i] != 0 && !err && -1 != ind);
    }
  }
  else {
    MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
    MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                              sharedh_tag, sharedhs_tag, pstatus_tag);
  
      // get single-proc destination handles and shared procs
    std::vector<int> sharing_procs(num_ents);
    result = mbImpl->tag_get_data(sharedh_tag, from_vec, num_ents,
                                  to_vec);
    RRA("Failed to get shared handle tag for remote_handles.");
    result = mbImpl->tag_get_data(sharedp_tag, from_vec, num_ents, &sharing_procs[0]);
    RRA("Failed to get sharing proc tag in remote_handles.");
    for (int j = 0; j < num_ents; j++) {
      if (to_vec[j] && sharing_procs[j] != to_proc)
        to_vec[j] = 0;
    }
    
    MBEntityHandle tmp_handles[MAX_SHARING_PROCS];
    int tmp_procs[MAX_SHARING_PROCS];
    int i;
      // go through results, and for 0-valued ones, look for multiple shared proc
    MBEntityHandle *tmp_eh;
    for (tmp_eh = to_vec, i = 0; i < num_ents; i++) {
      if (!to_vec[i]) {
        result = mbImpl->tag_get_data(sharedps_tag, from_vec+i, 1, tmp_procs);
        if (MB_SUCCESS == result) {
          for (int j = 0; j < MAX_SHARING_PROCS; j++) {
            if (-1 == tmp_procs[j]) break;
            else if (tmp_procs[j] == to_proc) {
              result = mbImpl->tag_get_data(sharedhs_tag, from_vec+i, 1, tmp_handles);
              RRA("Trouble getting sharedhs tag.");
              to_vec[i] = tmp_handles[j];
              assert(to_vec[i]);
              break;
            }
          }
        }
        if (!to_vec[i]) {
          int j = new_ents.index(from_vec[i]);
          if (-1 == j) {
            result = MB_FAILURE;
            RRA("Failed to find new entity in send list.");
          }
          int err;
          to_vec[i] = CREATE_HANDLE(MBMAXTYPE, j, err);
          if (err) {
            result = MB_FAILURE;
            RRA("Failed to create handle in remote_handles.");
          }
        }
      }
    }
  }
  
    // memcpy over results if from_vec and to_vec are the same
  if (to_vec_tmp == from_vec) 
    memcpy(from_vec, to_vec, num_ents * sizeof(MBEntityHandle));
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_remote_handles(const bool store_remote_handles,
                                               const MBRange &from_range, 
                                               MBEntityHandle *to_vec,
                                               int to_proc,
                                               const MBRange &new_ents) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE VECTOR-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  if (from_range.empty()) return MB_SUCCESS;
  
  if (!store_remote_handles) {
    int err;
      // in this case, substitute position in new_ents list
    MBRange::iterator rit;
    unsigned int i;
    for (rit = from_range.begin(), i = 0; rit != from_range.end(); rit++, i++) {
      int ind = new_ents.index(*rit);
      to_vec[i] = CREATE_HANDLE(MBMAXTYPE, ind, err);
      assert(to_vec[i] != 0 && !err && -1 != ind);
    }
  }
  else {
    MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
    MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                              sharedh_tag, sharedhs_tag, pstatus_tag);
  
      // get single-proc destination handles and shared procs
    std::vector<int> sharing_procs(from_range.size());
    result = mbImpl->tag_get_data(sharedh_tag, from_range, to_vec);
    RRA("Failed to get shared handle tag for remote_handles.");
    result = mbImpl->tag_get_data(sharedp_tag, from_range, &sharing_procs[0]);
    RRA("Failed to get sharing proc tag in remote_handles.");
    for (unsigned int j = 0; j < from_range.size(); j++) {
      if (to_vec[j] && sharing_procs[j] != to_proc)
        to_vec[j] = 0;
    }
    
    MBEntityHandle tmp_handles[MAX_SHARING_PROCS];
    int tmp_procs[MAX_SHARING_PROCS];
      // go through results, and for 0-valued ones, look for multiple shared proc
    MBRange::iterator rit;
    unsigned int i;
    for (rit = from_range.begin(), i = 0; rit != from_range.end(); rit++, i++) {
      if (!to_vec[i]) {
        result = mbImpl->tag_get_data(sharedhs_tag, &(*rit), 1, tmp_handles);
        if (MB_SUCCESS == result) {
          result = mbImpl->tag_get_data(sharedps_tag, &(*rit), 1, tmp_procs);
          RRA("Trouble getting sharedps tag.");
          for (int j = 0; j < MAX_SHARING_PROCS; j++)
            if (tmp_procs[j] == to_proc) {
              to_vec[i] = tmp_handles[j];
              break;
            }
        }
      
        if (!to_vec[i]) {
          int j = new_ents.index(*rit);
          if (-1 == j) {
            result = MB_FAILURE;
            RRA("Failed to find new entity in send list.");
          }
          int err;
          to_vec[i] = CREATE_HANDLE(MBMAXTYPE, j, err);
          if (err) {
            result = MB_FAILURE;
            RRA("Failed to create handle in remote_handles.");
          }
        }
      }
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_remote_handles(const bool store_remote_handles,
                                               const MBRange &from_range, 
                                               MBRange &to_range,
                                               int to_proc,
                                               const MBRange &new_ents) 
{
  std::vector<MBEntityHandle> to_vector(from_range.size());

  MBErrorCode result =
    get_remote_handles(store_remote_handles, from_range, &to_vector[0],
                       to_proc, new_ents);
  RRA("Trouble getting remote handles.");
  std::copy(to_vector.begin(), to_vector.end(), mb_range_inserter(to_range));
  return result;
}

MBErrorCode MBParallelComm::unpack_entities(unsigned char *&buff_ptr,
                                            const bool store_remote_handles,
                                            const unsigned int from_proc,
                                            const int from_ind,
                                            MBRange &entities,
                                            MBRange *nonowned_ghosts,
                                            std::set<unsigned int> &my_addl_procs) 
{
  MBErrorCode result;
  bool done = false;
  MBReadUtilIface *ru = NULL;
  result = mbImpl->query_interface(std::string("MBReadUtilIface"), 
                                   reinterpret_cast<void**>(&ru));
  RRA("Failed to get MBReadUtilIface.");

  MBRange remote_range;

    // unpack any nonowned entities
  int no_size = 0;
  std::vector<unsigned int> no_indices, no_procs, o_indices;
  std::vector<MBEntityHandle> no_handles, no_local_handles,
      o_local_handles, o_handles;
  if (store_remote_handles) {
    UNPACK_INT(buff_ptr, no_size);

    no_indices.resize(no_size);
    no_procs.resize(no_size);
    no_handles.resize(no_size);

      // no_local_handles[i] = local entity equivalent to new entity i in
      //     this message
    for (int i = 0; i < no_size; i++) {
      UNPACK_INT(buff_ptr, no_indices[i]);
      UNPACK_INT(buff_ptr, no_procs[i]);
      UNPACK_EH(buff_ptr, &no_handles[i], 1);
    }

      // unpack source handles; will be processed later
    UNPACK_RANGE(buff_ptr, remote_range);

    no_local_handles.resize(remote_range.size(), 0);

    for (int i = 0; i < no_size; i++) {
      no_local_handles[no_indices[i]] = 
          find_ghost_entity(no_handles[i], no_procs[i]);
    }

      // look for existing ghost entities from nonowned procs
    if (!sharedEnts[from_ind].noHandles.empty()) {
      MBRange tmp_range = sharedEnts[from_ind].noHandles.intersect(remote_range);
      for (MBRange::iterator rit = tmp_range.begin(); 
           rit != tmp_range.end(); rit++) {
        int idx = remote_range.index(*rit);
        no_local_handles[idx] = find_ghost_entity(*rit, from_proc);
      }
    }
  }

  MBRange new_ents;
  
  while (!done) {
    MBEntityType this_type = MBMAXTYPE;
    MBEntityHandle actual_start;
    UNPACK_INT(buff_ptr, this_type);
    assert(this_type >= MBVERTEX && 
           (this_type == MBMAXTYPE || this_type < MBENTITYSET));

      // MBMAXTYPE signifies end of entities data
    if (MBMAXTYPE == this_type) break;
    
      // get the number of ents
    int num_ents;
    UNPACK_INT(buff_ptr, num_ents);

    if (MBVERTEX == this_type) {
        // unpack coords
      if (num_ents) {
        std::vector<double*> coords(3);
        result = ru->get_node_arrays(3, num_ents, 0, 
                                     actual_start, coords);
        RRA("Failed to allocate node arrays.");

        new_ents.insert(actual_start, actual_start+num_ents-1);
      
          // unpack the buffer data directly into coords
        for (int i = 0; i < 3; i++) 
          memcpy(coords[i], buff_ptr+i*num_ents*sizeof(double), 
                 num_ents*sizeof(double));
        buff_ptr += 3*num_ents * sizeof(double);
        UPC(3*num_ents, " doubles");
      }
    }

    else {
      int verts_per_entity;
      
        // unpack the nodes per entity
      UNPACK_INT(buff_ptr, verts_per_entity);
      
      MBEntityHandle *connect;
      result = ru->get_element_array(num_ents, verts_per_entity, this_type,
                                     0, actual_start,
                                     connect);
      RRA("Failed to allocate element arrays.");

        // unpack the connectivity
      UNPACK_EH(buff_ptr, connect, (num_ents*verts_per_entity));
      new_ents.insert(actual_start, actual_start+num_ents-1);

        // convert to local handles
      result = get_local_handles(connect, num_ents*verts_per_entity,
                                 new_ents, &no_local_handles);
      RRA("Couldn't get local handles.");

        // update adjacencies
      result = ru->update_adjacencies(actual_start, num_ents, 
                                      verts_per_entity, connect);
      RRA("Failed to update adjacencies.");
    }

#ifdef DEBUG_PACKING
      std::cerr << "Unpacked " << num_ents << " ents of type " 
                << MBCN::EntityTypeName(TYPE_FROM_HANDLE(actual_start)) << std::endl;
#endif      

  }

  if (store_remote_handles) {
      // unpack source handles
    result = set_remote_data(new_ents, remote_range, from_proc);
    RRA("Couldn't set sharing data");

      // account for entities which existed before (i.e. delete duplicate 
      // new ones) then save remaining new ones in appropriate list

      // find ones to delete (existing nonowned and owned entities)
    MBRange to_delete;
    for (unsigned int i = 0; i < no_local_handles.size(); i++)
      if (no_local_handles[i]) to_delete.insert(new_ents[i]);

      // put non-existing non-owned entities on noHandles and localHandles list
    int ind = -1;
    for (int i = 0; i < no_size; i++) {
      if (!no_local_handles[no_indices[i]]) {
        if (-1 == ind || no_procs[i] != buffProcs[ind])
          ind = get_buffers(no_procs[i]);
        sharedEnts[ind].noHandles.insert(no_handles[i]);
        sharedEnts[ind].localHandles.insert(new_ents[no_indices[i]]);
          // remove remote handles for these entities; leaves 
          // remote handles for only new owned entities in remote range
        remote_range.erase(new_ents[no_indices[i]]);
      }
    }

    result = mbImpl->delete_entities(to_delete);
    RRA("Failed to delete non-owned duplicate entities.");

    entities = new_ents.subtract(to_delete);
  
      // put new ents on local and remote handles list for from_proc
    sharedEnts[from_ind].localHandles.merge(entities);
    sharedEnts[from_ind].remoteHandles.merge(remote_range);

    result = set_pstatus_entities(entities, PSTATUS_NOT_OWNED, false, false);
    RRA("Couldn't set pstatus values for received entities.");
  }
  else 
    entities.swap(new_ents);
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking entities." << std::endl;
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_iface_entities(unsigned char *&buff_ptr,
                                                  const int from_proc, 
                                                  const int ind,
                                                  std::vector<MBEntityHandle> &recd_ents) 
{
    // for all the entities in the received buffer, save
    // entities in this instance which match connectivity, or zero if none found
  MBErrorCode result;
  bool done = false;
  std::vector<MBEntityHandle> connect;
  unsigned int num_zero = 0;

  MBRange remote_range;
  UNPACK_RANGE(buff_ptr, remote_range);

  while (!done) {
    MBEntityType this_type;
    UNPACK_INT(buff_ptr, this_type);
    assert(this_type > MBVERTEX && 
           (this_type == MBMAXTYPE || this_type < MBENTITYSET));

      // MBMAXTYPE signifies end of entities data
    if (MBMAXTYPE == this_type) break;
    
      // get the number of ents
    int num_ents;
    UNPACK_INT(buff_ptr, num_ents);
    
    int verts_per_entity;
      
      // unpack the verts per entity
    UNPACK_INT(buff_ptr, verts_per_entity);

    connect.resize(verts_per_entity * num_ents);
    
      // unpack the connectivity
    UNPACK_EH(buff_ptr, &connect[0], (num_ents*verts_per_entity));

      // unpack source handles
    MBRange tmp_range;

      // save matching local entities, or 0 if none found
    for (int i = 0; i < num_ents; i++) {
        // test for existing entity
      result = mbImpl->get_adjacencies(&connect[0]+i*verts_per_entity, verts_per_entity, 
                                       MBCN::Dimension(this_type), false, tmp_range);
      if (MB_MULTIPLE_ENTITIES_FOUND == result) {
        RRA("Multiple entities found.");
      }
      if (!tmp_range.empty()) {
          // found a corresponding entity - save target
        assert(1 == tmp_range.size());
        recd_ents.push_back(*tmp_range.begin());
        tmp_range.clear();
      }
      else {
        recd_ents.push_back(0);
        num_zero++;
      }
    }

#ifdef DEBUG_PACKING
    std::cerr << "Unpacked " << num_ents << " ents of type " 
              << MBCN::EntityTypeName(this_type) << std::endl;
#endif      

  }
  
  if (num_zero) {
    result = MB_FAILURE;
    RRA("Didn't find an iface entity.");
  }
  std::vector<MBEntityHandle> tmp_remote;
  tmp_remote.reserve(recd_ents.size());
  std::copy(remote_range.begin(), remote_range.end(), 
            std::back_inserter(tmp_remote));
  result = set_remote_data(&recd_ents[0], &tmp_remote[0], recd_ents.size(),
                           from_proc);
  RRA("Failed to set remote handles for iface entities.");

    // add iface entities to sharedEnts list; since they were sent to us, assume we
    // don't own them
  std::copy(recd_ents.begin(), recd_ents.end(), 
            mb_range_inserter(sharedEnts[ind].localHandles));
  sharedEnts[ind].remoteHandles.merge(remote_range);
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking iface entities." << std::endl;
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_local_handles(const MBRange &remote_handles,
                                              MBRange &local_handles,
                                              const MBRange &new_ents) 
{
  std::vector<MBEntityHandle> rh_vec;
  rh_vec.reserve(remote_handles.size());
  std::copy(remote_handles.begin(), remote_handles.end(), std::back_inserter(rh_vec));
  MBErrorCode result = get_local_handles(&rh_vec[0], remote_handles.size(), new_ents);
  std::copy(rh_vec.begin(), rh_vec.end(), mb_range_inserter(local_handles));
  return result;
}
  
MBErrorCode MBParallelComm::get_local_handles(MBEntityHandle *from_vec, 
                                              int num_ents,
                                              const MBRange &new_ents,
                                              std::vector<MBEntityHandle> *no_ents) 
{
  for (int i = 0; i < num_ents; i++) {
    if (TYPE_FROM_HANDLE(from_vec[i]) == MBMAXTYPE) {
      assert(ID_FROM_HANDLE(from_vec[i]) < (int) new_ents.size());
      if (no_ents && i < (int)no_ents->size() && (*no_ents)[i])
        from_vec[i] = (*no_ents)[i];
      else
        from_vec[i] = new_ents[ID_FROM_HANDLE(from_vec[i])];
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::set_remote_data(MBRange &local_range,
                                            MBRange &remote_range,
                                            int other_proc) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE VECTOR-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!

  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);

    // get remote procs tag, if any
  MBRange tmp_range, tmp_local_range = local_range;
  std::vector<int> remote_proc(local_range.size());
  int remote_procs[MAX_SHARING_PROCS];
  std::vector<MBEntityHandle> remote_handle(local_range.size());
  MBEntityHandle remote_handles[MAX_SHARING_PROCS];
  std::fill(remote_procs, remote_procs+MAX_SHARING_PROCS, -1);
  std::fill(remote_handles, remote_handles+MAX_SHARING_PROCS, 0);
  result = mbImpl->tag_get_data(sharedp_tag, local_range,
                                &remote_proc[0]);
  RRA("Couldn't get sharedp tag (range).");
  result = mbImpl->tag_get_data(sharedh_tag, local_range,
                                &remote_handle[0]);
  RRA("Couldn't get sharedh tag (range).");
  MBRange::iterator rit, rit2;
  int i = 0;

    // for each pair of local/remote handles:
  for (rit = tmp_local_range.begin(), rit2 = remote_range.begin(); 
       rit != tmp_local_range.end(); rit++, rit2++, i++) {

      // get existing remote proc(s), handle(s) for this local handle
    remote_procs[0] = remote_proc[i];
    if (-1 != remote_procs[0]) remote_handles[0] = remote_handle[i];
    else {
      result = mbImpl->tag_get_data(sharedps_tag, &(*rit), 1,
                                    remote_procs);
      if (MB_SUCCESS == result) {
        result = mbImpl->tag_get_data(sharedhs_tag, &(*rit), 1,
                                      remote_handles);
        RRA("Couldn't get sharedhs tag (range).");
      }
    }

    result = add_remote_proc(*rit, remote_procs, remote_handles,
                             other_proc, *rit2);
    RRA(" ");
  }

    // also update shared flag for these ents
  result = set_pstatus_entities(local_range, PSTATUS_SHARED, false, false);
  RRA("Couldn't set pstatus tag (range)");  
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::set_remote_data(MBEntityHandle *local_ents,
                                            MBEntityHandle *remote_ents,
                                            int num_ents,
                                            int other_proc) 
{
  if (0 == num_ents)
    return MB_SUCCESS;

    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE RANGE-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);

    // get remote procs tag, if any
  MBRange tmp_range;
  std::vector<int> remote_proc(num_ents);
  int remote_procs[MAX_SHARING_PROCS];
  std::vector<MBEntityHandle> remote_handle(num_ents);
  MBEntityHandle remote_handles[MAX_SHARING_PROCS];
  std::vector<unsigned char> pstatus_vals(num_ents);
  result = mbImpl->tag_get_data(sharedp_tag, local_ents, num_ents,
                                &remote_proc[0]);
  RRA("Couldn't get sharedp tag (vector)");
  result = mbImpl->tag_get_data(sharedh_tag, local_ents, num_ents,
                                &remote_handle[0]);
  RRA("Couldn't get sharedh tag (vector)");
  result = mbImpl->tag_get_data(pstatus_tag, local_ents, num_ents,
                                &pstatus_vals[0]);
  RRA("Couldn't get sharedh tag (vector)");

    // for each local/remote handle pair
  for (int i = 0; i != num_ents; i++) {

      // get existing remote proc(s), handle(s) for this local handle
    if (!(pstatus_vals[i] & PSTATUS_SHARED)) {
      std::fill(remote_procs, remote_procs+MAX_SHARING_PROCS, -1);
      std::fill(remote_handles, remote_handles+MAX_SHARING_PROCS, 0);
    }
    else if (-1 != remote_proc[i]) {
      remote_procs[0] = remote_proc[i];
      remote_handles[0] = remote_handle[i];
      std::fill(remote_procs+1, remote_procs+MAX_SHARING_PROCS-1, -1);
      std::fill(remote_handles+1, remote_handles+MAX_SHARING_PROCS-1, 0);
    }
    else {
      result = mbImpl->tag_get_data(sharedps_tag, local_ents+i, 1,
                                    remote_procs);
      RRA("Couldn't get sharedps tag (vector)");
      
      result = mbImpl->tag_get_data(sharedhs_tag, local_ents+i, 1,
                                    remote_handles);
      RRA("Couldn't get sharedhs tag (vector)");
    }

      // now either insert other_proc, handle into these, or remove if
      // remote handle is 0
    if (0 == remote_ents[i]) {
      result = rmv_remote_proc(local_ents[i], remote_procs, remote_handles,
                               other_proc);
      RRA(" ");
    }
    else {
      result = add_remote_proc(local_ents[i], remote_procs, remote_handles,
                               other_proc, remote_ents[i]);
      RRA(" ");
    }
  }

    // also update shared flag for these ents
  result = set_pstatus_entities(local_ents, num_ents, PSTATUS_SHARED,
                                false, false);
  RRA("Couldn't set pstatus tag (range)");  

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::rmv_remote_proc(MBEntityHandle ent,
                                            int *remote_procs,
                                            MBEntityHandle *remote_hs,
                                            int remote_proc) 
{
  int i = std::find(remote_procs, remote_procs+MAX_SHARING_PROCS, -1) - remote_procs;
  int j = std::find(remote_procs, remote_procs+i, remote_proc) - remote_procs;
  MBErrorCode result;

    // remote_proc must be in this list
  assert(j != i);
  
  for (; j != i; j++) {
    remote_procs[j] = remote_procs[j+1];
    remote_hs[j] = remote_hs[j+1];
  }

  if (i == 1) {
      // only had one proc, and now have none; have to set the values,
      // since they're sparse tags
    result = mbImpl->tag_set_data(sharedp_tag(), &ent, 1, remote_procs);
    result = mbImpl->tag_set_data(sharedh_tag(), &ent, 1, remote_hs);
      // if it's not shared, it's also not a ghost, interface, or !owned entity
    unsigned char not_shared = 0;
    result = mbImpl->tag_set_data(pstatus_tag(), &ent, 1, &not_shared);
  }
  else if (i == 2) {
      // went from 2 to 1, need to unset 1 and set the other
    result = mbImpl->tag_set_data(sharedp_tag(), &ent, 1, remote_procs);
    RRA("Couldn't set sharedp tag");
    result = mbImpl->tag_set_data(sharedh_tag(), &ent, 1, remote_hs);
    RRA("Couldn't set sharedh tag");
    result = mbImpl->tag_delete_data(sharedps_tag(), &ent, 1);
    RRA("Couldn't remove sharedps tag");
    result = mbImpl->tag_delete_data(sharedhs_tag(), &ent, 1);
    RRA("Couldn't remove sharedhs tag");
  }
  else {
    result = mbImpl->tag_set_data(sharedps_tag(), &ent, 1, remote_procs);
    RRA("Couldn't set sharedps tag");
    result = mbImpl->tag_set_data(sharedhs_tag(), &ent, 1, remote_hs);
    RRA("Couldn't set sharedhs tag");
  }
  
  return MB_SUCCESS;
}

template <typename T> void
insert_in_array( T* array, size_t array_size, size_t location, T value )
{
  assert( location+1 < array_size );
  for (size_t i = array_size-1; i > location; --i)
    array[i] = array[i-1];
  array[location] = value;
}

MBErrorCode MBParallelComm::add_remote_proc(MBEntityHandle ent,
                                            int *remote_procs,
                                            MBEntityHandle *remote_hs,
                                            int remote_proc,
                                            MBEntityHandle remote_handle) 
{
  int* ptr = std::find( remote_procs, remote_procs+MAX_SHARING_PROCS, -1 );
  const size_t n = ptr - remote_procs;
  ptr = std::lower_bound( remote_procs, remote_procs+n, remote_proc );
  const size_t i = ptr - remote_procs;
  
  MBErrorCode result;
  const int invalid_proc = -1;
  const MBEntityHandle invalid_handle = 0;
  if (i == n || remote_procs[i] != remote_proc) {
    if (0 == n) {
      remote_procs[0] = remote_proc;
      remote_hs[0] = remote_handle;
    }
    else {
      insert_in_array( remote_procs, MAX_SHARING_PROCS, i, remote_proc );
      insert_in_array( remote_hs, MAX_SHARING_PROCS, i, remote_handle );
        // also insert this proc/handle if it's not already there
      ptr = std::lower_bound( remote_procs, remote_procs+n+1, procConfig.proc_rank() );
      const size_t i = ptr - remote_procs;
      if (i == n+1 || remote_procs[i] != (int)procConfig.proc_rank()) {
        insert_in_array( remote_procs, MAX_SHARING_PROCS, i, (int)procConfig.proc_rank());
        insert_in_array( remote_hs, MAX_SHARING_PROCS, i, ent);
      }
    }
    
    switch (n) {
      case 0:
        result = mbImpl->tag_set_data( sharedp_tag(), &ent, 1, remote_procs );
        RRA("Couldn't set sharedp tag");
        result = mbImpl->tag_set_data( sharedh_tag(), &ent, 1, remote_hs );
        RRA("Couldn't set sharedh tag");
        break;
      case 1:
          // going from 1 -> many, so clear single-value tag
        result = mbImpl->tag_set_data(  sharedp_tag(), &ent, 1, &invalid_proc );
        RRA("Couldn't set sharedp tag");
        result = mbImpl->tag_set_data( sharedh_tag(), &ent, 1, &invalid_handle );
        RRA("Couldn't set sharedh tag");
          // NO BREAK: fall through to next block to set many-valued tags
      default:
        result = mbImpl->tag_set_data(  sharedps_tag(), &ent, 1, remote_procs );
        RRA("Couldn't set sharedps tag");
        result = mbImpl->tag_set_data( sharedhs_tag(), &ent, 1, remote_hs );
        RRA("Couldn't set sharedhs tag");
        break;
    }
  }
  else if (remote_hs[i] != remote_handle) {
    assert(remote_hs[i] == invalid_handle);
    remote_hs[i] = remote_handle;
    if (n == 1) {
      result = mbImpl->tag_set_data( sharedh_tag(), &ent, 1, remote_hs );
      RRA("Couldn't set sharedh tag");
    }
    else {
      result = mbImpl->tag_set_data( sharedhs_tag(), &ent, 1, remote_hs );
      RRA("Couldn't set sharedhs tag");
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_range_map(MBRange &key_range, MBEntityHandle val_start,
                                           HandleMap &handle_map) 
{
  for (MBRange::const_pair_iterator key_it = key_range.const_pair_begin(); 
       key_it != key_range.const_pair_end(); key_it++) {
    int tmp_num = (*key_it).second - (*key_it).first + 1;
    handle_map.insert((*key_it).first, val_start, tmp_num);
    val_start += tmp_num;
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_sets(MBRange &entities,
                                      std::vector<unsigned char> &buff,
                                      unsigned char *&buff_ptr,
                                      const bool store_remote_handles,
                                      const int to_proc)
{
    // SETS:
    // . #sets
    // . for each set:
    //   - options[#sets] (unsigned int)
    //   - if (unordered) set range 
    //   - else if ordered
    //     . #ents in set
    //     . handles[#ents]
    //   - #parents
    //   - if (#parents) handles[#parents]
    //   - #children
    //   - if (#children) handles[#children]
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  MBErrorCode result;
  MBRange all_sets = entities.subset_by_type(MBENTITYSET);

  int buff_size = estimate_sets_buffer_size(all_sets, store_remote_handles);
  CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);

    // number of sets
  PACK_INT(buff_ptr, all_sets.size());

    // options for all sets
  std::vector<unsigned int> options(all_sets.size());
  MBRange::iterator rit;
  std::vector<MBEntityHandle> members;
  int i;
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
      result = mbImpl->get_meshset_options(*rit, options[i]);
      RRA("Failed to get meshset options.");
  }
  CHECK_BUFF_SPACE(buff, buff_ptr, all_sets.size()*sizeof(unsigned int));
  PACK_VOID(buff_ptr, &options[0], all_sets.size()*sizeof(unsigned int));
  
    // vectors/ranges
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
      MBRange set_range;
      if (options[i] & MESHSET_SET) {
        MBRange set_range;
        result = mbImpl->get_entities_by_handle(*rit, set_range);
        RRA("Failed to get set entities.");

        buff_size = RANGE_SIZE(set_range);
        CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);
        PACK_RANGE(buff_ptr, set_range);
      }
      else if (options[i] & MESHSET_ORDERED) {
        members.clear();
        result = mbImpl->get_entities_by_handle(*rit, members);
        RRA("Failed to get entities in ordered set.");
        
        CHECK_BUFF_SPACE(buff, buff_ptr,
                         members.size()*sizeof(MBEntityHandle)+sizeof(int));
        PACK_INT(buff_ptr, members.size());
        PACK_EH(buff_ptr, &members[0], members.size());
      }
  }
    // pack numbers of parents/children
  unsigned int tot_pch = 0;
  int num_pch;
  CHECK_BUFF_SPACE(buff, buff_ptr, 2*all_sets.size()*sizeof(int));
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
      // pack parents
    result = mbImpl->num_parent_meshsets(*rit, &num_pch);
    RRA("Failed to get num parents.");
    PACK_INT(buff_ptr, num_pch);
    tot_pch += num_pch;
    result = mbImpl->num_child_meshsets(*rit, &num_pch);
    RRA("Failed to get num children.");
    PACK_INT(buff_ptr, num_pch);
    tot_pch += num_pch;
  }

    // now pack actual parents/children
  members.clear();
  members.reserve(tot_pch);
  std::vector<MBEntityHandle> tmp_pch;
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
    result = mbImpl->get_parent_meshsets(*rit, tmp_pch);
    RRA("Failed to get parents.");
    std::copy(tmp_pch.begin(), tmp_pch.end(), std::back_inserter(members));
    tmp_pch.clear();
    result = mbImpl->get_child_meshsets(*rit, tmp_pch);
    RRA("Failed to get children.");
    std::copy(tmp_pch.begin(), tmp_pch.end(), std::back_inserter(members));
    tmp_pch.clear();
  }
  assert(members.size() == tot_pch);
  if (!members.empty()) {
    result = get_remote_handles(store_remote_handles,
                                &members[0], &members[0], 
                                members.size(), to_proc,
                                entities);
    RRA("Trouble getting remote handles for set parent/child sets.");
#ifndef NDEBUG
      // check that all handles are either sets or maxtype
    for (unsigned int __j = 0; __j < members.size(); __j++)
      assert((TYPE_FROM_HANDLE(members[__j]) == MBMAXTYPE &&
              ID_FROM_HANDLE(members[__j]) < (int)entities.size()) ||
             TYPE_FROM_HANDLE(members[__j]) == MBENTITYSET);
#endif        
    CHECK_BUFF_SPACE(buff, buff_ptr, members.size()*sizeof(MBEntityHandle));
    PACK_EH(buff_ptr, &members[0], members.size());
  }
    
    // pack the handles
  if (store_remote_handles && !all_sets.empty()) {
    buff_size = RANGE_SIZE(all_sets);
    CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);
    PACK_RANGE(buff_ptr, all_sets);
  }
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done packing sets." << std::endl;
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_sets(unsigned char *&buff_ptr,
                                        MBRange &entities,
                                        const bool store_remote_handles,
                                        const int from_proc)
{
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  MBErrorCode result;

  MBRange new_sets;
  int num_sets;
  UNPACK_INT(buff_ptr, num_sets);

  if (!num_sets) return MB_SUCCESS;
         
  std::vector<MBEntityHandle> members;
  int num_ents;
  std::vector<unsigned int> options_vec(num_sets);
      // option value
  if (num_sets)
    UNPACK_VOID(buff_ptr, &options_vec[0], num_sets*sizeof(unsigned int));

    // create sets
  int i;
  MBRange::const_iterator rit;
  for (i = 0; i < num_sets; i++) {
    
      // create the set
    MBEntityHandle set_handle;
    result = mbImpl->create_meshset(options_vec[i], set_handle);
    RRA("Failed to create set in unpack.");

    // make sure new sets handles are monotonically increasing
    assert(set_handle > *new_sets.rbegin());

    new_sets.insert(set_handle);
  }

  entities.merge(new_sets);
  
  for (rit = new_sets.begin(), i = 0; rit != new_sets.end(); rit++, i++) {
    if (options_vec[i] & MESHSET_SET) {
        // unpack entities as a range
      MBRange set_range, tmp_range;
      UNPACK_RANGE(buff_ptr, tmp_range);
      result = get_local_handles(tmp_range, set_range, entities);      
      RRA("Failed to get local handles for unordered set contents.");
      result = mbImpl->add_entities(*rit, set_range);
      RRA("Failed to add ents to unordered set in unpack.");
    }
    else if (options_vec[i] & MESHSET_ORDERED) {
        // unpack entities as vector, with length
      UNPACK_INT(buff_ptr, num_ents);
      members.resize(num_ents);
      if (num_ents) UNPACK_EH(buff_ptr, &members[0], num_ents);
      result = get_local_handles(&members[0], num_ents, entities);
      RRA("Failed to get local handles for ordered set contents.");
      result = mbImpl->add_entities(*rit, &members[0], num_ents);
      RRA("Failed to add ents to ordered set in unpack.");
    }
  }

  std::vector<int> num_pch(2*new_sets.size());
  std::vector<int>::iterator vit;
  int tot_pch = 0;
  for (vit = num_pch.begin(); vit != num_pch.end(); vit++) {
    UNPACK_INT(buff_ptr, *vit);
    tot_pch += *vit;
  }
  
  members.resize(tot_pch);
  UNPACK_EH(buff_ptr, &members[0], tot_pch);
  result = get_local_handles(&members[0], tot_pch, entities);
  RRA("Couldn't get local handle for parent/child sets.");

  int num = 0;
  MBEntityHandle *mem_ptr = &members[0];
  for (rit = new_sets.begin(); rit != new_sets.end(); rit++) {
      // unpack parents/children
    int num_par = num_pch[num++], num_child = num_pch[num++];
    if (num_par+num_child) {
      for (i = 0; i < num_par; i++) {
        assert(0 != mem_ptr[i]);
        result = mbImpl->add_parent_meshset(*rit, mem_ptr[i]);
        RRA("Failed to add parent to set in unpack.");
      }
      mem_ptr += num_par;
      for (i = 0; i < num_child; i++) {
        assert(0 != mem_ptr[i]);
        result = mbImpl->add_child_meshset(*rit, mem_ptr[i]);
        RRA("Failed to add child to set in unpack.");
      }
      mem_ptr += num_child;
    }
  }

    // unpack source handles
  MBRange dum_range;
  if (store_remote_handles && !new_sets.empty()) {
    UNPACK_RANGE(buff_ptr, dum_range);
    result = set_remote_data(new_sets, dum_range, from_proc);
    RRA("Couldn't set sharing data for sets");
  }

#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking sets." << std::endl;
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_adjacencies(MBRange &entities,
                                             MBRange::const_iterator &start_rit,
                                             MBRange &whole_range,
                                             unsigned char *&buff_ptr,
                                             int &count,
                                             const bool just_count,
                                             const bool store_handles,
                                             const int to_proc)
{
  return MB_FAILURE;
}

MBErrorCode MBParallelComm::unpack_adjacencies(unsigned char *&buff_ptr,
                                               MBRange &entities,
                                               const bool store_handles,
                                               const int from_proc)
{
  return MB_FAILURE;
}

MBErrorCode MBParallelComm::pack_tags(MBRange &entities,
                                      const std::vector<MBTag> &all_tags,
                                      const std::vector<MBRange> &tag_ranges,
                                      std::vector<unsigned char> &buff,
                                      unsigned char *&buff_ptr,
                                      const bool store_remote_handles,
                                      const int to_proc)
{
  

  MBErrorCode result;
  std::vector<MBTag>::const_iterator tag_it;
  std::vector<MBRange>::const_iterator rit;
  int count = 0;
  
  for (tag_it = all_tags.begin(), rit = tag_ranges.begin(); 
       tag_it != all_tags.end(); tag_it++, rit++) {

    result = packed_tag_size( *tag_it, *rit, count );
    if (MB_SUCCESS != result)
      return result;
  }
    
    // number of tags
  count += sizeof(int);

  CHECK_BUFF_SPACE(buff, buff_ptr, count);
  
  PACK_INT(buff_ptr, all_tags.size());
    
  for (tag_it = all_tags.begin(), rit = tag_ranges.begin(); 
       tag_it != all_tags.end(); tag_it++, rit++) {
    
    result = pack_tag( *tag_it, *tag_it, *rit, entities, buff, buff_ptr, 
                       store_remote_handles, to_proc );
    if (MB_SUCCESS != result)
      return result;
  }
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done packing tags." << std::endl;
#endif

  return MB_SUCCESS;
}
         

MBErrorCode MBParallelComm::packed_tag_size( MBTag tag,
                                             const MBRange &tagged_entities,
                                             int &count )
{
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  std::vector<int> var_len_sizes;
  std::vector<const void*> var_len_values;
    
  const TagInfo *tinfo = tagServer->get_tag_info(tag);
    // default value
  count += sizeof(int);
  if (NULL != tinfo->default_value()) 
    count += tinfo->default_value_size();

    // size, type, data type
  count += 3*sizeof(int);

    // name
  count += sizeof(int);
  count += tinfo->get_name().size();

    // range of tag
  count += sizeof(int) + tagged_entities.size() * sizeof(MBEntityHandle);

  if (tinfo->get_size() == MB_VARIABLE_LENGTH) {
    const int num_ent = tagged_entities.size();
      // send a tag size for each entity
    count += num_ent * sizeof(int);
      // send tag data for each entity
    var_len_sizes.resize( num_ent );
    var_len_values.resize( num_ent );
    MBErrorCode result = tagServer->get_data( tag,
                                              tagged_entities, 
                                              &var_len_values[0], 
                                              &var_len_sizes[0] );
    RRA("Failed to get lenghts of variable-length tag values.");
    count += std::accumulate( var_len_sizes.begin(), var_len_sizes.end(), 0 );
  }
  else {
      // tag data values for range or vector
    count += tagged_entities.size() * tinfo->get_size();
  }
  
  return MB_SUCCESS;
}


MBErrorCode MBParallelComm::pack_tag( MBTag src_tag,
                                      MBTag dst_tag,
                                      const MBRange &tagged_entities,
                                      const MBRange &whole_range,
                                      std::vector<unsigned char> &buff,
                                      unsigned char *&buff_ptr,
                                      const bool store_remote_handles,
                                      const int to_proc )
{
  MBErrorCode result;
  std::vector<int> var_len_sizes;
  std::vector<const void*> var_len_values;

  const TagInfo* tinfo = tagServer->get_tag_info(src_tag);
  if (!tinfo)
    return MB_TAG_NOT_FOUND;
    
  const TagInfo* dst_tinfo;
  if (src_tag == dst_tag) {
    dst_tinfo = tinfo;
  }
  else {
    dst_tinfo = tagServer->get_tag_info(dst_tag);
    if (!dst_tinfo)
      return MB_TAG_NOT_FOUND;
    if (dst_tinfo->get_size() != tinfo->get_size())
      return MB_TYPE_OUT_OF_RANGE;
    if (dst_tinfo->get_data_type() != tinfo->get_data_type() && 
        dst_tinfo->get_data_type() != MB_TYPE_OPAQUE &&
            tinfo->get_data_type() != MB_TYPE_OPAQUE)
      return MB_TYPE_OUT_OF_RANGE;
  }
    
    // size, type, data type
  CHECK_BUFF_SPACE(buff, buff_ptr, 3*sizeof(int));
  PACK_INT(buff_ptr, tinfo->get_size());
  MBTagType this_type;
  result = mbImpl->tag_get_type(dst_tag, this_type);
  PACK_INT(buff_ptr, (int)this_type);
  PACK_INT(buff_ptr, (int)(tinfo->get_data_type()));

    // default value
  if (NULL == tinfo->default_value()) {
    CHECK_BUFF_SPACE(buff, buff_ptr, sizeof(int));
    PACK_INT(buff_ptr, 0);
  }
  else {
    CHECK_BUFF_SPACE(buff, buff_ptr, tinfo->default_value_size());
    PACK_BYTES(buff_ptr, tinfo->default_value(), tinfo->default_value_size());
  }

    // name
  CHECK_BUFF_SPACE(buff, buff_ptr, tinfo->get_name().size());
  PACK_BYTES(buff_ptr, dst_tinfo->get_name().c_str(), dst_tinfo->get_name().size());

#ifdef DEBUG_PACKING
  std::cerr << "Packing tag \"" << tinfo->get_name() << "\"";
  if (tinfo != dst_tinfo)
    std::cerr << " (as tag \"" << dst_tinfo->get_name() << "\")";
  std::cerr << std::endl;
#endif    
    // pack entities
  CHECK_BUFF_SPACE(buff, buff_ptr, tagged_entities.size()*sizeof(MBEntityHandle)+sizeof(int));
  PACK_INT(buff_ptr, tagged_entities.size());
  result = get_remote_handles(store_remote_handles,
                              tagged_entities, (MBEntityHandle*)buff_ptr, to_proc,
                              whole_range);
#ifdef DEBUG_PACKING
  if (MB_SUCCESS != result) {
    std::cerr << "Trouble getting remote handles for tagged entities:" << std::endl;
    tagged_entities.print("  ");
  }
#else
  RRA("Trouble getting remote handles for tagged entities.");
#endif

  buff_ptr += tagged_entities.size() * sizeof(MBEntityHandle);

  const size_t num_ent = tagged_entities.size();
  if (tinfo->get_size() == MB_VARIABLE_LENGTH) {
    var_len_sizes.resize( num_ent, 0 );
    var_len_values.resize( num_ent, 0 );
    result = mbImpl->tag_get_data(src_tag, tagged_entities, &var_len_values[0], 
                                  &var_len_sizes[0] );
    RRA("Failed to get variable-length tag data in pack_tags.");
    CHECK_BUFF_SPACE(buff, buff_ptr, num_ent*sizeof(int));
    PACK_INTS(buff_ptr, &var_len_sizes[0], num_ent);
    for (unsigned int i = 0; i < num_ent; ++i) {
      CHECK_BUFF_SPACE(buff, buff_ptr, var_len_sizes[i]);
      PACK_VOID(buff_ptr, var_len_values[i], var_len_sizes[i]);
    }
  }
  else {
    CHECK_BUFF_SPACE(buff, buff_ptr, num_ent * tinfo->get_size());
    result = mbImpl->tag_get_data(src_tag, tagged_entities, buff_ptr);
    RRA("Failed to get tag data in pack_tags.");
    buff_ptr += num_ent * tinfo->get_size();
    PC(num_ent*tinfo->get_size(), " void");
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_tag_send_list( const MBRange& whole_range,
                                               std::vector<MBTag>& all_tags,
                                               std::vector<MBRange>& tag_ranges )
{
  std::vector<MBTag> tmp_tags;
  MBErrorCode result = tagServer->get_tags(tmp_tags);
  RRA("Failed to get tags in pack_tags.");

  std::vector<MBTag>::iterator tag_it;
  for (tag_it = tmp_tags.begin(); tag_it != tmp_tags.end(); tag_it++) {
    std::string tag_name;
    result = mbImpl->tag_get_name(*tag_it, tag_name);
    if (tag_name.c_str()[0] == '_' && tag_name.c_str()[1] == '_')
      continue;

    MBRange tmp_range;
    result = tagServer->get_entities(*tag_it, tmp_range);
    RRA("Failed to get entities for tag in pack_tags.");
    tmp_range = tmp_range.intersect(whole_range);

    if (tmp_range.empty()) continue;
        
      // ok, we'll be sending this tag
    all_tags.push_back( *tag_it );
    tag_ranges.push_back( MBRange() );
    tag_ranges.back().swap( tmp_range );
  }
  
  return MB_SUCCESS;
}



MBErrorCode MBParallelComm::unpack_tags(unsigned char *&buff_ptr,
                                        MBRange &entities,
                                        const bool store_remote_handles,
                                        const int from_proc)
{
    // tags
    // get all the tags
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  MBErrorCode result;
  
  int num_tags;
  UNPACK_INT(buff_ptr, num_tags);
  std::vector<MBEntityHandle> tag_ents;
  std::vector<const void*> var_len_vals;
  std::vector<int> var_lengths;

  for (int i = 0; i < num_tags; i++) {
    
        // tag handle
    MBTag tag_handle;

      // size, data type
    int tag_size, tag_data_type, tag_type;
    UNPACK_INT(buff_ptr, tag_size);
    UNPACK_INT(buff_ptr, tag_type);
    UNPACK_INT(buff_ptr, tag_data_type);
      
      // default value
    int def_val_size;
    UNPACK_INT(buff_ptr, def_val_size);
    void *def_val_ptr = NULL;
    if (def_val_size) {
      def_val_ptr = buff_ptr;
      buff_ptr += def_val_size;
      UPC(tag_size, " void");
    }
    
      // name
    int name_len;
    UNPACK_INT(buff_ptr, name_len);
    std::string tag_name( reinterpret_cast<char*>(buff_ptr), name_len );
    buff_ptr += name_len;
    UPC(64, " chars");
#ifdef DEBUG_PACKING
    std::cerr << "Unpacking tag " << tag_name << std::endl;
#endif    

      // create the tag
    if (tag_size == MB_VARIABLE_LENGTH) 
      result = mbImpl->tag_create_variable_length( tag_name.c_str(), (MBTagType)tag_type,
                                                   (MBDataType)tag_data_type, tag_handle,
                                                   def_val_ptr, def_val_size );
    else
      result = mbImpl->tag_create(tag_name.c_str(), tag_size, (MBTagType) tag_type, 
                                  (MBDataType) tag_data_type, tag_handle,
                                  def_val_ptr);
    if (MB_ALREADY_ALLOCATED == result) {
        // already allocated tag, check to make sure it's the same size, type, etc.
      const TagInfo *tag_info = tagServer->get_tag_info(tag_name.c_str());
      MBTagType this_type;
      result = mbImpl->tag_get_type(tag_handle, this_type);
      if (tag_size != tag_info->get_size() ||
          tag_type != this_type ||
          tag_data_type != tag_info->get_data_type() ||
          (def_val_ptr && !tag_info->default_value()) ||
          (!def_val_ptr && tag_info->default_value())) {
        RRA("Didn't get correct tag info when unpacking tag.");
      }
    }
    else if (MB_SUCCESS != result) return result;

      // go through handle vec (in buffer) and convert to local handles in-place
    int num_ents;
    UNPACK_INT(buff_ptr, num_ents);
    MBEntityHandle *handle_vec = (MBEntityHandle*)buff_ptr;
    result = get_local_handles(handle_vec, num_ents, entities);
    RRA("Failed to get local handles for tagged entities.");
    buff_ptr += num_ents * sizeof(MBEntityHandle);

      // if it's a handle type, also convert tag vals in-place in buffer
    if (MB_TYPE_HANDLE == tag_type) {
      MBEntityHandle *val_vec = (MBEntityHandle*)buff_ptr;
      result = get_local_handles(val_vec, num_ents, entities);
      RRA("Failed to get local handles for tag vals.");
    }

    if (tag_size == MB_VARIABLE_LENGTH) {
        // Be careful of alignment here.  If the integers are aligned
        // in the buffer, we can use them directly.  Otherwise we must
        // copy them.
      const int* size_arr;
      if (((size_t)buff_ptr)%4) {
        var_lengths.resize( num_ents );
        memcpy( &var_lengths[0], buff_ptr, num_ents*sizeof(int) );
        size_arr = &var_lengths[0];
      }
      else {
        size_arr = reinterpret_cast<const int*>(buff_ptr);
      }
      buff_ptr += sizeof(int) * num_ents;
      UPC(sizeof(int) * num_ents, " void");
      
        // get pointers into buffer for each tag value
      var_len_vals.resize(num_ents);
      for (std::vector<MBEntityHandle>::size_type i = 0; 
           i < (std::vector<MBEntityHandle>::size_type) num_ents; ++i) {
        var_len_vals[i] = buff_ptr;
        buff_ptr += size_arr[i];
        UPC(size_arr[i], " void");
      }
      result = mbImpl->tag_set_data( tag_handle, handle_vec, num_ents,
                                     &var_len_vals[0], size_arr );
      RRA("Trouble setting tag data when unpacking variable-length tag.");
    }
    else {
      result = mbImpl->tag_set_data(tag_handle, handle_vec,
                                    num_ents, buff_ptr);
      RRA("Trouble setting range-based tag data when unpacking tag.");
      buff_ptr += num_ents * tag_size;
      UPC(num_ents * tag_size, " void");
    }
  }
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking tags." << std::endl;
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::resolve_shared_ents(MBEntityHandle this_set,
                                                int resolve_dim,
                                                int shared_dim) 
{
  MBErrorCode result;
  MBRange proc_ents;
      // get the entities in the partition sets
  for (MBRange::iterator rit = partitionSets.begin(); rit != partitionSets.end(); rit++) {
    MBRange tmp_ents;
    result = mbImpl->get_entities_by_handle(*rit, tmp_ents, true);
    if (MB_SUCCESS != result) return result;
    proc_ents.merge(tmp_ents);
  }

    // resolve dim is maximal dim of entities in proc_ents
  if (-1 == resolve_dim) {
    resolve_dim = mbImpl->dimension_from_handle(*proc_ents.rbegin()); 
    RRA("Couldn't get dimension.");
    
  }

    // proc_ents should all be of same dimension
  if (resolve_dim > shared_dim &&
      mbImpl->dimension_from_handle(*proc_ents.rbegin()) !=
      mbImpl->dimension_from_handle(*proc_ents.begin())) {
    MBRange::iterator lower = proc_ents.lower_bound(MBCN::TypeDimensionMap[0].first),
      upper = proc_ents.upper_bound(MBCN::TypeDimensionMap[resolve_dim-1].second);
    proc_ents.erase(lower, upper);
  }
  
    // must call even if we don't have any entities, to make sure
    // collective comm'n works
  return resolve_shared_ents(this_set, proc_ents, resolve_dim, shared_dim);
}
  
MBErrorCode MBParallelComm::resolve_shared_ents(MBEntityHandle this_set,
                                                MBRange &proc_ents,
                                                int resolve_dim,
                                                int shared_dim) 
{
  MBErrorCode result;
  if (debug) std::cerr << "Resolving shared entities." << std::endl;

  if (-1 == shared_dim) {
    if (0 == resolve_dim) {
      result = mbImpl->get_dimension(shared_dim); 
      RRA("Couldn't get dimension.");
    }
    else shared_dim = mbImpl->dimension_from_handle(*proc_ents.begin())-1;
  }
  assert(shared_dim >= 0 && resolve_dim >= 0);
  
    // get the skin entities by dimension
  MBRange skin_ents[4];
  std::vector<int> gid_data;
  std::vector<MBEntityHandle> handle_vec;
  int skin_dim;

    // get the entities to be skinned
  if (resolve_dim < shared_dim) {
      // for vertex-based partition, it's the elements adj to the vertices
    result = mbImpl->get_adjacencies(proc_ents, shared_dim,
                                     false, skin_ents[resolve_dim],
                                     MBInterface::UNION);
    RRA("Failed getting skinned entities.");
    skin_dim = shared_dim-1;
  }
  else {
      // for element-based partition, it's just the elements
    skin_ents[resolve_dim] = proc_ents;
    skin_dim = resolve_dim-1;
  }

    // find the skin
  MBSkinner skinner(mbImpl);
  result = skinner.find_skin(skin_ents[skin_dim+1], skin_ents[skin_dim],
                             skin_ents[skin_dim], true);
  RRA("Failed to find skin.");
  if (debug) std::cerr << "Found skin, now resolving." << std::endl;

    // get entities adjacent to skin ents from shared_dim down to
    // zero; don't create them if they don't exist already
  for (int this_dim = skin_dim-1; this_dim >= 0; this_dim--) {
    result = mbImpl->get_adjacencies(skin_ents[skin_dim], this_dim,
                                     false, skin_ents[this_dim],
                                     MBInterface::UNION);
    RRA("Failed getting skin adjacencies.");
  }

    // resolve shared vertices first

    // global id tag
  MBTag gid_tag; int def_val = -1;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &def_val, true);
  if (MB_FAILURE == result) return result;

  else if (MB_ALREADY_ALLOCATED != result) {
      // just created it, so we need global ids
    result = assign_global_ids(0, skin_dim+1);
    RRA("Failed assigning global ids.");
  }

    // store index in temp tag; reuse gid_data 
  gid_data.resize(2*skin_ents[0].size());
  int idx = 0;
  for (MBRange::iterator rit = skin_ents[0].begin(); 
       rit != skin_ents[0].end(); rit++) 
    gid_data[idx] = idx, idx++;
  MBTag idx_tag;
  result = mbImpl->tag_create("__idx_tag", sizeof(int), MB_TAG_DENSE,
                              MB_TYPE_INTEGER, idx_tag, &def_val, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  result = mbImpl->tag_set_data(idx_tag, skin_ents[0], &gid_data[0]);
  RRA("Couldn't assign index tag.");

    // get gids for skin ents in a vector, to pass to gs
  result = mbImpl->tag_get_data(gid_tag, skin_ents[0], &gid_data[0]);
  RRA("Couldn't get gid tag for skin vertices.");

    // put handles in vector for passing to gs setup
  std::copy(skin_ents[0].begin(), skin_ents[0].end(), 
            std::back_inserter(handle_vec));
  
    // get a crystal router
  crystal_data *cd = procConfig.crystal_router();

/*  
    // get total number of entities; will overshoot highest global id, but
    // that's ok
  int num_total[2] = {0, 0}, num_local[2] = {0, 0};
  result = mbImpl->get_number_entities_by_dimension(0, 0, num_local);
  if (MB_SUCCESS != result) return result;
  int failure = MPI_Allreduce(num_local, num_total, 1,
                              MPI_INTEGER, MPI_SUM, procConfig.proc_comm());
  if (failure) {
    result = MB_FAILURE;
    RRA("Allreduce for total number of shared ents failed.");
  }
  
*/
    // call gather-scatter to get shared ids & procs
  gs_data *gsd;
  assert(sizeof(ulong_) == sizeof(MBEntityHandle));
  if (sizeof(int) != sizeof(ulong_)) {
    std::vector<long> lgid_data(gid_data.size());
    std::copy(gid_data.begin(), gid_data.end(), lgid_data.begin());
    gsd = gs_data_setup(skin_ents[0].size(), &lgid_data[0], 
                        (ulong_*)&handle_vec[0], 2, 1, 1, cd);
  }
  else {
    gsd = gs_data_setup(skin_ents[0].size(), (long*)&gid_data[0], 
                        (ulong_*)&handle_vec[0], 2, 1, 1, cd);
  }
  
  if (NULL == gsd) {
    result = MB_FAILURE;
    RRA("Couldn't create gs data.");
  }

    // get shared proc tags
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Couldn't get shared proc tags.");
  
    // load shared verts into a tuple, then sort by index
  tuple_list shared_verts;
  tuple_list_init_max(&shared_verts, 2, 0, 1, 0, 
                      skin_ents[0].size()*(MAX_SHARING_PROCS+1));
  unsigned int i = 0, j = 0;
  for (unsigned int p = 0; p < gsd->nlinfo->np; p++) 
    for (unsigned int np = 0; np < gsd->nlinfo->nshared[p]; np++) {
      shared_verts.vi[i++] = gsd->nlinfo->sh_ind[j];
      shared_verts.vi[i++] = gsd->nlinfo->target[p];
      shared_verts.vul[j] = gsd->nlinfo->ulabels[j];
      j++;
      shared_verts.n++;
    }
  
  int max_size = skin_ents[0].size()*(MAX_SHARING_PROCS+1);
  buffer sort_buffer;
  buffer_init(&sort_buffer, max_size);
  tuple_list_sort(&shared_verts, 0, &sort_buffer);
  buffer_free(&sort_buffer);

    // set sharing procs and handles tags on skin ents
  int maxp = -1;
  std::vector<int> sharing_procs(MAX_SHARING_PROCS);
  std::fill(sharing_procs.begin(), sharing_procs.end(), maxp);
  j = 0; i = 0;

    // get ents shared by 1 or n procs
  std::map<std::vector<int>, MBRange> proc_nranges;
  MBRange proc_verts;
  result = mbImpl->get_adjacencies(proc_ents, 0, false, proc_verts,
                                   MBInterface::UNION);
  RRA("Couldn't get proc_verts.");
  
  result = tag_shared_verts(shared_verts, skin_ents,
                            proc_nranges, proc_verts);
  RRA("Trouble tagging shared verts.");

    // get entities shared by 1 or n procs
  result = tag_shared_ents(resolve_dim, shared_dim, shared_verts, skin_ents,
                           proc_nranges);
  RRA("Trouble tagging shared entities.");

  tuple_list_free(&shared_verts);
  
  if (debug) {
    for (std::map<std::vector<int>, MBRange>::const_iterator mit = proc_nranges.begin();
         mit != proc_nranges.end(); mit++) {
      std::cout << "Iface: ";
      for (std::vector<int>::const_iterator vit = (mit->first).begin();
           vit != (mit->first).end(); vit++) std::cout << " " << *vit;
      std::cout << std::endl;
    }
  }
  
    // create the sets for each interface; store them as tags on
    // the interface instance
  MBRange iface_sets;
  result = create_interface_sets(proc_nranges, this_set, resolve_dim, shared_dim);
  RRA("Trouble creating iface sets.");

    // resolve shared entity remote handles; implemented in ghost cell exchange
    // code because it's so similar
  result = exchange_ghost_cells(-1, -1, 0, true, true);
  RRA("Trouble resolving shared entity remote handles.");

    // now set the shared/interface tag on non-vertex entities on interface
  result = tag_iface_entities();
  RRA("Failed to tag iface entities.");

    // now build parent/child links for interface sets
  result = create_iface_pc_links();
  RRA("Trouble creating interface parent/child links.");

  gs_data_free(gsd);
  
    // done
  return result;
}

MBErrorCode MBParallelComm::tag_iface_entities() 
{
  MBRange all_ents, if_ents;
  MBErrorCode result;
  for (MBRange::iterator if_it = interfaceSets.begin(); if_it != interfaceSets.end(); if_it++) {
      // get all ents
    result = mbImpl->get_entities_by_handle(*if_it, if_ents);
    RRA("Trouble getting iface entities.");
  }

  return set_pstatus_entities(if_ents, PSTATUS_INTERFACE, false, false);
}

MBErrorCode MBParallelComm::set_pstatus_entities(MBRange &pstatus_ents,
                                                 unsigned char pstatus_val,
                                                 bool lower_dim_ents,
                                                 bool verts_too,
                                                 int operation) 
{
  std::vector<unsigned char> pstatus_vals(pstatus_ents.size());
  MBRange all_ents, *range_ptr = &pstatus_ents;
  MBErrorCode result;
  if (lower_dim_ents || verts_too) {
    all_ents = pstatus_ents;
    range_ptr = &all_ents;
    int start_dim = (lower_dim_ents ? mbImpl->dimension_from_handle(*pstatus_ents.rbegin())-1 : 0);
    for (; start_dim >= 0; start_dim--) {
      result = mbImpl->get_adjacencies(all_ents, start_dim, true, all_ents,
                                       MBInterface::UNION);
      RRA(" ");
    }
  }
  if (MBInterface::UNION == operation) {
    result = mbImpl->tag_get_data(pstatus_tag(), *range_ptr, &pstatus_vals[0]);
    RRA("Couldn't get pstatus tag value.");
    for (unsigned int i = 0; i < pstatus_vals.size(); i++)
      pstatus_vals[i] |= pstatus_val;
  }
  else {
    for (unsigned int i = 0; i < pstatus_vals.size(); i++)
      pstatus_vals[i] = pstatus_val;
  }
  result = mbImpl->tag_set_data(pstatus_tag(), *range_ptr, &pstatus_vals[0]);
  RRA("Couldn't set pstatus tag value.");
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::set_pstatus_entities(MBEntityHandle *pstatus_ents,
                                                 int num_ents,
                                                 unsigned char pstatus_val,
                                                 bool lower_dim_ents,
                                                 bool verts_too,
                                                 int operation) 
{
  std::vector<unsigned char> pstatus_vals(num_ents);
  MBErrorCode result;
  if (lower_dim_ents || verts_too) {
      // in this case, call the range-based version
    MBRange tmp_range;
    std::copy(pstatus_ents, pstatus_ents+num_ents, mb_range_inserter(tmp_range));
    return set_pstatus_entities(tmp_range, pstatus_val, lower_dim_ents, 
                                verts_too, operation);
  }

  if (MBInterface::UNION == operation) {
    result = mbImpl->tag_get_data(pstatus_tag(), pstatus_ents, num_ents, &pstatus_vals[0]);
    RRA("Couldn't get pstatus tag value.");
    for (unsigned int i = 0; i < (unsigned int) num_ents; i++)
      pstatus_vals[i] |= pstatus_val;
  }
  else {
    for (unsigned int i = 0; i < (unsigned int) num_ents; i++)
      pstatus_vals[i] = pstatus_val;
  }
  result = mbImpl->tag_set_data(pstatus_tag(), pstatus_ents, num_ents, &pstatus_vals[0]);
  RRA("Couldn't set pstatus tag value.");
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::create_interface_sets(std::map<std::vector<int>, MBRange> &proc_nranges,
                                                  MBEntityHandle this_set,
                                                  int resolve_dim, int shared_dim) 
{
  if (proc_nranges.empty()) return MB_SUCCESS;
  
  int proc_ids[MAX_SHARING_PROCS];
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag,
                                            pstatus_tag);
  RRA("Trouble getting shared proc tags in create_interface_sets.");
  MBRange::iterator rit;

    // create interface sets, tag them, and tag their contents with iface set tag
  std::vector<MBEntityHandle> tag_vals;
  std::vector<unsigned char> pstatus;
  for (std::map<std::vector<int>,MBRange>::iterator mit = proc_nranges.begin();
       mit != proc_nranges.end(); mit++) {
      // create the set
    MBEntityHandle new_set;
    result = mbImpl->create_meshset(MESHSET_SET, new_set); 
    RRA("Failed to create interface set.");
    interfaceSets.insert(new_set);

      // add entities
    result = mbImpl->add_entities(new_set, mit->second); 
    RRA("Failed to add entities to interface set.");
      // tag set with the proc rank(s)
    if (mit->first.size() == 1)
      result = mbImpl->tag_set_data(sharedp_tag, &new_set, 1, 
                                    &(mit->first)[0]); 
    else {
      // pad tag data out to MAX_SHARING_PROCS with -1
      assert( mit->first.size() <= MAX_SHARING_PROCS );
      std::copy( mit->first.begin(), mit->first.end(), proc_ids );
      std::fill( proc_ids + mit->first.size(), proc_ids + MAX_SHARING_PROCS, -1 );
      result = mbImpl->tag_set_data(sharedps_tag, &new_set, 1, proc_ids );
    }
    RRA("Failed to tag interface set with procs.");
    
      // get the owning proc, then set the pstatus tag on iface set
    int min_proc = (mit->first)[0];
    unsigned char pval = (PSTATUS_SHARED | PSTATUS_INTERFACE);
    if (min_proc < (int) procConfig.proc_rank()) pval |= PSTATUS_NOT_OWNED;
    result = mbImpl->tag_set_data(pstatus_tag, &new_set, 1, &pval); 
    RRA("Failed to tag interface set with pstatus.");

      // tag the entities with the same thing
    pstatus.clear();
    pstatus.resize(mit->second.size(), pval);
    result = mbImpl->tag_set_data(pstatus_tag, mit->second, &pstatus[0]); 
    RRA("Failed to tag interface set entities with pstatus.");

      // finally, if I own this interface, add all the entities to my 
      // sharedEnts.ownedShared for the other procs
    if ((mit->first)[0] == (int)procConfig.proc_rank()) {
      std::vector<int>::const_iterator vit = (mit->first).begin();
      for (++vit; vit != (mit->first).end(); vit++) {
        int ind = get_buffers(*vit);
        assert(-1 != ind);
        sharedEnts2.merge(mit->second);
      }
    }
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::create_iface_pc_links() 
{
    // now that we've resolved the entities in the iface sets, 
    // set parent/child links between the iface sets

    // first tag all entities in the iface sets
  MBTag tmp_iface_tag;
  MBEntityHandle tmp_iface_set = 0;
  MBErrorCode result = mbImpl->tag_create("__tmp_iface", sizeof(MBEntityHandle),
                                          MB_TAG_DENSE, MB_TYPE_HANDLE,
                                          tmp_iface_tag, &tmp_iface_set);
  if (MB_ALREADY_ALLOCATED != result && MB_SUCCESS != result) 
    RRA("Failed to create temporary iface set tag.");

  MBRange iface_ents;
  std::vector<MBEntityHandle> tag_vals;
  MBRange::iterator rit;
  
  for (rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
      // tag entities with interface set
    iface_ents.clear();
    result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    RRA("Couldn't get entities in iface set.");
    
    if (iface_ents.empty()) continue;
    
    tag_vals.resize(iface_ents.size());
    std::fill(tag_vals.begin(), tag_vals.end(), *rit);
    result = mbImpl->tag_set_data(tmp_iface_tag, iface_ents, &tag_vals[0]); 
    RRA("Failed to tag iface entities with interface set.");
  }
  
    // now go back through interface sets and add parent/child links
  MBRange tmp_ents2;
  for (int d = 2; d >= 0; d--) {
    for (rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
        // get entities on this interface
      iface_ents.clear();
      result = mbImpl->get_entities_by_handle(*rit, iface_ents, true);
      RRA("Couldn't get entities by dimension.");
      if (iface_ents.empty() ||
          mbImpl->dimension_from_handle(*iface_ents.rbegin()) != d) continue;

        // get higher-dimensional entities and their interface sets
      result = mbImpl->get_adjacencies(&(*iface_ents.begin()), 1, d+1,
                                       false, tmp_ents2);
      RRA("Couldn't get adjacencies for interface sets.");
      tag_vals.resize(tmp_ents2.size());
      result = mbImpl->tag_get_data(tmp_iface_tag, tmp_ents2, &tag_vals[0]);
      RRA("Couldn't get iface set tag for interface sets.");
      
        // go through and for any on interface make it a parent
      MBEntityHandle last_set = 0;
      for (unsigned int i = 0; i < tag_vals.size(); i++) {
        if (tag_vals[i] && tag_vals[i] != last_set) {
          result = mbImpl->add_parent_child(tag_vals[i], *rit);
          RRA("Couldn't add parent/child link for interface set.");
          last_set = tag_vals[i];
        }
      }
    }
  }
  
    // delete the temporary tag
  result = mbImpl->tag_delete(tmp_iface_tag);
  RRA("Couldn't delete tmp iface tag.");

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::tag_shared_ents(int resolve_dim,
                                            int shared_dim,
                                            tuple_list &shared_verts,
                                            MBRange *skin_ents,
                                            std::map<std::vector<int>, MBRange> &proc_nranges) 
{
    // set sharing procs tags on other skin ents
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Trouble getting shared proc tags in tag_shared_ents.");
  const MBEntityHandle *connect; int num_connect;
  std::vector<int> sharing_procs, sharing_procs1, sharing_procs2;
  std::vector<int>::iterator vii;
  std::vector<unsigned char> pstatus_flags;
  std::vector<MBEntityHandle> dum_connect;

  for (int d = 3; d > 0; d--) {
    if (resolve_dim == d) continue;
    
    for (MBRange::iterator rit = skin_ents[d].begin();
         rit != skin_ents[d].end(); rit++) {
        // get connectivity
      result = mbImpl->get_connectivity(*rit, connect, num_connect, true,
                                        &dum_connect);
      RRA("Failed to get connectivity on non-vertex skin entities.");
 
        // if any vertices not shared, this entity isn't
      bool is_shared = true;
      pstatus_flags.resize( num_connect );
      result = mbImpl->tag_get_data(pstatus_tag, connect, num_connect,
                                    &pstatus_flags[0]);
      RRA("Couldn't get pstatus flag.");
      for (int nc = 0; nc < num_connect; nc++) {
        if (!(pstatus_flags[nc] & PSTATUS_SHARED)) {
          is_shared = false;
          break;
        }
      }
      if (!is_shared) continue;

      for (int nc = 0; nc < num_connect; nc++) {
        sharing_procs2.clear();
        
          // get sharing procs
        sharing_procs2.resize(1);
        result = mbImpl->tag_get_data(sharedp_tag, connect+nc, 1, &sharing_procs2[0]);
        RRA("Couldn't get sharedp_tag on skin vertices in entity.");
        if (sharing_procs2[0] == -1) {
          sharing_procs2.resize(MAX_SHARING_PROCS);
          result = mbImpl->tag_get_data(sharedps_tag, connect+nc, 1, &sharing_procs2[0]);
          RRA("Couldn't get sharedps_tag on skin vertices in entity.");
        }
        assert(-1 != sharing_procs2[0]);
          // remove any unnecessary entries
        vii = std::find( sharing_procs2.begin(), sharing_procs2.end(), -1 );
        sharing_procs2.erase( vii, sharing_procs2.end() );
        
          // build range of sharing procs for this vertex
          // intersect with range for this skin ent
        if (0 == nc) {
          sharing_procs.swap( sharing_procs2 );
        }
        else if (resolve_dim < shared_dim) {
          sharing_procs1.clear();
          set_union( sharing_procs.begin(), sharing_procs.end(), 
                     sharing_procs2.begin(), sharing_procs2.end(),
                     std::back_inserter( sharing_procs1 ) );
          sharing_procs.swap( sharing_procs1 );
        }
        else {
          sharing_procs1.clear();
          set_intersection( sharing_procs.begin(), sharing_procs.end(), 
                            sharing_procs2.begin(), sharing_procs2.end(),
                            std::back_inserter( sharing_procs1 ) );
          sharing_procs.swap( sharing_procs1 );
        }
      }

      if (sharing_procs.empty() && resolve_dim < shared_dim) continue;

        // intersection is the owning proc(s) for this skin ent
      if (sharing_procs.empty()) continue;

      proc_nranges[sharing_procs].insert(*rit);

      if (sharing_procs.size() < 2) {
        result = mbImpl->tag_set_data(sharedp_tag, &(*rit), 1,
                                      &sharing_procs[0]);
        RRA("Failed to set sharedp_tag on non-vertex skin entity.");
      }
      else {
          // fill extra entries with -1
        assert(sharing_procs.size() <= MAX_SHARING_PROCS);
        sharing_procs.resize( MAX_SHARING_PROCS, -1 );
        result = mbImpl->tag_set_data(sharedps_tag, &(*rit), 1,
                                      &sharing_procs[0]);
        RRA("Failed to set sharedps_tag on non-vertex skin entity.");
          // set hs tag as placeholder, will get filled in later
        static std::vector<MBEntityHandle> sharing_hs(MAX_SHARING_PROCS, 0);
        result = mbImpl->tag_set_data(sharedhs_tag, &(*rit), 1,
                                      &sharing_hs[0]);
        RRA("Failed to set sharedhs_tag on non-vertex skin entity.");
      }

        // reset sharing proc(s) tags
      sharing_procs.clear();
    }
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::tag_shared_verts(tuple_list &shared_ents,
                                             MBRange *skin_ents,
                                             std::map<std::vector<int>, MBRange> &proc_nranges,
                                             MBRange &proc_verts) 
{
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Trouble getting shared proc tags in tag_shared_verts.");
  
  unsigned int j = 0, i = 0;
  std::vector<int> sharing_procs, sharing_procs2;
  std::vector<MBEntityHandle> sharing_handles, sharing_handles2;
  
  while (j < 2*shared_ents.n) {
      // count & accumulate sharing procs
    int this_idx = shared_ents.vi[j];
    MBEntityHandle this_ent = skin_ents[0][this_idx];
    while (j < 2*shared_ents.n && shared_ents.vi[j] == this_idx) {
      j++;
      sharing_procs.push_back( shared_ents.vi[j++] );
      sharing_handles.push_back( shared_ents.vul[i++] );
    }

    if (sharing_procs.size() > 1) {
        // add current proc/handle to list
      sharing_procs.push_back(procConfig.proc_rank());
      sharing_handles.push_back(this_ent);
    }
      
      // sort sharing_procs and sharing_handles such that
      // sharing_procs is in ascending order.  Use temporary
      // lists and binary search to re-order sharing_handles.
    sharing_procs2 = sharing_procs;
    std::sort( sharing_procs2.begin(), sharing_procs2.end() );
    sharing_handles2.resize( sharing_handles.size() );
    for (size_t k = 0; k < sharing_handles.size(); ++k) {
      size_t idx = std::lower_bound( sharing_procs2.begin(), 
                                     sharing_procs2.end(), 
                                     sharing_procs[k] ) - sharing_procs2.begin();
      sharing_handles2[idx] = sharing_handles[k];
    }
    sharing_procs.swap( sharing_procs2 );
    sharing_handles.swap( sharing_handles2 );
    
    
    proc_nranges[sharing_procs].insert(this_ent);

    if (sharing_procs.size() == 1) {
      result = mbImpl->tag_set_data(sharedp_tag, &this_ent, 1,
                                    &sharing_procs[0]);
      result = mbImpl->tag_set_data(sharedh_tag, &this_ent, 1,
                                    &sharing_handles[0]);

    }
    else {
        // pad lists 
      assert( sharing_procs.size() <= MAX_SHARING_PROCS );
      sharing_procs.resize( MAX_SHARING_PROCS, -1 );
      sharing_handles.resize( MAX_SHARING_PROCS, 0 );
      result = mbImpl->tag_set_data(sharedps_tag, &this_ent, 1,
                                    &sharing_procs[0]);
      result = mbImpl->tag_set_data(sharedhs_tag, &this_ent, 1,
                                    &sharing_handles[0]);
    }
    RRA("Failed setting shared_procs tag on skin vertices.");

      // tag the entity as shared & interface too (this function only
      // ever called for interface vertices)
    unsigned char share_flag = PSTATUS_SHARED;
    result = mbImpl->tag_set_data(pstatus_tag, &this_ent, 1, &share_flag);
    RRA("Couldn't set shared tag on shared vertex.");

      // if this entity isn't owned here, add it to a list of ghosts
      // by owner
    if (sharing_procs[0] < (int)procConfig.proc_rank()) {
      for (unsigned int j = 0; j < sharing_procs.size() && sharing_procs[j] != -1; j++) {
        if (sharing_procs[j] == (int)procConfig.proc_rank()) continue;
        int ind = get_buffers(sharing_procs[j]);
        sharedEnts2.insert(this_ent);
      }
    }
    
      // otherwise it is owned here and shared with others
    else {
      for (unsigned int j = 0; j < sharing_procs.size() && sharing_procs[j] != -1; j++) {
        if (sharing_procs[j] == (int)procConfig.proc_rank()) continue;
        assert(sharing_procs[j] > (int)procConfig.proc_rank());
        int ind = get_buffers(sharing_procs[j]);
        sharedEnts2.insert(this_ent);
      }
    }

      // reset sharing proc(s) tags
    sharing_procs.clear();
    sharing_handles.clear();
  }

  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::get_iface_entities(int other_proc,
                                               MBRange &iface_ents,
                                               int dim) 
{
  MBRange iface_sets;
  MBErrorCode result = MB_SUCCESS;
  
  for (MBRange::iterator rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
    if (-1 != other_proc && !is_iface_proc(*rit, other_proc)) continue;
    
    if (-1 == dim) result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    else result = mbImpl->get_entities_by_dimension(*rit, dim, iface_ents);
    RRA(" Failed to get entities in iface set.");
  }
  
  return MB_SUCCESS;
}

  //! get processors with which this processor communicates; sets are sorted by processor
MBErrorCode MBParallelComm::get_interface_procs(std::set<unsigned int> &procs_set)
{
    // make sure the sharing procs vector is empty
  procs_set.clear();

    // pre-load vector of single-proc tag values
  unsigned int i, j;
  std::vector<int> iface_proc(interfaceSets.size());
  MBErrorCode result = mbImpl->tag_get_data(sharedp_tag(), interfaceSets, &iface_proc[0]);
  RRA("Failed to get iface_proc for iface sets.");

    // get sharing procs either from single-proc vector or by getting
    // multi-proc tag value
  int tmp_iface_procs[MAX_SHARING_PROCS];
  std::fill(tmp_iface_procs, tmp_iface_procs+MAX_SHARING_PROCS, -1);
  MBRange::iterator rit;
  for (rit = interfaceSets.begin(), i = 0; rit != interfaceSets.end(); rit++, i++) {
    if (-1 != iface_proc[i]) procs_set.insert((unsigned int) iface_proc[i]);
    else {
        // get the sharing_procs tag
      result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1,
                                    tmp_iface_procs);
      RRA("Failed to get iface_procs for iface set.");
      for (j = 0; j < MAX_SHARING_PROCS; j++) {
        if (-1 != tmp_iface_procs[j] && tmp_iface_procs[j] != (int)procConfig.proc_rank()) 
          procs_set.insert((unsigned int) tmp_iface_procs[j]);
        else {
          std::fill(tmp_iface_procs, tmp_iface_procs+j, -1);
          break;
        }
      }
    }
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::get_pstatus_entities(int dim,
                                                 unsigned char pstatus_val,
                                                 MBRange &pstatus_ents)
{
  MBRange ents;
  MBErrorCode result;
  
  if (-1 == dim) result = mbImpl->get_entities_by_handle(0, ents);
  else result = mbImpl->get_entities_by_dimension(0, dim, ents);
  RRA(" ");
  
  std::vector<unsigned char> pstatus(ents.size());
  result = mbImpl->tag_get_data(pstatus_tag(), ents, &pstatus[0]);
  RRA("Couldn't get pastatus tag.");
  MBRange::iterator rit = ents.begin();
  int i = 0;
  if (pstatus_val) {
    for (; rit != ents.end(); i++, rit++)
      if (pstatus[i]&pstatus_val &&
          (-1 == dim || mbImpl->dimension_from_handle(*rit) == dim)) 
        pstatus_ents.insert(*rit);
  }
  else {
    for (; rit != ents.end(); i++, rit++)
      if (!pstatus[i] &&
          (-1 == dim || mbImpl->dimension_from_handle(*rit) == dim)) 
        pstatus_ents.insert(*rit);
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::check_global_ids(MBEntityHandle this_set,
                                             const int dimension, 
                                             const int start_id,
                                             const bool largest_dim_only,
                                             const bool parallel)
{
    // global id tag
  MBTag gid_tag; int def_val = -1;
  MBErrorCode result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                                          MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                                          &def_val, true);
  if (MB_ALREADY_ALLOCATED != result &&
      MB_SUCCESS != result) {
    RRA("Failed to create/get gid tag handle.");
  }

  MBRange dum_range;
  if (MB_ALREADY_ALLOCATED == result) {
    void *tag_ptr = &def_val;
    MBErrorCode tmp_result = mbImpl->get_entities_by_type_and_tag(this_set, MBVERTEX, 
                                                                  &gid_tag, &tag_ptr, 1,
                                                                  dum_range);
    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
      RRA("Failed to get gid tag.");
    }
  }
  
  if (MB_ALREADY_ALLOCATED != result || !dum_range.empty()) {
      // just created it, so we need global ids
    result = assign_global_ids(this_set, dimension, start_id, largest_dim_only,
                               parallel);
    RRA("Failed assigning global ids.");
  }

  return MB_SUCCESS;
}

bool MBParallelComm::is_iface_proc(MBEntityHandle this_set,
                                   int to_proc) 
{
  int sharing_procs[MAX_SHARING_PROCS];
  std::fill(sharing_procs, sharing_procs+MAX_SHARING_PROCS, -1);
  MBErrorCode result = mbImpl->tag_get_data(sharedp_tag(), &this_set, 1,
                                            sharing_procs);
  if (to_proc == sharing_procs[0]) return true;
  
  result = mbImpl->tag_get_data(sharedps_tag(), &this_set, 1,
                                sharing_procs);
  for (int i = 0; i < MAX_SHARING_PROCS; i++) {
    if (to_proc == sharing_procs[i]) return true;
    else if (-1 == sharing_procs[i]) return false;
  }
  
  return false;
}

MBErrorCode MBParallelComm::filter_owned_shared( MBRange &ents,
                                                 bool owned_test,
                                                 bool owned_val,
                                                 bool shared_test,
                                                 bool shared_val,
                                                 int to_proc,
                                                 MBRange *returned_ents)
{
  if (!owned_test && !shared_test) return MB_FAILURE;

  MBRange tmp_ents;

    // Put into tmp_ents any entities which are not owned locally or
    // who are already shared with to_proc
  std::vector<unsigned char> shared_flags(ents.size());
  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), ents,
                                            &shared_flags[0]);
  RRA("Failed to get pstatus flag.");
  MBRange::const_iterator rit;
  int sharing_procs[MAX_SHARING_PROCS];
  std::fill(sharing_procs, sharing_procs+MAX_SHARING_PROCS, -1);
  int i;
  for (rit = ents.begin(), i = 0; rit != ents.end(); rit++, i++) {
    bool owned = !(PSTATUS_NOT_OWNED & shared_flags[i]),
        shared = (PSTATUS_SHARED & shared_flags[i]);

    bool owned_passed = !owned_test ||
        (owned_test && (owned_val == owned));
    bool shared_passed = !shared_test ||
        (shared_val == shared && (!shared_val || -1 == to_proc));
    
    if (owned_passed && shared_passed)
      tmp_ents.insert(*rit);
      
    else if (owned_passed && shared_test && -1 != to_proc &&
             shared_val == shared) {
      // we need to check sharing procs
      result = mbImpl->tag_get_data(sharedp_tag(), &(*rit), 1,
                                    sharing_procs);
      RRA(" ");
      if (-1 == sharing_procs[0]) {
        result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1,
                                      sharing_procs);
        assert(-1 != sharing_procs[0]);
        RRA(" ");
      }
      unsigned int j;
      for (j = 0; j < MAX_SHARING_PROCS; j++) {
          // if to_proc shares this entity, add it to list
        if (-1 != to_proc && sharing_procs[j] == to_proc) 
          tmp_ents.insert(*rit);
        
          // if we get here, no more sharing procs, and it's not shared
          // with to_proc
        else if (-1 == sharing_procs[j])
          break;
      }
      std::fill(sharing_procs, sharing_procs+j, -1);
    }
  }

  if (returned_ents)
    returned_ents->swap(tmp_ents);
  else
    ents.swap(tmp_ents);
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::populate_shared_ents() 
{
    // take sharedEnts2 range and separate into appropriate places in sharedEnts

  sharedEnts.clear();
  sharedEnts.resize(buffProcs.size());

  std::vector<int> sharedp(sharedEnts2.size());
  std::vector<MBEntityHandle> sharedh(sharedEnts2.size());
  MBErrorCode result = mbImpl->tag_get_data(sharedp_tag(), sharedEnts2, &sharedp[0]);
  RRA("Couldn't get sharedp tag.");
  result = mbImpl->tag_get_data(sharedh_tag(), sharedEnts2, &sharedh[0]);
  RRA("Couldn't get sharedh tag.");
  int sharedps[MAX_SHARING_PROCS];
  MBEntityHandle sharedhs[MAX_SHARING_PROCS];
  MBRange::iterator rit;
  unsigned int i, j;
  int ind;
  for (i = 0, rit = sharedEnts2.begin(); rit != sharedEnts2.end(); rit++, i++) {
    if (sharedp[i] != -1) {
        // just pairwise shared
      ind = get_buffers(sharedp[i]);
      assert(-1 != ind);
      sharedEnts[ind].remoteHandles.insert(sharedh[i]);
      if (sharedp[i] > (int)procConfig.proc_rank())
        sharedEnts[ind].ownedShared.insert(*rit);
      else
        sharedEnts[ind].localHandles.insert(*rit);
    }
    else {
      result = mbImpl->tag_get_data(sharedps_tag(), &*rit, 1, sharedps);
      RRA("Couldn't get sharedps tag.");
      result = mbImpl->tag_get_data(sharedhs_tag(), sharedEnts2, sharedhs);
      RRA("Couldn't get sharedhs tag.");
      for (j = 0; j < MAX_SHARING_PROCS && sharedps[j] != -1; j++) {
        ind = get_buffers(sharedps[j]);
        assert(-1 != ind);
        sharedEnts[ind].remoteHandles.insert(sharedhs[j]);
        if (sharedps[0] == (int)procConfig.proc_rank())
          sharedEnts[ind].ownedShared.insert(*rit);
        else
          sharedEnts[ind].localHandles.insert(*rit);
      }
      std::fill(sharedps, sharedps+j, -1);
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::exchange_ghost_cells(int ghost_dim, int bridge_dim,
                                                 int num_layers,
                                                 bool store_remote_handles,
                                                 bool wait_all)
{
    // if we're only finding out about existing ents, we have to be storing
    // remote handles too
  assert(num_layers > 0 || store_remote_handles);
  
    // re-populate the sharedEnts data structure
  MBErrorCode result = populate_shared_ents();
  RRA("Unable to populate shared ents during ghost exchange.");
  
    // get the b-dimensional interface(s) with with_proc, where b = bridge_dim
  
  int success;
  unsigned char *buff_ptr;

    // when this function is called, buffProcs should already have any 
    // communicating procs

    //===========================================
    // post ghost irecv's for ghost entities from all communicating procs
    //===========================================
    // index reqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(buffProcs.size(), MPI_REQUEST_NULL),
      send_reqs(buffProcs.size(), MPI_REQUEST_NULL);
  std::vector<unsigned int>::iterator proc_it;
  int ind;
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
      // for interface, only get ents from higher-rank procs
    if (0 == num_layers && procConfig.proc_rank() < *proc_it) continue;
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, buffProcs[ind],
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    //===========================================
    // get entities to be sent to neighbors
    //===========================================

    // done in a separate loop over procs because sometimes later procs 
    // need to add info to earlier procs' messages
  MBRange sent_ents[MAX_SHARING_PROCS];
  std::set<unsigned int> addl_procs[MAX_SHARING_PROCS], my_addl_procs;
  int num_higher = 0;
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {

    if (num_layers) {
      result = get_ghosted_entities(bridge_dim, ghost_dim, buffProcs[ind],
                                    num_layers, sent_ents[ind], addl_procs);
      RRA("Failed to get ghost layers.");

        // if we're exchanging ghosts, remove already-shared entities
      result = filter_owned_shared(sent_ents[ind], 
                                   false, false, true, false, *proc_it);
      RRA("Failed to filter already-shared entities in ghost exchange.");
    }
    else {
      result = get_iface_entities(buffProcs[ind], sent_ents[ind]);
      RRA("Failed to get interface layers.");

        // remove vertices here, for efficiency in the filter
      std::pair<MBRange::const_iterator,MBRange::const_iterator> vert_it =
          sent_ents[ind].equal_range(MBVERTEX);
      sent_ents[ind].erase(vert_it.first, vert_it.second);
    }

  }
  
    //===========================================
    // pack and send ents from this proc to others
    //===========================================
    // initialize sendReqs
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {

      // if resolving interface, only send to higher-rank procs
    if (0 == num_layers && procConfig.proc_rank() > *proc_it) continue;
    num_higher++;
    
      // get an estimate of the buffer size & pre-allocate buffer size
    unsigned int buff_size = estimate_ents_buffer_size(sent_ents[ind], 
                                                       store_remote_handles);
    RRA("Failed to estimate buffer size.");
    ownerSBuffs[ind].clear();
    ownerSBuffs[ind].reserve(buff_size);
    
      // buff_ptr points to the END (one past last occupied byte) of buffer
    buff_ptr = &ownerSBuffs[ind][0];

      // entities (including info about non-owned ents)
    result = pack_entities(sent_ents[ind], ownerSBuffs[ind], buff_ptr,
                           store_remote_handles, buffProcs[ind],
                           num_layers, sharedEnts[ind].ownedShared); 
    RRA("Packing entities failed.");

      // addl procs, only if not iface layer
    if (num_layers) {
      result = pack_addl_procs(addl_procs[ind], ownerSBuffs[ind], buff_ptr);
      RRA("Failed to pack addl procs.");
    }
  
      // now we're ready to send the buffer
    result = send_buffer(*proc_it, &ownerSBuffs[ind][0], 
                         buff_ptr-&ownerSBuffs[ind][0], MB_MESG_ENTS,
                         send_reqs[ind]);
    RRA("Failed to Isend in ghost exchange.");
  }

    //===========================================
    // receive/unpack new entities
    //===========================================
    // number of incoming messages for ghosts is the number of procs we 
    // communicate with; for iface, it's the number of those with lower rank
  int num_incoming = buffProcs.size();
  if (!num_layers) num_incoming -= num_higher;
  
  std::vector<MPI_Status> status(buffProcs.size());
  std::vector<std::vector<MBEntityHandle> > recd_ents(num_incoming);
  MBRange new_ghosts[MAX_SHARING_PROCS];
  
  while (num_incoming) {
      // wait for all recvs of ghost ents before proceeding,
      // b/c some procs may have sent to a 3rd proc ents owned by me;
    success = MPI_Waitany(buffProcs.size(), &recv_reqs[0], &ind, &status[0]);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
      // branch on message type
    if (MB_MESG_SIZE == status[0].MPI_TAG) {
        // incoming message just has size; resize buffer and re-call recv,
        // then re-increment incoming count
      int new_size = *((int*)&ghostRBuffs[ind][0]);
      assert(0 > new_size);
      result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                              MB_MESG_ENTS);
      RRA("Failed to resize recv buffer.");
      num_incoming++;
    }
    else if (MB_MESG_ENTS == status[0].MPI_TAG) {
      
        // incoming ghost entities; unpack; returns entities received
        // both from sending proc and from owning proc (which may be different)
      unsigned char *buff_ptr = &ghostRBuffs[ind][0];
      if (num_layers) {
        result = unpack_entities(buff_ptr,
                                 store_remote_handles,
                                 buffProcs[ind],
                                 ind, 
                                 new_ghosts[ind], new_ghosts,
                                 my_addl_procs);
        RRA("Failed to unpack entities.");
          // unpack addl procs I need to receive handles from
        result = unpack_addl_procs(buff_ptr, my_addl_procs);
        RRA("Failed to unpack addl procs.");
      }
      else {
        result = unpack_iface_entities(buff_ptr, buffProcs[ind], ind, recd_ents[ind]);
        RRA("Failed to unpack iface entities.");
      }
    }
  }

    // allocate buffs for any new entities, then add buffers 
    // for any new addl procs
  for (std::set<unsigned int>::iterator sit = my_addl_procs.begin();
       sit != my_addl_procs.end(); sit++) {
    ind = get_buffers(*sit);
    assert(-1 != ind);
  }
  recv_reqs.resize(buffProcs.size());
  send_reqs.resize(buffProcs.size());
    
    //===========================================
    // post recvs for remote handles of my sent ents
    //===========================================
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
      // skip if iface layer and lower-rank proc
    if (!num_layers && procConfig.proc_rank() > *proc_it) continue;
    
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, buffProcs[ind],
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }

    //===========================================
    // send local handles for new ghosts to owner, then add
    // those to ghost list for that owner
    //===========================================
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
      // skip if iface layer and higher-rank proc
    if (!num_layers && procConfig.proc_rank() < *proc_it) continue;
    buff_ptr = &ghostSBuffs[ind][0];
    if (num_layers)
      result = pack_remote_handles(new_ghosts[ind], *proc_it,
                                   ghostSBuffs[ind], buff_ptr);
    else
      result = pack_remote_handles(recd_ents[ind], *proc_it,
                                   ghostSBuffs[ind], buff_ptr);
    RRA("Failed to pack remote handles.");
    result = send_buffer(buffProcs[ind], &ghostSBuffs[ind][0], 
                         buff_ptr - &ghostSBuffs[ind][0], 
                         MB_MESG_REMOTE_HANDLES, send_reqs[ind]);
    RRA("Failed to send remote handles.");
  }
  
    //===========================================
    // process remote handles of my ghosteds
    //===========================================
  MBRange new_ghosted;
  num_incoming = (num_layers ? buffProcs.size() : num_higher);
  while (num_incoming) {
    success = MPI_Waitany(buffProcs.size(), &recv_reqs[0], &ind, &status[0]);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
      // branch on message type
    if (MB_MESG_SIZE == status[0].MPI_TAG) {
        // incoming message just has size; resize buffer and re-call recv,
        // then re-increment incoming count
      int new_size = *((int*)&ghostRBuffs[ind][0]);
      assert(0 > new_size);
      result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                              MB_MESG_REMOTE_HANDLES);
      RRA("Failed to resize recv buffer.");
      num_incoming++;
    }
    else if (MB_MESG_REMOTE_HANDLES == status[0].MPI_TAG) {
        // incoming remote handles
      result = unpack_remote_handles(buffProcs[ind], &ghostRBuffs[ind][0], num_layers, ind,
                                     new_ghosted);
      RRA("Failed to unpack remote handles.");

        // should receive back at least as many handles as entities we sent
        // (at least 'cuz other procs may have sent ents owned by me to another proc)
      assert(new_ghosted.size() >= sent_ents[ind].size());
    }
    else assert(false);
  }
    
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_sharedhps_tag(MBRange &new_ghosted, 
                                               MBRange &ghosted_ents, 
                                               unsigned int to_proc,
                                               std::vector<unsigned char> &buff, 
                                               unsigned char *&buff_ptr) 
{
    // only need to send entities which have neg sharedp tag
    // use buffer as tmp space to count how many
  MBRange tmp_range = new_ghosted.intersect(ghosted_ents);
  std::vector<int> shp_tag(tmp_range.size());
  MBErrorCode result = mbImpl->tag_get_data(sharedp_tag(), tmp_range,
                                            &shp_tag[0]);
  RRA("Failed to get sharedp_tag.");
  int num_neg = count_if(shp_tag.begin(), shp_tag.end(),
                         bind2nd(std::equal_to<int>(), -1));
  
  CHECK_BUFF_SPACE(buff, buff_ptr, sizeof(int) + num_neg*(sizeof(int)+sizeof(MBEntityHandle))*
                   MAX_SHARING_PROCS);
  
    // now pack the buffer
  PACK_INT(buff_ptr, num_neg);
  MBRange::iterator rit = tmp_range.begin();
  std::vector<int>::iterator vit = shp_tag.begin();
  for (; rit != tmp_range.end(); rit++, vit++) {
      // don't need to send for ents shared only with this proc
    if (-1 != *vit) continue;
    
      // pack directly into buffer, one past start to leave room for
      // # shared procs/handles
    int *tmp_sharedps = (int*) (buff_ptr + sizeof(int));
    result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1, tmp_sharedps);
    assert(MB_SUCCESS == result);
    
      // cast to int so we can compare to proc starting point
    int *tmp_sharedhs = (int*) std::find(tmp_sharedps, tmp_sharedps+MAX_SHARING_PROCS, -1);
    int num_pr = tmp_sharedhs-tmp_sharedps;
      // pack # shared handles/procs
    PACK_INT(buff_ptr, num_pr);
      // skip past shared procs space, then pack handles directly
    buff_ptr += num_pr*sizeof(int);
    result = mbImpl->tag_get_data(sharedhs_tag(), &(*rit), 1, buff_ptr);
    RRA("Failed to get sharedhs tag.");
    buff_ptr += num_pr * sizeof(MBEntityHandle);
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_sharedhps_tag(unsigned int from_proc,
                                                 unsigned char *&buff_ptr) 
{
    // get # entities
  int num_ents, num_procs;
  UNPACK_INT(buff_ptr, num_ents);
  int tmp_sharedps[MAX_SHARING_PROCS];
  MBEntityHandle tmp_sharedhs[MAX_SHARING_PROCS];
  num_procs = MAX_SHARING_PROCS;
  for (int i = 0; i < num_ents; i++) {
      // initialize tmp proc/handle space from last iterate
    std::fill(tmp_sharedps, tmp_sharedps+num_procs, -1);
    std::fill(tmp_sharedhs, tmp_sharedhs+num_procs, 0);

    UNPACK_INT(buff_ptr, num_procs);
    assert(num_procs <= MAX_SHARING_PROCS && num_procs >= 2);
    int *proc_ptr = (int*) buff_ptr;
    int *this_proc = std::find(proc_ptr, proc_ptr+num_procs, 
                               (int)procConfig.proc_rank());
    assert(this_proc - proc_ptr < num_procs);
    MBEntityHandle my_handle = (MBEntityHandle)(proc_ptr + num_procs)[this_proc-proc_ptr];
    memcpy(tmp_sharedps, proc_ptr, num_procs*sizeof(int));
    memcpy(tmp_sharedhs, proc_ptr+num_procs, num_procs*sizeof(MBEntityHandle));
    MBErrorCode result = mbImpl->tag_set_data(sharedps_tag(), &my_handle, 1, tmp_sharedps);
    RRA("Failed to set sharedps tag.");
    result = mbImpl->tag_set_data(sharedhs_tag(), &my_handle, 1, tmp_sharedhs);
    RRA("Failed to set sharedhs tag.");

    buff_ptr += num_procs*(sizeof(int) + sizeof(MBEntityHandle));
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_remote_handles(MBRange &entities,
                                                unsigned int to_proc,
                                                std::vector<unsigned char> &buff,
                                                unsigned char *&buff_ptr) 
{
    // 2 vectors of handles
  CHECK_BUFF_SPACE(buff, buff_ptr, sizeof(int) + 2*entities.size()*sizeof(MBEntityHandle));
  
  PACK_INT(buff_ptr, entities.size());
  MBEntityHandle *local_handles = (MBEntityHandle*)buff_ptr,
      *remote_handles = local_handles + entities.size();
  std::copy(entities.begin(), entities.end(), local_handles);
  MBRange dum_range;
  MBErrorCode result = get_remote_handles(true, local_handles, remote_handles, 
                                          entities.size(), to_proc, dum_range);
  RRA("Trouble setting remote data range on sent entities in ghost exchange.");

  buff_ptr += 2*sizeof(MBEntityHandle)*entities.size();
  
  return result;
}

MBErrorCode MBParallelComm::pack_remote_handles(std::vector<MBEntityHandle> &entities,
                                                unsigned int to_proc,
                                                std::vector<unsigned char> &buff,
                                                unsigned char *&buff_ptr) 
{
    // IMPLEMENTATION IDENTICAL TO OTHER VERSION; COULD PROBABLY DO WITH TEMPLATES...
    // 2 vectors of handles
  CHECK_BUFF_SPACE(buff, buff_ptr, sizeof(int) + 2*entities.size()*sizeof(MBEntityHandle));
  
  PACK_INT(buff_ptr, entities.size());
  MBEntityHandle *local_handles = (MBEntityHandle*)buff_ptr,
      *remote_handles = local_handles + entities.size();
  std::copy(entities.begin(), entities.end(), local_handles);
  MBRange dum_range;
  MBErrorCode result = get_remote_handles(true, local_handles, remote_handles, 
                                          entities.size(), to_proc, dum_range);
  RRA("Trouble setting remote data range on sent entities in ghost exchange.");

  buff_ptr += 2*sizeof(MBEntityHandle)*entities.size();
  
  return result;
}

MBErrorCode MBParallelComm::unpack_remote_handles(unsigned int from_proc,
                                                  unsigned char *&buff_ptr,
                                                  const int num_layers,
                                                  const int ind,
                                                  MBRange &new_ghosted) 
{
    // incoming remote handles; use to set remote handles
  int num_eh;
  UNPACK_INT(buff_ptr, num_eh);
  MBEntityHandle *remote_handles = (MBEntityHandle*)buff_ptr,
      *local_handles = remote_handles + num_eh;
  MBErrorCode result = set_remote_data(local_handles, remote_handles, num_eh, from_proc);
  RRA("Trouble setting remote data range on sent entities in ghost exchange.");

  std::copy(local_handles, local_handles+num_eh, mb_range_inserter(new_ghosted));
  std::copy(local_handles, local_handles+num_eh, 
            mb_range_inserter(sharedEnts[ind].localHandles));
  std::copy(remote_handles, remote_handles+num_eh, 
            mb_range_inserter(sharedEnts[ind].remoteHandles));

  buff_ptr += 2*num_eh*sizeof(MBEntityHandle);

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_addl_procs(MBRange &ghosted_ents, 
                                           unsigned int to_proc,
                                           std::set<unsigned int> *addl_procs)  
{
    // for all the procs owning entities ghosted here, if I'm sending those
    // ghosts to 3rd proc, tell the owner it needs to communicate with that 3rd
    // proc too
  std::vector<int> tmp_procs(ghosted_ents.size()),
      tmp_procss(MAX_SHARING_PROCS);
  MBErrorCode result = mbImpl->tag_get_data(sharedp_tag(), ghosted_ents, 
                                            &tmp_procs[0]);
  RRA("Couldn't get sharedp tag.");

  MBRange::iterator rit;
  unsigned int i;
  for (i = 0, rit = ghosted_ents.begin(); 
       rit != ghosted_ents.end(); rit++, i++) {
    if (-1 == tmp_procs[i]) {
      result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1,
                                    &tmp_procss[0]);
        // could be it's not yet shared
      if (MB_TAG_NOT_FOUND == result) continue;
      RRA("Couldn't get sharedps tag.");
    }
    else tmp_procss[0] = tmp_procs[i];
    
    if (tmp_procss[0] != (int)procConfig.proc_rank() &&
        tmp_procss[0] != (int) to_proc) 
      addl_procs[tmp_procss[0]].insert(to_proc);
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::get_ghosted_entities(int bridge_dim,
                                                 int ghost_dim,
                                                 int to_proc, 
                                                 int num_layers,
                                                 MBRange &ghosted_ents,
                                                 std::set<unsigned int> *addl_procs) 
{
    // get bridge ents on interface(s)
  MBRange from_ents;
  MBErrorCode result = MB_SUCCESS;
  assert(0 < num_layers);
  for (MBRange::iterator rit = interfaceSets.begin(); rit != interfaceSets.end();
       rit++) {
    if (!is_iface_proc(*rit, to_proc)) continue;
      
      // get starting "from" entities
    if (bridge_dim == -1)
      result = mbImpl->get_entities_by_handle(*rit, from_ents);
    else
      result = mbImpl->get_entities_by_dimension(*rit, bridge_dim, from_ents);
    RRA("Couldn't get bridge ents in the set.");

      // need to get layers of bridge-adj entities
    result = MeshTopoUtil(mbImpl).get_bridge_adjacencies(from_ents, bridge_dim,
                                                         ghost_dim, ghosted_ents, 
                                                         num_layers);
    RRA("Couldn't get bridge adjacencies.");
  }
  
  result = add_verts(ghosted_ents);
  RRA("Couldn't add verts.");

  result = filter_owned_shared(ghosted_ents,
                               false, false, true, false, to_proc);
  RRA("Trouble filtering.");

  result = get_addl_procs(ghosted_ents, to_proc, addl_procs);
  RRA("Failed to get addl procs.");
  return result;
}

MBErrorCode MBParallelComm::add_verts(MBRange &sent_ents) 
{
      // get the verts adj to these entities, since we'll have to send those too

    // first check sets
  std::pair<MBRange::const_iterator, MBRange::const_iterator>
      set_range = sent_ents.equal_range(MBENTITYSET);
  MBErrorCode result = MB_SUCCESS, tmp_result;
  for (MBRange::const_iterator rit = set_range.first; rit != set_range.second; rit++) {
    tmp_result = mbImpl->get_entities_by_type(*rit, MBVERTEX, sent_ents);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
  }
  RRA("Failed to get contained verts.");
  
    // now non-sets
  MBRange tmp_ents;
  std::copy(sent_ents.begin(), set_range.first, mb_range_inserter(tmp_ents));
  result = mbImpl->get_adjacencies(tmp_ents, 0, false, sent_ents,
                                   MBInterface::UNION);
  RRA("Couldn't get vertices adj to ghosted ents.");

  return result;
}


MBErrorCode MBParallelComm::exchange_tags(std::vector<MBTag> &tags)
{
  MBErrorCode result;
  int success;

    // get all procs interfacing to this proc
  std::set<unsigned int> exch_procs;
  result = get_comm_procs(exch_procs);  

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::vector<MPI_Status> gstatus(MAX_SHARING_PROCS);
  std::vector<unsigned int>::iterator sit;
  int ind;
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, *sit,
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    // pack and send tags from this proc to others
    // make sendReqs vector to simplify initialization
  std::fill(sendReqs, sendReqs+MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    
    MBRange tag_ents;
    
      // get bridge ents on interface(s)
    for (MBRange::iterator rit = interfaceSets.begin(); rit != interfaceSets.end();
         rit++) {
      if (!is_iface_proc(*rit, *sit)) continue;

      int owner;
      result = get_owner(*rit, owner);
      if (MB_SUCCESS != result || owner != (int)proc_config().proc_rank()) 
        continue;
      
      result = mbImpl->get_entities_by_handle(*rit, tag_ents);
      RRA("Failed to get tag ents for exchange.");
    }

      // also get ghosted entities for this proc
    if (!sharedEnts[ind].ownedShared.empty())
      tag_ents.merge(sharedEnts[ind].ownedShared);

      // pack-send; this also posts receives if store_remote_handles is true
    std::vector<MBRange> tag_ranges;
    for (std::vector<MBTag>::iterator vit = tags.begin(); vit != tags.end(); vit++) {
      const void* ptr;
      int size;
      if (tagServer->get_default_data_ref( *vit, ptr, size ) != MB_SUCCESS) {
        MBRange tagged_ents;
        tagServer->get_entities( *vit, tagged_ents );
        tag_ranges.push_back(tag_ents.intersect(tagged_ents));
      } 
      else {
        tag_ranges.push_back(tag_ents);
      }
    }
    
      // pack the data
    unsigned char *buff_ptr = &ownerSBuffs[ind][0];
    result = pack_tags(tag_ents,
                       tags, tag_ranges, 
                       ownerSBuffs[ind], buff_ptr, true, *sit);
    RRA("Failed to count buffer in pack_send_tag.");

      // now send it
    result = send_buffer(*sit, &ownerSBuffs[ind][0], 
                         buff_ptr-&ownerSBuffs[ind][0], 
                         MB_MESG_TAGS, sendReqs[ind]);
    RRA("Failed to send buffer.");
                         
  }
  
    // receive/unpack tags
  int num_incoming = exch_procs.size();
  
  while (num_incoming) {
    int ind;
    MPI_Status status;
    success = MPI_Waitany(MAX_SHARING_PROCS, &recv_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
    int new_size;
    unsigned char *buff_ptr;
    MBRange dum_range;
    
      // branch on message type
    switch (status.MPI_TAG) {
      case MB_MESG_SIZE:
          // incoming message just has size; resize buffer and re-call recv,
          // then re-increment incoming count
        assert(ind < MAX_SHARING_PROCS);
        new_size = *((int*)&ghostRBuffs[ind][0]);
        assert(0 > new_size);
        result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                                MB_MESG_TAGS);
        RRA("Failed to resize recv buffer.");
        num_incoming++;
        break;
      case MB_MESG_TAGS:
          // incoming ghost entities; process
          buff_ptr = &ghostRBuffs[ind][0];
          result = unpack_tags(buff_ptr, dum_range, true,
                               buffProcs[ind]);
        RRA("Failed to recv-unpack-tag message.");
        break;
      default:
        result = MB_FAILURE;
        RRA("Failed to get message of correct type in exch_tags.");
        break;
    }
  }
  
    // ok, now wait
  MPI_Status status[MAX_SHARING_PROCS];
  success = MPI_Waitall(MAX_SHARING_PROCS, &sendReqs[0], status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Failure in waitall in tag exchange.");
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::exchange_tags( MBTag src_tag, 
                                           MBTag dst_tag, 
                                           const MBRange& entities )
{
  MBErrorCode result;
  int success;

    // get all procs interfacing to this proc
  std::set<unsigned int> exch_procs;
  result = get_comm_procs(exch_procs);  

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::vector<MPI_Status> gstatus(MAX_SHARING_PROCS);
  std::vector<unsigned int>::iterator sit;
  int ind;
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, *sit,
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    // figure out which entities are shared with which processors
  std::map<int,MBRange> proc_ents;
  int other_procs[MAX_SHARING_PROCS], num_sharing;
  for (MBRange::const_iterator i = entities.begin(); i != entities.end(); ++i) {
    int owner;
    result = get_owner( *i, owner );
    RRA("Failed to get entity owner.");

      // only send entities that this proc owns
    if ((unsigned)owner != proc_config().proc_rank()) 
      continue;
    
    result = get_sharing_parts( *i, other_procs, num_sharing );
    RRA("Failed to get procs sharing entity.");
    if (num_sharing == 0) // keep track of non-shared entities for later
      proc_ents[proc_config().proc_rank()].insert( *i );
    for (int j = 0; j < num_sharing; ++j)
      proc_ents[other_procs[j]].insert( *i );
  }
  
    // pack and send tags from this proc to others
    // make sendReqs vector to simplify initialization
  std::fill(sendReqs, sendReqs+MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::map<unsigned int,MBRange>::const_iterator mit;
  
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    
      // count first
      // buffer needs to begin with the number of tags (one)
    int buff_size = sizeof(int);
    result = packed_tag_size( src_tag, proc_ents[*sit], buff_size );
    RRA("Failed to count buffer in pack_send_tag.");

    unsigned char *buff_ptr = &ownerSBuffs[ind][0];
    CHECK_BUFF_SPACE(ownerSBuffs[ind], buff_ptr, buff_size);
    PACK_INT( buff_ptr, 1 ); // number of tags
    result = pack_tag( src_tag, dst_tag, proc_ents[*sit], proc_ents[*sit],
                       ownerSBuffs[ind], buff_ptr, true, *sit );
    RRA("Failed to pack buffer in pack_send_tag.");

      // if the message is large, send a first message to tell how large
    if (INITIAL_BUFF_SIZE < buff_size) {
      int tmp_buff_size = -buff_size;
      int success = MPI_Send(&tmp_buff_size, sizeof(int), MPI_UNSIGNED_CHAR, 
                             *sit, MB_MESG_SIZE, procConfig.proc_comm());
      if (success != MPI_SUCCESS) return MB_FAILURE;
    }
    
      // send the buffer
    success = MPI_Isend(&ownerSBuffs[ind][0], buff_size, MPI_UNSIGNED_CHAR, *sit, 
                        MB_MESG_TAGS, procConfig.proc_comm(), &sendReqs[ind]);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
  
    // receive/unpack tags
  int num_incoming = exch_procs.size();
  
  while (num_incoming) {
    int ind;
    MPI_Status status;
    success = MPI_Waitany(MAX_SHARING_PROCS, &recv_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
    int new_size;
    unsigned char *buff_ptr;
    MBRange dum_range;
    
      // branch on message type
    switch (status.MPI_TAG) {
      case MB_MESG_SIZE:
          // incoming message just has size; resize buffer and re-call recv,
          // then re-increment incoming count
        assert(ind < MAX_SHARING_PROCS);
        new_size = *((int*)&ghostRBuffs[ind][0]);
        assert(0 > new_size);
        result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                                MB_MESG_TAGS);
        RRA("Failed to resize recv buffer.");
        num_incoming++;
        break;
      case MB_MESG_TAGS:
          // incoming ghost entities; process
          buff_ptr = &ghostRBuffs[ind][0];
          result = unpack_tags(buff_ptr, dum_range, true,
                               buffProcs[ind]);
        RRA("Failed to recv-unpack-tag message.");
        break;
      default:
        result = MB_FAILURE;
        RRA("Failed to get message of correct type in exch_tags.");
        break;
    }
  }
  
    // ok, now wait
  MPI_Status status[MAX_SHARING_PROCS];
  success = MPI_Waitall(MAX_SHARING_PROCS, &sendReqs[0], status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Failure in waitall in tag exchange.");
  }
  
    // if src and destination tags aren't the same, need to copy 
    // values for local entities
  if (src_tag != dst_tag) {
    const MBRange& myents = proc_ents[proc_config().proc_rank()];
    std::vector<const void*> data_ptrs(myents.size());
    std::vector<int> data_sizes(myents.size());
    result = get_moab()->tag_get_data( src_tag, myents, &data_ptrs[0], &data_sizes[0] );
    RRA("Failure to get pointers to local data.");
    result = get_moab()->tag_set_data( dst_tag, myents, &data_ptrs[0], &data_sizes[0] );
    RRA("Failure to get pointers to local data.");
  }  
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::update_shared_mesh()
{
  MBErrorCode result;
  int success;

    // ,,,
    /*

    // get all procs interfacing to this proc
  std::set<unsigned int> iface_procs;
  result = get_interface_procs(iface_procs);
  RRA("Failed to get iface sets, procs");

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(2*MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::vector<MPI_Status> gstatus(MAX_SHARING_PROCS);
  std::set<unsigned int>::iterator sit;
  for (sit = iface_procs.begin(); sit != iface_procs.end(); sit++) {
    int ind = get_buffers(*sit);
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, *sit,
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    // pack and send vertex coordinates from this proc to others
    // make sendReqs vector to simplify initialization
  std::fill(sendReqs, sendReqs+2*MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  MBRange recd_ents[MAX_SHARING_PROCS];
  
  for (sit = iface_procs.begin(); sit != iface_procs.end(); sit++) {
    int ind = get_buffers(*sit);
    
    MBRange vertices;
    for (MBRange::iterator rit = interfaceSets.begin(); rit != interfaceSets.end();
         rit++) {
      if (!is_iface_proc(*rit, *sit)) 
        continue;
      
      result = mbImpl->get_entities_by_type( *rit, MBVERTEX, vertices );
      RRA("Bad interface set.");
    }
    std::map<unsigned int,MBRange>::iterator ghosted = ghostedEnts.find(*sit);
    if (ghosted != ghostedEnts.end()) {
      MBRange::iterator e = ghosted->second.upper_bound(MBVERTEX);
      vertices.merge( ghosted->second.begin(), e );
    }

      // pack-send; this also posts receives if store_remote_handles is true
    MBRange sent;
    result = pack_send_entities(*sit, vertices, false, false, 
                                false, true,
                                ownerSBuffs[ind], ownerRBuffs[MAX_SHARING_PROCS+ind], 
                                sendReqs[ind], recv_reqs[MAX_SHARING_PROCS+ind], 
                                sent);
    RRA("Failed to pack-send in mesh update exchange.");
  }
  
    // receive/unpack entities
    // number of incoming messages depends on whether we're getting back
    // remote handles
  int num_incoming = iface_procs.size();
  
  while (num_incoming) {
    int ind;
    MPI_Status status;
    success = MPI_Waitany(2*MAX_SHARING_PROCS, &recv_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
    std::vector<MBEntityHandle> remote_handles_v, sent_ents_tmp;
    MBRange remote_handles_r;
    int new_size;
    
      // branch on message type
    switch (status.MPI_TAG) {
      case MB_MESG_SIZE:
          // incoming message just has size; resize buffer and re-call recv,
          // then re-increment incoming count
        assert(ind < MAX_SHARING_PROCS);
        new_size = *((int*)&ghostRBuffs[ind][0]);
        assert(0 > new_size);
        result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                                MB_MESG_ENTS);
        RRA("Failed to resize recv buffer.");
        num_incoming++;
        break;
      case MB_MESG_ENTS:
          // incoming ghost entities; process
        result = recv_unpack_entities(buffProcs[ind], true,
                                      false, 
                                      ghostRBuffs[ind], ghostSBuffs[ind], 
                                      sendReqs[ind], recd_ents[ind]);
        RRA("Failed to recv-unpack message.");
        break;
    }
  }
  
    // ok, now wait if requested
  MPI_Status status[2*MAX_SHARING_PROCS];
  success = MPI_Waitall(2*MAX_SHARING_PROCS, &sendReqs[0], status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Failure in waitall in ghost exchange.");
  }
  
  return MB_SUCCESS;
}
MBErrorCode MBParallelComm::update_iface_sets(MBRange &sent_ents,
                                              std::vector<MBEntityHandle> &remote_handles, 
                                              int from_proc) 
{
  std::vector<MBEntityHandle>::iterator remote_it = remote_handles.begin();
  MBRange::iterator sent_it = sent_ents.begin();
  MBRange ents_to_remove;
  for (; sent_it != sent_ents.end(); sent_it++, remote_it++) {
    if (!*remote_it) ents_to_remove.insert(*sent_it);
  }
  
  for (MBRange::iterator set_it = interfaceSets.begin(); set_it != interfaceSets.end(); set_it++) {
    if (!is_iface_proc(*set_it, from_proc)) continue;
    MBErrorCode result = mbImpl->remove_entities(*set_it, ents_to_remove);
    RRA("Couldn't remove entities from iface set in update_iface_sets.");
  }

*/
  
  return MB_SUCCESS;
}

  //! return sharedp tag
MBTag MBParallelComm::sharedp_tag()
{
  if (!sharedpTag) {
    int def_val = -1;
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_PROC_TAG_NAME, 
                                            sizeof(int), 
                                            MB_TAG_DENSE,
                                            MB_TYPE_INTEGER, sharedpTag, 
                                            &def_val, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }
  
  return sharedpTag;
}

  //! return sharedps tag
MBTag MBParallelComm::sharedps_tag()
{
  if (!sharedpsTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_PROCS_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(int), 
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, sharedpsTag, NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }
  
  return sharedpsTag;
}
  
  //! return sharedh tag
MBTag MBParallelComm::sharedh_tag()
{
  if (!sharedhTag) {
    MBEntityHandle def_val = 0;
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_HANDLE_TAG_NAME, 
                                            sizeof(MBEntityHandle), 
                                            MB_TAG_DENSE,
                                            MB_TYPE_HANDLE, sharedhTag, 
                                            &def_val, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return sharedhTag;
}
  
  //! return sharedhs tag
MBTag MBParallelComm::sharedhs_tag()
{  
  if (!sharedhsTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_HANDLES_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(MBEntityHandle), 
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, sharedhsTag, NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }

  return sharedhsTag;
}
  
  //! return pstatus tag
MBTag MBParallelComm::pstatus_tag()
{  
  if (!pstatusTag) {
    unsigned char tmp_pstatus = 0;
    MBErrorCode result = mbImpl->tag_create(PARALLEL_STATUS_TAG_NAME, 
                                            sizeof(unsigned char),
                                            MB_TAG_DENSE,
                                            MB_TYPE_OPAQUE, pstatusTag, 
                                            &tmp_pstatus, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return pstatusTag;
}
  
  //! return partition set tag
MBTag MBParallelComm::partition_tag()
{  
  if (!partitionTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_PARTITION_TAG_NAME, 
                                            sizeof(int),
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, 
                                            partitionTag, 
                                            NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return partitionTag;
}
  
  //! return pcomm tag; passes in impl 'cuz this is a static function
MBTag MBParallelComm::pcomm_tag(MBInterface *impl,
                                bool create_if_missing)
{
  MBTag this_tag = 0;
  MBErrorCode result;
  result = impl->tag_get_handle(PARALLEL_COMM_TAG_NAME, this_tag);
  if ((MB_TAG_NOT_FOUND == result || 0 == this_tag) &&
      create_if_missing) {
    result = impl->tag_create(PARALLEL_COMM_TAG_NAME, 
                              MAX_SHARING_PROCS*sizeof(MBParallelComm*),
                              MB_TAG_SPARSE,
                              MB_TYPE_OPAQUE, this_tag,
                              NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return this_tag;
}

    //! get the indexed pcomm object from the interface
MBParallelComm *MBParallelComm::get_pcomm(MBInterface *impl, const int index) 
{
  MBTag pc_tag = pcomm_tag(impl, false);
  if (0 == pc_tag) return NULL;
  
  MBParallelComm *pc_array[MAX_SHARING_PROCS];
  MBErrorCode result = impl->tag_get_data(pc_tag, 0, 0, (void*)pc_array);
  if (MB_SUCCESS != result) return NULL;
  
  return pc_array[index];
}

MBErrorCode MBParallelComm::get_all_pcomm( MBInterface* impl, std::vector<MBParallelComm*>& list )
{
  MBTag pc_tag = pcomm_tag(impl, false);
  if (0 == pc_tag)
    return MB_TAG_NOT_FOUND;
  
  MBParallelComm *pc_array[MAX_SHARING_PROCS];
  MBErrorCode rval = impl->tag_get_data( pc_tag, 0, 0, pc_array );
  if (MB_SUCCESS != rval)
    return rval;
  
  for (int i = 0; i < MAX_SHARING_PROCS; ++i)
    if (pc_array[i])
      list.push_back( pc_array[i] );
  
  return MB_SUCCESS;
}
  

    //! get the indexed pcomm object from the interface
MBParallelComm *MBParallelComm::get_pcomm( MBInterface *impl, 
                                           MBEntityHandle prtn,
                                           const MPI_Comm* comm ) 
{
  MBErrorCode rval;
  MBParallelComm* result = 0;
  
  MBTag prtn_tag;
  rval = impl->tag_create( PARTITIONING_PCOMM_TAG_NAME, 
                           sizeof(int),
                           MB_TAG_SPARSE,
                           MB_TYPE_INTEGER,
                           prtn_tag,
                           0, true );
  if (MB_SUCCESS != rval)
    return 0;
  
  int pcomm_id;
  rval = impl->tag_get_data( prtn_tag, &prtn, 1, &pcomm_id );
  if (MB_SUCCESS == rval) {
    result= get_pcomm( impl, pcomm_id );
  }
  else if (MB_TAG_NOT_FOUND == rval && comm) {
    result = new MBParallelComm( impl, *comm, &pcomm_id );
    if (!result)
      return 0;
    result->set_partitioning( prtn );
    
    rval = impl->tag_set_data( prtn_tag, &prtn, 1, &pcomm_id );
    if (MB_SUCCESS != rval) {
      delete result;
      result = 0;
    }
  }
  
  return result;
}

MBErrorCode MBParallelComm::set_partitioning( MBEntityHandle set) 
{
  MBErrorCode rval;
  MBTag prtn_tag;
  rval = mbImpl->tag_create( PARTITIONING_PCOMM_TAG_NAME, 
                           sizeof(int),
                           MB_TAG_SPARSE,
                           MB_TYPE_INTEGER,
                           prtn_tag,
                           0, true );
  if (MB_SUCCESS != rval)
    return rval;

    // get my id
  MBParallelComm* pcomm_arr[MAX_SHARING_PROCS];
  MBTag pc_tag = pcomm_tag(mbImpl, false);
  if (0 == pc_tag) 
    return MB_FAILURE;
  MBErrorCode result = mbImpl->tag_get_data(pc_tag, 0, 0, pcomm_arr);
  if (MB_SUCCESS != result) 
    return MB_FAILURE;  
  int id = std::find(pcomm_arr,pcomm_arr+MAX_SHARING_PROCS,this) - pcomm_arr;
  if (id == MAX_SHARING_PROCS)
    return MB_FAILURE;

  MBEntityHandle old = partitioningSet;
  if (old) {
    rval = mbImpl->tag_delete_data( prtn_tag, &old, 1 );
    if (MB_SUCCESS != rval)
      return rval;
    partitioningSet = 0;
  }
  
  if (!set) 
    return MB_SUCCESS;
  
  MBRange contents;
  if (old) {
    rval = mbImpl->get_entities_by_handle( old, contents );
    if (MB_SUCCESS != rval)
      return rval;
  }
  else {
    contents = partition_sets();
  }

  rval = mbImpl->add_entities( set, contents );
  if (MB_SUCCESS != rval)
    return rval;
  
    // store pcomm id on new partition set
  rval = mbImpl->tag_set_data( prtn_tag, &set, 1, &id );
  if (MB_SUCCESS != rval)
    return rval;
  
  partitioningSet = set;
  return MB_SUCCESS;
}
  

  //! return all the entities in parts owned locally
MBErrorCode MBParallelComm::get_part_entities(MBRange &ents, int dim) 
{
  MBErrorCode result;
  
  for (MBRange::iterator rit = partitionSets.begin(); 
       rit != partitionSets.end(); rit++) {
    MBRange tmp_ents;
    if (-1 == dim) 
      result = mbImpl->get_entities_by_handle(*rit, tmp_ents, true);
    else
      result = mbImpl->get_entities_by_dimension(*rit, dim, tmp_ents, true);

    if (MB_SUCCESS != result) return result;
    ents.merge(tmp_ents);
  }
  
  return MB_SUCCESS;
}

  /** \brief Return the rank of the entity owner
   */
MBErrorCode MBParallelComm::get_owner_handle(MBEntityHandle entity,
                                             int &owner,
                                             MBEntityHandle &handle) 
{
    // I'm sure there's a much more efficient logic to this,

    // but I'm tired...
  unsigned char pstat;
  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &entity, 1,
                                            &pstat);
  if (!(pstat & PSTATUS_NOT_OWNED)) {
    owner = proc_config().proc_rank();
    handle = entity;
    return MB_SUCCESS;
  }
  
  int sharing_procs[MAX_SHARING_PROCS];
  MBEntityHandle sharing_handles[MAX_SHARING_PROCS];
  result = mbImpl->tag_get_data(sharedp_tag(), &entity, 1,
                                sharing_procs);
  RRA(" ");
  if (-1 != sharing_procs[0]) {
    owner = sharing_procs[0];
    result = mbImpl->tag_get_data(sharedh_tag(), &entity, 1,
                                  sharing_handles);
    handle = sharing_handles[0];
    return MB_SUCCESS;
  }
  
  result = mbImpl->tag_get_data(sharedps_tag(), &entity, 1,
                                sharing_procs);
  if (MB_SUCCESS == result && -1 != sharing_procs[0]) {
    owner = sharing_procs[0];
    result = mbImpl->tag_get_data(sharedhs_tag(), &entity, 1,
                                  sharing_handles);
    handle = sharing_handles[0];
    return MB_SUCCESS;
  }

  owner = -1;
  handle = 0;
  return MB_FAILURE;
}

MBErrorCode MBParallelComm::get_global_part_count( int& count_out ) const
{
  count_out = globalPartCount;
  return count_out < 0 ? MB_FAILURE : MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_part_owner( int part_id, int& owner ) const
{
  // FIXME: assumes only 1 local part
  owner = part_id;
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_part_id( MBEntityHandle /*part*/, int& id_out ) const
{
  // FIXME: assumes only 1 local part
  id_out = proc_config().proc_rank();
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_part_handle( int id, MBEntityHandle& handle_out ) const
{
  // FIXME: assumes only 1 local part
  if ((unsigned)id != proc_config().proc_rank())
    return MB_ENTITY_NOT_FOUND;
  handle_out = partition_sets().front();
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::create_part( MBEntityHandle& set_out )
{
    // mark as invalid so we know that it needs to be updated
  globalPartCount = -1;
  
    // create set representing part
  MBErrorCode rval = mbImpl->create_meshset( MESHSET_SET, set_out );
  if (MB_SUCCESS != rval)
    return rval;
  
    // set tag on set
    // FIXME: need to assign valid global id
  int val = 0;
  rval = mbImpl->tag_set_data( part_tag(), &set_out, 1, &val );
  if (MB_SUCCESS != rval) {
    mbImpl->delete_entities( &set_out, 1 );
    return rval;
  }
  
  if (get_partitioning()) {
    rval = mbImpl->add_entities( get_partitioning(), &set_out, 1 );
    if (MB_SUCCESS != rval) {
      mbImpl->delete_entities( &set_out, 1 );
      return rval;
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::destroy_part( MBEntityHandle part_id )
{
    // mark as invalid so we know that it needs to be updated
  globalPartCount = -1;
  
  MBErrorCode rval;
  if (get_partitioning()) {
    rval = mbImpl->remove_entities( get_partitioning(), &part_id, 1 );
    if (MB_SUCCESS != rval)
      return rval;
  }
  return mbImpl->delete_entities( &part_id, 1 );
}

MBErrorCode MBParallelComm::collective_sync_partition()
{
  int count = partition_sets().size();
  globalPartCount = 0;
  int err = MPI_Allreduce( &count, &globalPartCount, 1, MPI_INT, MPI_SUM, 
                           proc_config().proc_comm() );
  return err ? MB_FAILURE : MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_part_neighbor_ids( MBEntityHandle part,
                                                   int neighbors_out[MAX_SHARING_PROCS],
                                                   int& num_neighbors_out )
{
  MBErrorCode rval;
  MBRange iface;
  rval = get_interface_sets( part, iface );
  if (MB_SUCCESS != rval)
    return rval;
  
  num_neighbors_out = 0;
  int n, j = 0;
  int tmp[MAX_SHARING_PROCS], curr[MAX_SHARING_PROCS];
  int *parts[2] = { neighbors_out, tmp };
  for (MBRange::iterator i = iface.begin(); i != iface.end(); ++i) {
    rval = get_sharing_parts( *i, curr, n );
    if (MB_SUCCESS != rval)
      return rval;
    std::sort( curr, curr+n );
    assert( num_neighbors_out < MAX_SHARING_PROCS );
    int* k = std::set_union( parts[j], parts[j]+num_neighbors_out,
                             curr, curr + n, parts[1-j] );
    j = 1-j;
    num_neighbors_out = k - parts[j];
  }
  if (parts[j] != neighbors_out)
    std::copy( parts[j], parts[j]+num_neighbors_out, neighbors_out );
    
    
    // remove input part from list
  int id;
  rval = get_part_id( part, id );
  if (MB_SUCCESS == rval) 
    num_neighbors_out = std::remove( neighbors_out, neighbors_out+num_neighbors_out, id ) - neighbors_out;
  return rval;
}

MBErrorCode MBParallelComm::get_interface_sets( MBEntityHandle ,
                                                MBRange& iface_sets_out,
                                                int* adj_part_id )
{
    // FIXME : assumes one part per processor.
    // Need to store part iface sets as children to implement
    // this correctly.
  iface_sets_out = interface_sets();

  if (adj_part_id) {
    int part_ids[MAX_SHARING_PROCS], num_parts;
    MBRange::iterator i = iface_sets_out.begin();
    while (i != iface_sets_out.end()) {
      MBErrorCode rval = get_sharing_parts( *i, part_ids, num_parts );
      if (MB_SUCCESS != rval)
        return rval;
      
      if (std::find(part_ids, part_ids+num_parts, *adj_part_id) - part_ids != num_parts)
        ++i;
      else
        i = iface_sets_out.erase( i );
    }
  }
    
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_owning_part( MBEntityHandle handle,
                                             int& owning_part_id,
                                             MBEntityHandle* remote_handle )
{

  // FIXME : assumes one part per proc, and therefore part_id == rank
  
    // If entity is not shared, then we're the owner.
  unsigned char pstat;
  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &handle, 1,
                                            &pstat);
  if (!(pstat & PSTATUS_NOT_OWNED)) {
    owning_part_id = proc_config().proc_rank();
    if (remote_handle)
      *remote_handle = handle;
    return MB_SUCCESS;
  }
  
    // If entity is shared with one other proc, then
    // sharedp_tag will contain a positive value.
  result = mbImpl->tag_get_data( sharedp_tag(), &handle, 1, &owning_part_id );
  if (MB_SUCCESS != result)
    return result;
  if (owning_part_id != -1) {
      // done?
    if (!remote_handle)
      return MB_SUCCESS;
      
      // get handles on remote processors (and this one)
    return mbImpl->tag_get_data( sharedh_tag(), &handle, 1, remote_handle );
  }
  
    // If here, then the entity is shared with at least two other processors.
    // Get the list from the sharedps_tag
  const void* part_id_list = 0;
  result = mbImpl->tag_get_data( sharedps_tag(), &handle, 1, &part_id_list );
  if (MB_SUCCESS != result)
    return result;
  owning_part_id = ((const int*)part_id_list)[0];
 
    // done?
  if (!remote_handle)
    return MB_SUCCESS;
  
    // get remote handles
  const void* handle_list = 0;
  result = mbImpl->tag_get_data( sharedhs_tag(), &handle, 1, &handle_list );
  if (MB_SUCCESS != result)
    return result;
  
  *remote_handle = ((const MBEntityHandle*)handle_list)[0];
  return MB_SUCCESS;
}    

MBErrorCode MBParallelComm::get_sharing_parts( MBEntityHandle entity,
                                               int part_ids_out[MAX_SHARING_PROCS],
                                               int& num_part_ids_out,
                                               MBEntityHandle remote_handles[MAX_SHARING_PROCS] )
{

  // FIXME : assumes one part per proc, and therefore part_id == rank
  
    // If entity is not shared, then we're the owner.
  unsigned char pstat;
  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &entity, 1,
                                            &pstat);
  if (!(pstat & PSTATUS_SHARED)) {
    part_ids_out[0] = proc_config().proc_rank();
    if (remote_handles)
      remote_handles[0] = entity;
    num_part_ids_out = 1;
    return MB_SUCCESS;
  }
  
    // If entity is shared with one other proc, then
    // sharedp_tag will contain a positive value.
  result = mbImpl->tag_get_data( sharedp_tag(), &entity, 1, part_ids_out );
  if (MB_SUCCESS != result)
    return result;
  if (part_ids_out[0] != -1) {
    
    num_part_ids_out = 2;
    part_ids_out[1] = proc_config().proc_rank();

      // done?
    if (!remote_handles)
      return MB_SUCCESS;
      
      // get handles on remote processors (and this one)
    remote_handles[1] = entity;
    return mbImpl->tag_get_data( sharedh_tag(), &entity, 1, remote_handles );
  }
  
    // If here, then the entity is shared with at least two other processors.
    // Get the list from the sharedps_tag
  result = mbImpl->tag_get_data( sharedps_tag(), &entity, 1, part_ids_out );
  if (MB_SUCCESS != result)
    return result;
    // Count number of valid (positive) entries in sharedps_tag
  for (num_part_ids_out = 0; num_part_ids_out < MAX_SHARING_PROCS &&
       part_ids_out[num_part_ids_out] >= 0; ++num_part_ids_out);
  part_ids_out[num_part_ids_out++] = proc_config().proc_rank();
  
    // done?
  if (!remote_handles)
    return MB_SUCCESS;
  
    // get remote handles
  result = mbImpl->tag_get_data( sharedhs_tag(), &entity, 1, remote_handles );
  remote_handles[num_part_ids_out-1] = entity;
  return result;
}

  /** \brief Find a local ghost entity corresponding to owning proc and handle
   * Given the handle on the owner and the proc, find the local handle to
   * the corresponding entity.  Looks in both the nonowned and owned ghost
   * lists.
   */
MBEntityHandle MBParallelComm::find_ghost_entity(MBEntityHandle owner_handle,
                                                 int owner_proc) 
{
  int ind = get_buffers(owner_proc);
  MBRange sub_range;
  struct GhostStruct &ostruct = sharedEnts[ind];
  if ((!ostruct.remoteHandles.empty() &&
      ostruct.remoteHandles.find(owner_handle) != ostruct.remoteHandles.end()) ||
      (!ostruct.noHandles.empty() &&
       ostruct.noHandles.find(owner_handle) != ostruct.noHandles.end()))
    // ok, it's an existing ghost entity; find local handle
       sub_range = ostruct.localHandles.subset_by_type(
           TYPE_FROM_HANDLE(owner_handle));

  if (sub_range.empty()) return 0;
  
  std::vector<MBEntityHandle> sharedhs(sub_range.size());
  std::vector<int> sharedps(sub_range.size());
  MBErrorCode result = mbImpl->tag_get_data(sharedh_tag(), sub_range,
                                            &sharedhs[0]);
  RRA("Failed to get shared handle tag.");
  result = mbImpl->tag_get_data(sharedp_tag(), sub_range,
                                &sharedps[0]);
  RRA("Failed to get shared proc tag.");

  for (int i = sub_range.size()-1; i >= 0; i--) {
    if (sharedps[i] == owner_proc && sharedhs[i] == owner_handle)
      return sub_range[i];
  }
  
    // not there; search sharedps/sharedhs tags
  sharedhs.resize(MAX_SHARING_PROCS);
  sharedps.resize(MAX_SHARING_PROCS);
  for (MBRange::iterator rit = sub_range.begin();
       rit != sub_range.end(); rit++) {
    result = mbImpl->tag_get_data(sharedhs_tag(), &(*rit), 1, &sharedhs[0]);
    if (MB_TAG_NOT_FOUND == result) continue;
    RRA("Failed to get shared handles tag.");
    result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1, &sharedps[0]);
    RRA("Failed to get shared procs tag.");
    for (int i = 0; i < MAX_SHARING_PROCS; i++) {
      if (sharedhs[i] == owner_handle && sharedps[i] == owner_proc) return *rit;
      else if (-1 == sharedps[i]) break;
    }
  }
  
    // not found
  return 0;
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
  result = pcomm.pack_buffer(all_mesh, false, true, true, false, whole_range, buff_size);
  PM;


    // allocate space in the buffer
  pcomm.buffer_size(buff_size);

    // pack the actual buffer
  int actual_buff_size;
  result = pcomm.pack_buffer(whole_range, false, true, false, false, all_mesh, 
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
  result = pcomm2.unpack_buffer(all_mesh, store_remote_handles, from_proc);
  PM;
  
  std::cout << "ENTITIES UNPACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);
  
  std::cout << "Success, processor " << mbImpl->proc_rank() << "." << std::endl;
  
  return 1;
}
#endif
