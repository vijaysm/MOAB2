/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/**
 * \class MBParallelComm
 * \brief Parallel communications in MOAB
 * \author Tim Tautges
 *
 *  This class implements methods to communicate mesh between processors
 *
 */

#ifndef MB_PARALLEL_COMM_HPP
#define MB_PARALLEL_COMM_HPP

#include "MBForward.hpp"
#include "MBRange.hpp"
#include "MBProcConfig.hpp"
#include <map>
#include "math.h"

extern "C" {
  struct tuple_list;
}

class TagServer;
class SequenceManager;
template <typename KeyType, typename ValType, ValType NullVal> class RangeMap;
typedef RangeMap<MBEntityHandle, MBEntityHandle, 0> HandleMap;

#define MAX_SHARING_PROCS 10  

class MBParallelComm 
{
public:

    //! constructor
  MBParallelComm(MBInterface *impl,
                 MPI_Comm comm = MPI_COMM_WORLD);

    //! constructor taking packed buffer, for testing
  MBParallelComm(MBInterface *impl,
                 std::vector<unsigned char> &tmp_buff,
                 MPI_Comm comm = MPI_COMM_WORLD);

  static unsigned char PROC_SHARED, PROC_OWNER;
  
    //! assign a global id space, for largest-dimension or all entities (and
    //! in either case for vertices too)
  MBErrorCode assign_global_ids(MBEntityHandle this_set,
                                const int dimension,
                                const int start_id = 1,
                                const bool largest_dim_only = true,
                                const bool parallel = true);

    //! check for global ids; based only on tag handle being there or not;
    //! if it's not there, create them for the specified dimensions
  MBErrorCode check_global_ids(MBEntityHandle this_set,
                               const int dimension, 
                               const int start_id = 1,
                               const bool largest_dim_only = true,
                               const bool parallel = true);
  

    /** \brief send entities to another processor, optionally waiting until it's done
     *
     * Send entities to another processor, with adjs, sets, and tags.  
     * If store_remote_handles is true, this call receives back handles assigned to
     * entities sent to destination processor and stores them in sharedh_tag or 
     * sharedhs_tag.
     * \param to_proc Destination processor
     * \param orig_ents Entities requested to send
     * \param adjs If true, send adjacencies for equiv entities (currently unsupported)
     * \param tags If true, send tag values for all tags assigned to entities
     * \param store_remote_handles If true, also recv message with handles on destination processor
     * \param final_ents Range containing all entities sent
     * \param wait_all If true, wait until all messages received/sent complete
     */
  MBErrorCode send_entities(const int to_proc,
                            MBRange &orig_ents,
                            const bool adjs,
                            const bool tags,
                            const bool store_remote_handles,
                            MBRange &final_ents,
                            bool wait_all = true);
  
    /** \brief Receive entities from another processor, optionally waiting until it's done
     *
     * Receive entities from another processor, with adjs, sets, and tags.  
     * If store_remote_handles is true, this call sends back handles assigned to
     * the entities received.
     * \param from_proc Source processor
     * \param store_remote_handles If true, send message with new entity handles to source processor
     * \param final_ents Range containing all entities received
     * \param wait_all If true, wait until all messages received/sent complete
     */
  MBErrorCode recv_entities(const int from_proc,
                            const bool store_remote_handles,
                            MBRange &final_ents,
                            bool wait_all = true);
  
    /** \brief Exchange ghost cells with neighboring procs
     * Neighboring processors are those sharing an interface 
     * with this processor.  All entities of dimension ghost_dim
     * within num_layers of interface, measured going through bridge_dim,
     * are exchanged.  See MeshTopoUtil::get_bridge_adjacencies for description
     * of bridge adjacencies.  If wait_all is false and store_remote_handles
     * is true, MPI_Request objects are available in the sendReqs[2*MAX_SHARING_PROCS] 
     * member array, with inactive requests marked as MPI_REQUEST_NULL.  If
     * store_remote_handles or wait_all is false, this function returns after 
     * all entities have been received and processed.
     * \param ghost_dim Dimension of ghost entities to be exchanged
     * \param bridge_dim Dimension of entities used to measure layers from interface
     * \param num_layers Number of layers of ghosts requested
     * \param store_remote_handles If true, send message with new entity handles to source processor
     * \param wait_all If true, function does not return until all send buffers
     *       are cleared.
     */
  MBErrorCode exchange_ghost_cells(int ghost_dim, int bridge_dim, 
                                   int num_layers,
                                   bool store_remote_handles,
                                   bool wait_all = true);

    /** \brief Broadcast all entities resident on from_proc to other processors
     * This function assumes remote handles are *not* being stored, since (usually)
     * every processor will know about the whole mesh.
     * \param from_proc Processor having the mesh to be broadcast
     * \param entities On return, the entities sent or received in this call
     * \param adjacencies If true, adjacencies are sent for equiv entities (currently unsupported)
     * \param tags If true, all non-default-valued tags are sent for sent entities
     */
  MBErrorCode broadcast_entities(const int from_proc,
                                 MBRange& entities,
                                 const bool adjacencies = false,
                                 const bool tags = true );

    /** \brief Resolve shared entities between processors
     *
     * Resolve shared entities between processors for entities in proc_ents,
     * by comparing global id tag values on vertices on skin of elements in
     * proc_ents.  Shared entities are assigned a tag that's either
     * PARALLEL_SHARED_PROC_TAG_NAME, which is 1 integer in length, or 
     * PARALLEL_SHARED_PROCS_TAG_NAME, whose length depends on the maximum
     * number of sharing processors.  Values in these tags denote the ranks
     * of sharing processors, and the list ends with the value -1.
     *
     * If shared_dim is input as -1 or not input, a value one less than the
     * maximum dimension of entities in proc_ents is used.
     *
     * \param proc_ents Entities for which to resolve shared entities
     * \param shared_dim Maximum dimension of shared entities to look for
     */
  MBErrorCode resolve_shared_ents(MBRange &proc_ents, 
                                  int resolve_dim = -1,
                                  int shared_dim = -1);
  
    /** \brief Resolve shared entities between processors
     *
     * Same as resolve_shared_ents(MBRange&), except works for
     * all entities in instance of dimension dim.  
     *
     * If shared_dim is input as -1 or not input, a value one less than the
     * maximum dimension of entities is used.

     * \param dim Dimension of entities in the partition
     * \param shared_dim Maximum dimension of shared entities to look for
     */
  MBErrorCode resolve_shared_ents(int resolve_dim = 3, 
                                  int shared_dim = -1);
  
    /** \brief Get entities with the given pstatus bit(s) set
     * Returns any entities whose pstatus tag value v satisfies (v & pstatus_val)
     *
     * \param dim Dimension of entities to be returned, or -1 if any
     * \param pstatus_val pstatus value of desired entities
     * \param pstatus_ents Entities returned from function
     */
  MBErrorCode get_pstatus_entities(int dim,
                                   unsigned char pstatus_val,
                                   MBRange &pstatus_ents);
  
    /** \brief Get entities on an inter-processor interface and of specified dimension
     * If other_proc is -1, any interface entities are returned.  If dim is -1,
     * entities of all dimensions on interface are returned.
     * \param other_proc Rank of processor for which interface entities are requested
     * \param dim Dimension of interface entities requested
     * \param iface_ents Entities returned from function
     */
  MBErrorCode get_iface_entities(int other_proc,
                                 int dim,
                                 MBRange &iface_ents);
  
    //! pack the buffer with ALL data for orig_ents; return entities actually
    //! packed (including reference sub-entities) in final_ents
  MBErrorCode pack_buffer(MBRange &orig_ents,
                          const bool adjacencies,
                          const bool tags,
                          const bool store_remote_handles,
                          const int to_proc,
                          MBRange &final_ents,
                          std::vector<unsigned char> &buff,
                          int &buff_size);
  
    //! unpack a buffer; assume information is already in myBuffer
  MBErrorCode unpack_buffer(unsigned char *buff_ptr,
                            const bool store_remote_handles, 
                            int from_proc,
                            MBRange &entities);

    //! Get proc config for this communication object
  const MBProcConfig &proc_config() const {return procConfig;}
  
    //! return the tags used to indicate shared procs and handles
  MBErrorCode get_shared_proc_tags(MBTag &sharedp_tag,
                                   MBTag &sharedps_tag,
                                   MBTag &sharedh_tag,
                                   MBTag &sharedhs_tag,
                                   MBTag &pstatus_tag);
  
    //! return sharedp tag
  MBTag sharedp_tag();
  
    //! return sharedps tag
  MBTag sharedps_tag();
  
    //! return sharedh tag
  MBTag sharedh_tag();
  
    //! return sharedhs tag
  MBTag sharedhs_tag();
  
    //! return pstatus tag
  MBTag pstatus_tag();
  
    //! return iface_set tag
  MBTag iface_sets_tag();
  
    //! return partitions set tag
  MBTag partition_tag();
  
private:

  int num_subranges(const MBRange &this_range);

    //! get (and possibly allocate) buffers for messages to/from to_proc; returns
    //! index of to_proc in buffProcs vector
  int get_buffers(int to_proc);

    //! pack entities (with adjs, tags too) and send to to_proc; also post a recv
    //! if store_remote_handles is true; pass back MPI_Request objects for these
    //! messages, along with entities actually sent; this function calls:
    //! - MPI_Irecv if store_remote_handles is true (posts recv for remote handles)
    //! - MPI_Send if send buffer is larger than INITIAL_BUFF_SIZE
    //! - MPI_Isend for sending entities
  MBErrorCode pack_send_entities(const int to_proc,
                                 MBRange &orig_ents,
                                 const bool adjacencies,
                                 const bool tags,
                                 const bool store_remote_handles,
                                 const bool iface_layer,
                                 std::vector<unsigned char> &send_buff,
                                 std::vector<unsigned char> &recv_buff,
                                 MPI_Request &send_req,
                                 MPI_Request &recv_req,
                                 MBRange &final_ents);
  
    //! use integer size in buffer to resize buffer, then post an
    //! Irecv to get message
  MBErrorCode recv_size_buff(const int from_proc,
                             std::vector<unsigned char> &recv_buff,
                             MPI_Request &recv_req);
  
    //! process contents of receive buffer to get new entities; if store_remote_handles
    //! is true, also Isend (using send_buff) handles for these entities back to 
    //! source proc, returning request handle in &send_req; if iface_layer is true,
    //! don't instantiate the entities, just check to see if they correspond to 
    //! existing entities, and if not, set corresponding recd_ents handle to zero
  MBErrorCode recv_unpack_entities(const int from_proc,
                                   const bool store_remote_handles,
                                   const bool iface_layer,
                                   std::vector<unsigned char> &recv_buff,
                                   std::vector<unsigned char> &send_buff,
                                   MPI_Request &send_req,
                                   MBRange &recd_ents);
  
    //! for all the entities in the received buffer; for each, save
    //! entities in this instance which match connectivity, or zero if none found
  MBErrorCode unpack_iface_entities(unsigned char *&buff_ptr, 
                                    const int from_proc,
                                    std::vector<MBEntityHandle> &recd_ents);
  
  MBErrorCode pack_entities(MBRange &entities,
                            MBRange::const_iterator &start_rit,
                            MBRange &whole_range,
                            unsigned char *&buff_ptr,
                            int &count,
                            const bool just_count,
                            const bool store_remote_handles,
                            const int to_proc,
                            std::vector<MBEntityType> &ent_types, 
                            std::vector<MBRange> &all_ranges, 
                            std::vector<int> &verts_per_entity);
  
  MBErrorCode unpack_entities(unsigned char *&buff_ptr,
                              const bool store_remote_handles,
                              const int from_proc,
                              MBRange &entities);
  
  MBErrorCode pack_sets(MBRange &entities,
                        MBRange::const_iterator &start_rit,
                        MBRange &whole_range,
                        unsigned char *&buff_ptr,
                        int &count,
                        const bool just_count,
                        const bool store_handles,
                        const int to_proc,
                        MBRange &set_range,
                        std::vector<MBRange> &set_ranges,
                        std::vector<int> &set_sizes,
                        std::vector<unsigned int> &options_vec);
  
  MBErrorCode unpack_sets(unsigned char *&buff_ptr,
                          MBRange &entities,
                          const bool store_handles,
                          const int to_proc);
  
  MBErrorCode pack_adjacencies(MBRange &entities,
                               MBRange::const_iterator &start_rit,
                               MBRange &whole_range,
                               unsigned char *&buff_ptr,
                               int &count,
                               const bool just_count,
                               const bool store_handles,
                               const int to_proc);

  MBErrorCode unpack_adjacencies(unsigned char *&buff_ptr,
                                 MBRange &entities,
                                 const bool store_handles,
                                 const int from_proc);
  
  MBErrorCode pack_tags(MBRange &entities,
                        MBRange::const_iterator &start_rit,
                        MBRange &whole_range,
                        unsigned char *&buff_ptr,
                        int &count,
                        const bool just_count,
                        const bool store_handles,
                        const int to_proc,
                        std::vector<MBTag> &all_tags,
                        std::vector<MBRange> &tag_ranges);

  MBErrorCode unpack_tags(unsigned char *&buff_ptr,
                          MBRange &entities,
                          const bool store_handles,
                          const int to_proc);
  
  MBErrorCode tag_shared_verts(tuple_list &shared_verts,
                               MBRange *skin_ents,
                               std::map<std::vector<int>, MBRange> &proc_nranges,
                               MBRange &proc_verts);
  
  MBErrorCode tag_shared_ents(int resolve_dim,
                              int shared_dim,
                              tuple_list &shared_verts,
                              MBRange *skin_ents,
                              std::map<std::vector<int>, MBRange> &proc_nranges);

    // each iterate in proc_nranges contains a set of procs and the entities *possibly*
    // on the interface between those procs; this function makes sets for each,
    // and tags the set with the procs sharing it; interface sets are optionally
    // returned; NOTE: a subsequent step is used to verify entities on the interface
    // and remove them if they're not shared
  MBErrorCode create_interface_sets(std::map<std::vector<int>, MBRange> &proc_nranges,
                                    int resolve_dim, int shared_dim,
                                    MBRange *iface_sets_ptr = NULL);

    // after verifying shared entities, now parent/child links between sets can be established
  MBErrorCode create_iface_pc_links(MBRange &iface_sets);
  
    //! resolve remote handles for shared non-vertex ents, assuming
    //! this has already been done for vertices
  MBErrorCode resolve_ent_remote_handles(MBRange &iface_sets);
  
    //! pack a range map with keys in this_range and values a contiguous series
    //! of handles starting at actual_start
  MBErrorCode pack_range_map(MBRange &this_range, MBEntityHandle actual_start,
                             HandleMap &handle_map);

    //! remove from the range all ents not owned by this proc or already
    //! shared with to_proc
  MBErrorCode remove_nonowned_shared(MBRange &ents,
                                     int to_proc,
                                     bool owned_test = true);
  
    //! for a given interface set, gets a number of layers of bridge entities
    //! of dimension to_dim going through bridge dimension bridge_dim
  MBErrorCode get_ghost_layers(MBEntityHandle iface_set,
                               int to_dim, int bridge_dim,
                               int num_layers,
                               MBRange &to_ents);
  
    //! returns true if the set is an interface shared with to_proc
  bool is_iface_proc(MBEntityHandle this_set, int to_proc);
  
    //! for any remote_handles set to zero, remove corresponding sent_ents from
    //! iface_sets corresponding to from_proc
  MBErrorCode update_iface_sets(MBRange &iface_sets, 
                                MBRange &sent_ents,
                                std::vector<MBEntityHandle> &remote_handles, 
                                int from_proc);
  
    //! replace handles in from_vec with corresponding handles on
    //! to_proc (by checking shared[p/h]_tag and shared[p/h]s_tag;
    //! if no remote handle and new_ents is non-null, substitute
    //! instead CREATE_HANDLE(MBMAXTYPE, index) where index is handle's
    //! position in new_ents
  MBErrorCode get_remote_handles(const bool store_remote_handles,
                                 MBEntityHandle *from_vec, 
                                 MBEntityHandle *to_vec_tmp,
                                 int num_ents, int to_proc,
                                 const MBRange &new_ents);
  
    //! same as other version, except from_range and to_range should be
    //! different here
  MBErrorCode get_remote_handles(const bool store_remote_handles,
                                 const MBRange &from_range, 
                                 MBRange &to_range,
                                 int to_proc,
                                 const MBRange &new_ents);
  
    //! same as other version, except packs range into vector
  MBErrorCode get_remote_handles(const bool store_remote_handles,
                                 const MBRange &from_range, 
                                 MBEntityHandle *to_vec,
                                 int to_proc,
                                 const MBRange &new_ents);
  
    //! goes through from_vec, and for any with type MBMAXTYPE, replaces with
    //! new_ents value at index corresponding to id of entity in from_vec
  MBErrorCode get_local_handles(MBEntityHandle *from_vec, 
                                int num_ents,
                                const MBRange &new_ents);

    //! same as above except puts results in range
  MBErrorCode get_local_handles(const MBRange &remote_handles,
                                MBRange &local_handles,
                                const MBRange &new_ents);
  
    //! adjust shared proc tags/handles to incude from_proc and remote_range
  MBErrorCode set_remote_data(MBRange &local_range,
                              MBRange &remote_range,
                              int from_proc);
  
    //! adjust shared proc tags/handles to incude from_proc and remote_range
  MBErrorCode set_remote_data(MBEntityHandle *local_ents,
                              MBEntityHandle *remote_ents,
                              int num_ents,
                              int other_proc);
  
    //! remove a remote processor and the entity's handle
  MBErrorCode rmv_remote_proc(MBEntityHandle ent,
                              int *remote_procs,
                              MBEntityHandle *remote_hs,
                              int remote_proc);
  
    //! add a remote processor and the entity's handle
  MBErrorCode add_remote_proc(MBEntityHandle ent,
                              int *remote_procs,
                              MBEntityHandle *remote_hs,
                              int remote_proc,
                              MBEntityHandle remote_handle);
  
    //! returns the interface sets and union of sharing procs
  MBErrorCode get_iface_sets_procs(MBRange &iface_sets,
                                   std::vector<int> &sharing_procs);
  
    //! MB interface associated with this writer
  MBInterface *mbImpl;

    //! Proc config object, keeps info on parallel stuff
  MBProcConfig procConfig;
  
    //! Tag server, so we can get more info about tags
  TagServer *tagServer;
  
    //! Sequence manager, to get more efficient access to entities
  SequenceManager *sequenceManager;
  
    //! data buffer used to communicate
  std::vector<unsigned char> myBuffer;

    //! more data buffers, proc-specific
  std::vector<unsigned char> ownerRBuffs[MAX_SHARING_PROCS],
    ownerSBuffs[MAX_SHARING_PROCS], ghostRBuffs[MAX_SHARING_PROCS],
    ghostSBuffs[MAX_SHARING_PROCS];

    //! request objects, may be used if store_remote_handles is used
  MPI_Request sendReqs[2*MAX_SHARING_PROCS];

  std::vector<int> buffProcs;
  
    //! tags used to save sharing procs and handles
  MBTag sharedpTag, sharedpsTag, sharedhTag, sharedhsTag, pstatusTag, 
    ifaceSetsTag, partitionTag;

    //! interface sets, one set per unique combination of procs
  MBRange ifaceSets;
  
    //! ghost sets (sets of ghost entities), one set per unique combination of procs
  MBRange ghostSets;
  
    //! ghosted sets (sets of ghosted entities), one set per unique combination of procs
  MBRange ghostedSets;
  
};

inline MBErrorCode MBParallelComm::get_shared_proc_tags(MBTag &sharedp,
                                                        MBTag &sharedps,
                                                        MBTag &sharedh,
                                                        MBTag &sharedhs,
                                                        MBTag &pstatus) 
{
  sharedp = sharedp_tag();
  sharedps = sharedps_tag();
  sharedh = sharedh_tag();
  sharedhs = sharedhs_tag();
  pstatus = pstatus_tag();
  
  return MB_SUCCESS;
}

#endif
