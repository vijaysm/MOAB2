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

class TagServer;
class SequenceManager;

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
  

    //! communicate entities from/to this range
  MBErrorCode communicate_entities(const int from_proc, const int to_proc,
                                   MBRange &entities,
                                   const bool adjacencies = false,
                                   const bool tags = true);
  
  MBErrorCode broadcast_entities( const int from_proc,
                                  MBRange& entities,
                                  const bool adjacencies = false,
                                  const bool tags = true );

    /** Resolve shared entities between processors
     * Resolve shared entities between processors for entities in proc_ents,
     * by comparing global id tag values on vertices on skin of elements in
     * proc_ents.  Shared entities are assigned a tag that's either
     * PARALLEL_SHARED_PROC_TAG_NAME, which is 2 integers in length, or 
     * PARALLEL_SHARED_PROCS_TAG_NAME, whose length depends on the maximum
     * number of sharing processors.  Values in these tags denote the ranks
     * of sharing processors, and the list ends with the value -10*#procs.
     *
     * If shared_dim is input as -1 or not input, a value one less than the
     * maximum dimension of entities in proc_ents is used.
     *
     * \param proc_ents Entities for which to resolve shared entities
     * \param shared_dim Maximum dimension of shared entities to look for
     */
  MBErrorCode resolve_shared_ents(MBRange &proc_ents, 
                                  int shared_dim = 0);
  
    /** Resolve shared entities between processors
     * Same as resolve_shared_ents(MBRange&), except works for
     * all entities in instance of dimension dim.  
     *
     * If dim = -1 or no dim input, uses entities of maximal
     * dimension (3) in the instance.
     *
     * If shared_dim is input as -1 or not input, a value one less than the
     * maximum dimension of entities is used.

     * \param dim Maximal dimension of entities to be resolved
     * \param shared_dim Maximum dimension of shared entities to look for
     */
  MBErrorCode resolve_shared_ents(int dim = -1, 
                                  int shared_dim = 0);
  
    /** Get entities shared with other processors, based on 
     * PARALLEL_SHARED_PROC_TAG_NAME and PARALLEL_SHARED_PROCS_TAG_NAME.
     *
     * \param dim Dimension of entities shared with other processors
     * \param shared_ents Shared entities returned in this range
     */
  MBErrorCode get_shared_entities(int dim,
                                  MBRange &shared_ents);

    //! pack a buffer (stored in this class instance) with ALL data for these entities
  MBErrorCode pack_buffer(MBRange &entities, 
                          const bool adjacencies,
                          const bool tags,
                          const bool just_count,
                          MBRange &whole_range,
                          int &buff_size);
  
    //! unpack a buffer; assume information is already in myBuffer
  MBErrorCode unpack_buffer(MBRange &entities);

    //! set the buffer size; return true if size actually changed
  bool buffer_size(const unsigned int new_size);

    //! take the buffer from this instance; switches with vector passed in
  void take_buffer(std::vector<unsigned char> &new_buff);

    //! Get proc config for this communication object
  const MBProcConfig &proc_config() const {return procConfig;}
  
      
private:

  int num_subranges(const MBRange &this_range);
  
  MBErrorCode pack_entities(MBRange &entities,
                            MBRange::const_iterator &start_rit,
                            MBRange &whole_range,
                            unsigned char *&buff_ptr,
                            int &count,
                            const bool just_count);
  
  MBErrorCode unpack_entities(unsigned char *&buff_ptr,
                              MBRange &entities);
  
  MBErrorCode pack_sets(MBRange &entities,
                        MBRange::const_iterator &start_rit,
                        MBRange &whole_range,
                        unsigned char *&buff_ptr,
                        int &count,
                        const bool just_count);
  
  MBErrorCode unpack_sets(unsigned char *&buff_ptr,
                          MBRange &entities);
  
  MBErrorCode pack_adjacencies(MBRange &entities,
                               MBRange::const_iterator &start_rit,
                               MBRange &whole_range,
                               unsigned char *&buff_ptr,
                               int &count,
                               const bool just_count);

  MBErrorCode unpack_adjacencies(unsigned char *&buff_ptr,
                                 MBRange &entities);
  
  MBErrorCode pack_tags(MBRange &entities,
                        MBRange::const_iterator &start_rit,
                        MBRange &whole_range,
                        unsigned char *&buff_ptr,
                        int &count,
                        const bool just_count);

  MBErrorCode unpack_tags(unsigned char *&buff_ptr,
                          MBRange &entities);
  

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

    //! types of ranges to be communicated
  std::vector<MBEntityType> entTypes;

    //! ranges to be communicated
  std::vector<MBRange> allRanges;
  
    //! vertices per entity in ranges
  std::vector<int> vertsPerEntity;

    //! sets to be communicated
  MBRange setRange;
  
    //! ranges from sets to be communicated
  std::vector<MBRange> setRanges;
  
    //! sizes of vector-based sets to be communicated
  std::vector<int> setSizes;

    //! tags to be communicated
  std::vector<MBTag> allTags;

    //! ranges from sparse tags to be communicated
  std::vector<MBRange> tagRanges;

    //! vector of set options for transferred sets
  std::vector<unsigned int> optionsVec;
  
    //! numbers of parents/children for transferred sets
  std::vector<int> setPcs;
};

#endif
