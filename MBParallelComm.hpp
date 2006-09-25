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

class TagServer;
class EntitySequenceManager;

class MBParallelComm 
{
public:

    //! constructor
  MBParallelComm(MBInterface *impl, TagServer *tag_server, 
                 EntitySequenceManager *sequence_manager);

    //! constructor taking packed buffer, for testing
  MBParallelComm(MBInterface *impl, TagServer *tag_server, 
                 EntitySequenceManager *sequence_manager,
                 std::vector<unsigned char> &tmp_buff);

    //! communicate entities from/to this range
  MBErrorCode communicate_entities(const int from_proc, const int to_proc,
                                   MBRange &entities,
                                   const bool adjacencies = false,
                                   const bool tags = true);
  
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
  
    //! Tag server, so we can get more info about tags
  TagServer *tagServer;
  
    //! Sequence manager, to get more efficient access to entities
  EntitySequenceManager *sequenceManager;
  
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
