#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

#include "MBMeshSet.hpp"
#include "AEntityFactory.hpp"


/*****************************************************************************************
 *                          Helper Function Declarations                                 *
 *****************************************************************************************/

/**\brief Insert into parent/child list */
static inline 
MBMeshSet::Count insert_in_vector( const MBMeshSet::Count count, 
                                   MBMeshSet::CompactList& list,
                                   const MBEntityHandle h,
                                   int &result );

/**\brief Remvoe from parent/child list */
static inline
MBMeshSet::Count remove_from_vector( const MBMeshSet::Count count, 
                                     MBMeshSet::CompactList& list,
                                     const MBEntityHandle h,
                                     int &result );


/**\brief Resize MBMeshSet::CompactList.  Returns pointer to storage */
static MBEntityHandle* resize_compact_list( MBMeshSet::Count& count,
                                            MBMeshSet::CompactList& clist,
                                            size_t new_list_size );
/**\brief Methods to insert/remove range-based data from contents list.
 *        Templatized to operate on both MBRange and set-based MBMeshSets.
 */
template <typename pair_iter_t> class range_tool
{
public:
  /** Insert range-based data into range-based MBMeshSet */
  inline static MBErrorCode ranged_insert_entities( MBMeshSet::Count& count, 
                                                    MBMeshSet::CompactList& clist, 
                                                    pair_iter_t begin, 
                                                    pair_iter_t end, 
                                                    MBEntityHandle my_handle, 
                                                    AEntityFactory* adj );
  
  /** Remove range-based data from range-based MBMeshSet */
  inline static MBErrorCode ranged_remove_entities( MBMeshSet::Count& count, 
                                                    MBMeshSet::CompactList& clist, 
                                                    pair_iter_t begin, 
                                                    pair_iter_t end, 
                                                    MBEntityHandle my_handle, 
                                                    AEntityFactory* adj );

  /** Insert range-based data into list-based MBMeshSet */
  inline static MBErrorCode vector_insert_entities( MBMeshSet::Count& count, 
                                                    MBMeshSet::CompactList& clist, 
                                                    pair_iter_t begin, 
                                                    pair_iter_t end, 
                                                    MBEntityHandle my_handle, 
                                                    AEntityFactory* adj );
};

/** Remove MBRange of handles fromr vector-based MBMeshSet */
static MBErrorCode vector_remove_range( MBMeshSet::Count& count, 
                                        MBMeshSet::CompactList& clist, 
                                        const MBRange& range, 
                                        MBEntityHandle my_handle, 
                                        AEntityFactory* adj );

/** Remove range-based MBMeshSet contents from vector-based MBMeshSet */
static MBErrorCode vector_remove_ranges( MBMeshSet::Count& count, 
                                         MBMeshSet::CompactList& clist, 
                                         const MBEntityHandle* pair_list,
                                         size_t num_pairs,
                                         MBEntityHandle my_handle, 
                                         AEntityFactory* adj );

/** Remove unsorted array of handles from vector-based MBMeshSet */
static MBErrorCode vector_remove_vector( MBMeshSet::Count& count, 
                                         MBMeshSet::CompactList& clist, 
                                         const MBEntityHandle* vect,
                                         size_t vect_size,
                                         MBEntityHandle my_handle, 
                                         AEntityFactory* adj );

/** Insert unsorted array of handles into vector-based MBMeshSet */
static MBErrorCode vector_insert_vector( MBMeshSet::Count& count, 
                                         MBMeshSet::CompactList& clist, 
                                         const MBEntityHandle* vect,
                                         size_t vect_size,
                                         MBEntityHandle my_handle, 
                                         AEntityFactory* adj );

/** Convert unsorted array of handles into array of ranged [begin,end] pairs */
static void convert_to_ranges( const MBEntityHandle* vect_in, size_t vect_in_len,
                               std::vector<MBEntityHandle>& vect_out );


/*****************************************************************************************
 *                             Parent/Child Operations                                   *
 *****************************************************************************************/

static inline 
MBMeshSet::Count insert_in_vector( const MBMeshSet::Count count, 
                                MBMeshSet::CompactList& list,
                                const MBEntityHandle h,
                                int &result )
{
  switch (count) {
    case MBMeshSet::ZERO:
      list.hnd[0] = h;
      result = true;
      return MBMeshSet::ONE;
    case MBMeshSet::ONE:
      if (list.hnd[0] == h) {
        result = false;
        return MBMeshSet::ONE;
      }
      else {
        result = true;
        list.hnd[1] = h;
        return MBMeshSet::TWO;
      }
    case MBMeshSet::TWO:
      if (list.hnd[0] == h || list.hnd[1] == h) {
        result = false;
        return MBMeshSet::TWO;
      }
      else {
        MBEntityHandle* ptr = (MBEntityHandle*)malloc(3*sizeof(MBEntityHandle));
        ptr[0] = list.hnd[0];
        ptr[1] = list.hnd[1];
        ptr[2] = h;
        list.ptr[0] = ptr;
        list.ptr[1] = ptr + 3;
        result = true;
        return MBMeshSet::MANY;
      }
    case MBMeshSet::MANY:
      if (std::find( list.ptr[0], list.ptr[1], h ) != list.ptr[1]) {
        result = false;
      }
      else {
        int size = list.ptr[1] - list.ptr[0];
        list.ptr[0] = (MBEntityHandle*)realloc( list.ptr[0], (size+1)*sizeof(MBEntityHandle) );
        list.ptr[0][size] = h;
        list.ptr[1] = list.ptr[0] + size + 1;
        result = true;
      }
      return MBMeshSet::MANY;
  }

  return MBMeshSet::ZERO;
}

static inline
MBMeshSet::Count remove_from_vector( const MBMeshSet::Count count, 
                                  MBMeshSet::CompactList& list,
                                  const MBEntityHandle h,
                                  int &result )
{
  switch (count) {
    case MBMeshSet::ZERO:
      result = false;
      return MBMeshSet::ZERO;
    case MBMeshSet::ONE:
      if (h == list.hnd[0]) {
        result = true;
        return MBMeshSet::ZERO;
      }
      else {
        result = false;
        return MBMeshSet::ONE;
      }
    case MBMeshSet::TWO:
      if (h == list.hnd[0]) {
        list.hnd[0] = list.hnd[1];
        result = true;
        return MBMeshSet::ONE;
      } 
      else if (h == list.hnd[1]) {
        result = true;
        return MBMeshSet::ONE;
      }
      else {
        result = false;
        return MBMeshSet::TWO;
      }
    case MBMeshSet::MANY: {
      MBEntityHandle *i, *j, *p;
      i = std::find( list.ptr[0], list.ptr[1], h );
      if (i == list.ptr[1]) {
        result = false;
        return MBMeshSet::MANY;
      }
      
      result = true;
      p = list.ptr[1] - 1;
      while (i != p) {
        j = i + 1;
        *i = *j;
        i = j;
      }
      int size = p - list.ptr[0];
      if (size == 2) {
        p = list.ptr[0];
        list.hnd[0] = p[0];
        list.hnd[1] = p[1];
        free( p );
        return MBMeshSet::TWO;
      }
      else {
        list.ptr[0] = (MBEntityHandle*)realloc( list.ptr[0], size*sizeof(MBEntityHandle) );
        list.ptr[1] = list.ptr[0] + size;
        return MBMeshSet::MANY;
      }
    }
  }

  return MBMeshSet::ZERO;
}


int MBMeshSet::add_parent( MBEntityHandle parent )
{ 
  int result = 0;
  mParentCount = insert_in_vector( (Count)mParentCount, parentMeshSets, parent, result );
  return result;
}
int MBMeshSet::add_child( MBEntityHandle child )
{ 
  int result = 0;
  mChildCount = insert_in_vector( (Count)mChildCount, childMeshSets, child, result );
  return result;
}

int MBMeshSet::remove_parent( MBEntityHandle parent )
{ 
  int result = 0;
  mParentCount = remove_from_vector( (Count)mParentCount, parentMeshSets, parent, result );
  return result;
}
int MBMeshSet::remove_child( MBEntityHandle child )
{ 
  int result = 0;
  mChildCount = remove_from_vector( (Count)mChildCount, childMeshSets, child, result );
  return result;
}


/*****************************************************************************************
 *                          Flag Conversion Operations                                   *
 *****************************************************************************************/

MBErrorCode MBMeshSet::convert( unsigned flags, MBEntityHandle my_handle, AEntityFactory* adj )
{
  MBErrorCode rval = MB_SUCCESS;
  if ((mFlags & MESHSET_TRACK_OWNER) && !(flags & MESHSET_TRACK_OWNER))
    rval = remove_adjacencies( my_handle, adj );
  else if (!(mFlags & MESHSET_TRACK_OWNER) && (flags & MESHSET_TRACK_OWNER))
    rval = create_adjacencies( my_handle, adj );
  if (MB_SUCCESS != rval)
    return rval;

  if ((mFlags & MESHSET_ORDERED) && !(flags & MESHSET_ORDERED)) {
    size_t datalen;
    MBEntityHandle* data = get_contents(datalen);
    if (datalen) {
      std::vector<MBEntityHandle> list( datalen );
      memcpy( &list[0], data, datalen*sizeof(MBEntityHandle) );
      int num_ents = num_entities();
      Count count = (Count)mContentCount;
      data = resize_compact_list( count, contentList, num_ents );
      mContentCount = count;
      assert( list.size() % 2 == 0 );
      std::vector<MBEntityHandle>::iterator i = list.begin();
      while (i != list.end()) {
        MBEntityHandle h = *i; ++i;
        MBEntityHandle e = *i; ++i;
        for (; h <= e; ++h) {
          *data = h; 
          ++data;
        }
      }
    }
  }
  else if (!(mFlags & MESHSET_ORDERED) && (flags & MESHSET_ORDERED)) {
    size_t datalen;
    MBEntityHandle* data = get_contents(datalen);
    if (datalen) {
      std::vector<MBEntityHandle> ranges;
      convert_to_ranges( data, datalen, ranges );
      Count count = (Count)mContentCount;
      data = resize_compact_list( count, contentList, ranges.size() );
      mContentCount = count;
      memcpy( data, &ranges[0], ranges.size()*sizeof(MBEntityHandle) );
    }
  }

  return MB_SUCCESS;
}

MBErrorCode MBMeshSet::create_adjacencies( MBEntityHandle my_handle, AEntityFactory* adj )
{
  MBErrorCode rval = MB_SUCCESS;;
  size_t count;
  const MBEntityHandle *const ptr = get_contents( count );
  const MBEntityHandle *const end = ptr + count;
  if (vector_based()) {
    for (const MBEntityHandle* i = ptr; i != end; ++i) {
      rval = adj->add_adjacency( *i, my_handle, false );
      if (MB_SUCCESS != rval) {
        for (const MBEntityHandle* j = ptr; j != i; ++j) 
          adj->remove_adjacency( *j, my_handle );
        return rval;
      }
    }
  }
  else {
    assert( 0 == count % 2 );
    for (const MBEntityHandle* i = ptr; i != end; i += 2) {
      for (MBEntityHandle h = i[0]; h <= i[1]; ++h) {
        rval = adj->add_adjacency( h, my_handle, false );
        if (MB_SUCCESS != rval) {
          for (MBEntityHandle j = i[0]; j < h; ++j)
            adj->remove_adjacency( j, my_handle );
          for (const MBEntityHandle* j = ptr; j != i; j += 2)
            for (MBEntityHandle k = j[0]; k <= j[1]; ++k)
              adj->remove_adjacency( k, my_handle );
          return rval;
        }
      }
    }
  }
  return MB_SUCCESS;
}

MBErrorCode MBMeshSet::remove_adjacencies( MBEntityHandle my_handle, AEntityFactory* adj )
{
  size_t count;
  const MBEntityHandle *const ptr = get_contents( count );
  const MBEntityHandle *const end = ptr + count;
  if (vector_based()) {
    for (const MBEntityHandle* i = ptr; i != end; ++i)
      adj->remove_adjacency( *i, my_handle );
  }
  else {
    assert( 0 == count % 2 );
    for (const MBEntityHandle* i = ptr; i != end; i += 2)
      for (MBEntityHandle h = i[0]; h <= i[1]; ++h)
        adj->remove_adjacency( h, my_handle );
  }
  return MB_SUCCESS;
}


/*****************************************************************************************
 *                          Contents Modifiction Methods                                 *
 *****************************************************************************************/

static MBEntityHandle* resize_compact_list( MBMeshSet::Count& count,
                                            MBMeshSet::CompactList& clist,
                                            size_t new_list_size )
{
  if (count <= 2) {
    if (new_list_size <= 2) {
      count = (MBMeshSet::Count)new_list_size;
      return clist.hnd;
    }
    else {
      MBEntityHandle* list = (MBEntityHandle*)malloc( new_list_size*sizeof(MBEntityHandle) );
      list[0] = clist.hnd[0];
      list[1] = clist.hnd[1];
      clist.ptr[0] = list;
      clist.ptr[1] = list + new_list_size;
      count = MBMeshSet::MANY;
      return list;
    }
  }
  else if (new_list_size > 2) {
    if (new_list_size > (size_t)(clist.ptr[1] - clist.ptr[0]))
      clist.ptr[0] = (MBEntityHandle*)realloc( clist.ptr[0], new_list_size*sizeof(MBEntityHandle) );
    clist.ptr[1] = clist.ptr[0] + new_list_size;
    count = MBMeshSet::MANY;
    return clist.ptr[0];
  }
  else {
    MBEntityHandle* list = clist.ptr[0];
    clist.hnd[0] = list[0];
    clist.hnd[1] = list[1];
    free(list);
    count = (MBMeshSet::Count)new_list_size;
    return clist.hnd;
  }
}

template <typename pair_iter_t> inline MBErrorCode
range_tool<pair_iter_t>::ranged_insert_entities( MBMeshSet::Count& count, 
                                                 MBMeshSet::CompactList& clist, 
                                                 pair_iter_t begin, 
                                                 pair_iter_t end, 
                                                 MBEntityHandle my_handle, 
                                                 AEntityFactory* adj )
{
    //first pass:
    // 1) merge existing ranges 
    // 2) count number of new ranges that must be inserted
  ptrdiff_t insert_count = 0;
  MBEntityHandle *list;
  size_t list_size;
  if (count < MBMeshSet::MANY) {
    list = clist.hnd;
    list_size = count;
  }
  else {
    list = clist.ptr[0];
    list_size = clist.ptr[1] - clist.ptr[0];
  }

  MBEntityHandle* list_write = list;
  MBEntityHandle *const list_end = list + list_size, *list_read = list;
  pair_iter_t i = begin;
  
  while(i != end) {
    
      // if there are holes in the current array, shuffle blocks 
      // down until we find the next block to merge with or insert before
    if (list_read != list_write) {
      while (list_read != list_end && i->second + 1 < list_read[0]) {
      	list_write[0] = list_read[0];
        list_write[1] = list_read[1];
        list_write += 2;
        list_read += 2;
      }
    }
      // otherwise do a binary search
    else {
      list_write = std::lower_bound( list_write, list_end, i->first - 1 );
      	// if in middle of range block (odd index), back up to start of block
      list_write -= (list_write - list)%2;
      list_read = list_write;
    }
    
      // handle any straight insertions of range blocks
    for ( ; i != end && (list_read == list_end || i->second+1 < list_read[0]); ++i) {
        // If we haven't removed any range pairs, we don't have space to
        // insert here.  Defer the insertion until later.
      if (list_read == list_write) {
        ++insert_count;
      }
      else {
        if (adj) 
          for (MBEntityHandle j = i->first; j <= i->second; ++j)
            adj->add_adjacency( j, my_handle, false );

        list_write[0] = i->first;
        list_write[1] = i->second;
        list_write += 2;
      }
    }
    if (i == end)
      break;
    
      // check if we need to prepend to the current range block
    if (i->first >= list_read[0]) 
      list_write[0] = list_read[0];
    else {
      if (adj)
        for (MBEntityHandle h = i->first; h < list_read[0]; ++h)
          adj->add_adjacency( h, my_handle, false );
      list_write[0] = i->first;
    }
    list_write[1] = list_read[1];
    list_read += 2;
    
      // discard any input blocks already in the set
    for (; i != end && i->second <= list_write[1]; ++i);
    if (i == end) {
      list_write += 2;
      break;
    }
    
      // merge subsequent blocks in meshset
    for (;list_read != list_end && list_read[0]+1 <= i->second; list_read += 2) {
      if (adj)
      	for (MBEntityHandle h = list_write[1]+1; h < list_read[0]; ++h)
      	  adj->add_adjacency( h, my_handle, false );
      list_write[1] = list_read[1];
    }
    
      // check if we need to append to current meshset block
    if (i->second > list_write[1]) {
      if (adj)
      	for (MBEntityHandle h = list_write[1]+1; h <= i->second; ++h)
      	  adj->add_adjacency( h, my_handle, false );
      list_write[1] = i->second;
    }
    
    ++i;
    list_write += 2;
  }

    // shuffle down entries to fill holes
  if (list_read == list_write) 
    list_read = list_write = list_end;
  else while(list_read < list_end) {
    list_write[0] = list_read[0];
    list_write[1] = list_read[1];
    list_read += 2;
    list_write += 2;
  }

    // adjust allocated array size
  const size_t occupied_size = list_write - list;
  const size_t new_list_size = occupied_size + 2*insert_count;
  list = resize_compact_list( count, clist, new_list_size );
    // done?
  if (!insert_count)
    return MB_SUCCESS;

    // Second pass: insert non-mergable range pairs
    // All range pairs in the input are either completely disjoint from
    // the ones in the mesh set and must be inserted or are entirely contained
    // within range pair in the mesh set.
  assert( begin != end ); // can't have items to insert if given empty input list
  pair_iter_t ri = end; --ri;
  list_write = list + new_list_size - 2;
  list_read = list + occupied_size - 2;
  for ( ; list_write >= list; list_write -= 2 ) {
    if (list_read >= list) {
      while (ri->first >= list_read[0] && ri->second <= list_read[1]) {
        assert(ri != begin);
        --ri;
      }
    
      if (list_read[0] > ri->second) {
        list_write[0] = list_read[0];
        list_write[1] = list_read[1];
        list_read -= 2;
        continue;
      }
    }
    
    assert( insert_count > 0 );
    if (adj) 
      for (MBEntityHandle j = ri->first; j <= ri->second; ++j) 
        adj->add_adjacency( j, my_handle, false );
    list_write[0] = ri->first;
    list_write[1] = ri->second;

      // don't have reverse iterator, so check before decrement
      // if insert_count isn't zero, must be more in range
    if (0 == --insert_count) {
      assert( list_read == list_write-2 );
      break;
    }
    else {
      --ri;
    }
  }

  assert(!insert_count);
  return MB_SUCCESS;
}
  
template <typename pair_iter_t> inline MBErrorCode
range_tool<pair_iter_t>::ranged_remove_entities( MBMeshSet::Count& count, 
                                                 MBMeshSet::CompactList& clist, 
                                                 pair_iter_t begin, 
                                                 pair_iter_t end, 
                                                 MBEntityHandle my_handle, 
                                                 AEntityFactory* adj )
{
    //first pass:
    // 1) remove (from) existing ranges 
    // 2) count number of ranges that must be split
  ptrdiff_t split_count = 0;
  MBEntityHandle *list;
  size_t list_size;
  if (count < MBMeshSet::MANY) {
    list = clist.hnd;
    list_size = count;
  }
  else {
    list = clist.ptr[0];
    list_size = clist.ptr[1] - clist.ptr[0];
  }

  MBEntityHandle* list_write = list;
  MBEntityHandle *const list_end = list + list_size, *list_read = list;
  pair_iter_t i = begin;
  
  while(list_read != list_end && i != end) {
    
    while (i != end && i->second < list_read[0])
      ++i;
    if (i == end)
      break;
    
      // if there are holes in the current array, shuffle blocks 
      // down until we find the next block to remove
    if (list_read != list_write) {
      while (list_read != list_end && i->second < list_read[0]) {
      	list_write[0] = list_read[0];
        list_write[1] = list_read[1];
        list_write += 2;
        list_read += 2;
      }
    }
      // otherwise do a binary search
    else {
      list_write = std::lower_bound( list_write, list_end, i->first );
      	// if in middle of range block (odd index), back up to start of block
      list_write -= (list_write - list)%2;
      list_read = list_write;
    }
    
      // if everything remaning is past end of set contents...
    if (list_read == list_end) 
      break;
      
      // skip any remove pairs that aren't in the list
    if (i->second < list_read[0]) {
      ++i;
      continue;
    }
    
      // Begin by assuming that we will keep the entire block
    list_write[0] = list_read[0];
    list_write[1] = list_read[1];
    list_read += 2;
    
    for (; i != end && i->first <= list_write[1]; ++i) {
      if (i->first <= list_write[0]) {
          // remove whole block
        if (i->second >= list_write[1]) {
          if (adj)
            for (MBEntityHandle h = list_write[0]; h <= list_write[1]; ++h)
              adj->remove_adjacency( h, my_handle );
          list_write -= 2;
          break;
        }
          // remove from start of block
        else if (i->second >= list_write[0]) {
          if (adj)
            for (MBEntityHandle h = list_write[0]; h <= i->second; ++h)
              adj->remove_adjacency( h, my_handle );
          list_write[0] = i->second + 1;
        }
      }
      else if (i->first <= list_write[1]) {
          // remove from end of block
        if (i->second >= list_write[1]) {
          if (adj)
            for (MBEntityHandle h = i->first; h <= list_write[1]; ++h)
              adj->remove_adjacency( h, my_handle );
          list_write[1] = i->first - 1;
          //list_write += 2;
          break;
        }
          // split block
        else {
          if (adj)
            for (MBEntityHandle h = i->first; h <= i->second; ++h)
              adj->remove_adjacency( h, my_handle );

          if (list_read - list_write <= 2) {
            ++split_count;
            continue;
          }
          else {
            list_write[3] = list_write[1];
            list_write[1] = i->first - 1;
            list_write[2] = i->second + 1;
            list_write += 2;
          }
        }
      }
    }
    list_write += 2;
  }

    // shuffle down entries to fill holes
  if (list_read == list_write) 
    list_read = list_write = list_end;
  else 
    while(list_read < list_end) {
      list_write[0] = list_read[0];
      list_write[1] = list_read[1];
      list_read += 2;
      list_write += 2;
    }

    // adjust allocated array size
  const size_t occupied_size = list_write - list;
  const size_t new_list_size = occupied_size + 2*split_count;
  list = resize_compact_list( count, clist, new_list_size );
    // done?
  if (!split_count)
    return MB_SUCCESS;

    // Second pass: split range pairs
    // All range pairs in the input are either already removed or
    // require one of the existing range pairs to be split
  assert( begin != end ); // can't have ranges to split if given empty input list
  pair_iter_t ri = end; --ri;
  list_write = list + new_list_size - 2;
  list_read = list + occupied_size - 2;
  for ( ; list_write >= list; list_write -= 2 ) {
    if (list_read >= list) {
      while (ri->second > list_read[1]) {
        assert(ri != begin);
        --ri;
      }
    
      if (list_read[0] > ri->second) {
        list_write[0] = list_read[0];
        list_write[1] = list_read[1];
        list_read -= 2;
        continue;
      }
    }
    
    assert( split_count > 0 );
    list_write[0] = ri->second + 1;
    list_write[1] = list_read[1];
    list_read[1] = ri->first - 1;

      // don't have reverse iterator, so check before decrement
      // if insert_count isn't zero, must be more in range
    if (0 == --split_count) {
      assert( list_read == list_write-2 );
      break;
    }
    else {
      --ri;
    }
  }

  assert(!split_count);
  return MB_SUCCESS;
}


template <typename pair_iter_t> inline MBErrorCode
range_tool<pair_iter_t>::vector_insert_entities( MBMeshSet::Count& count, 
                                                 MBMeshSet::CompactList& clist, 
                                                 pair_iter_t begin, 
                                                 pair_iter_t end, 
                                                 MBEntityHandle my_handle, 
                                                 AEntityFactory* adj )
{
  const size_t init_size = count < MBMeshSet::MANY ? count : clist.ptr[1] - clist.ptr[0];
  size_t add_size = 0;
  for (pair_iter_t i = begin; i != end; ++i)
    add_size += i->second - i->first + 1;
  MBEntityHandle* list = resize_compact_list( count, clist, init_size + add_size );
  MBEntityHandle* li = list + init_size;

  for (pair_iter_t i = begin; i != end; ++i) {
    for (MBEntityHandle h = i->first; h <= i->second; ++h) {
      if (adj)
        adj->add_adjacency( h, my_handle, false );
      *li = h;
      ++li;
    }
  }

  return MB_SUCCESS;
}

static MBErrorCode vector_remove_range( MBMeshSet::Count& count, 
                                        MBMeshSet::CompactList& clist, 
                                        const MBRange& range, 
                                        MBEntityHandle my_handle, 
                                        AEntityFactory* adj )
{
  MBEntityHandle *list;
  size_t list_size;
  if (count < MBMeshSet::MANY) {
    list = clist.hnd;
    list_size = count;
  }
  else {
    list = clist.ptr[0];
    list_size = clist.ptr[1] - clist.ptr[0];
  }

  const MBEntityHandle * const list_end = list + list_size;
  MBEntityHandle* list_write = list;
  for (const MBEntityHandle* list_read = list; list_read != list_end; ++list_read) {
    if (range.find(*list_read) == range.end()) { // keep
      *list_write = *list_read;
      ++list_write;
    }
    else if (adj) {    
      adj->remove_adjacency( *list_read, my_handle );
    }
  }

  resize_compact_list( count, clist, list_write - list );
  return MB_SUCCESS;
}

static MBErrorCode vector_remove_ranges( MBMeshSet::Count& count, 
                                         MBMeshSet::CompactList& clist, 
                                         const MBEntityHandle* pair_list,
                                         size_t num_pairs,
                                         MBEntityHandle my_handle, 
                                         AEntityFactory* adj )
{
  MBEntityHandle *list;
  size_t list_size;
  if (count < MBMeshSet::MANY) {
    list = clist.hnd;
    list_size = count;
  }
  else {
    list = clist.ptr[0];
    list_size = clist.ptr[1] - clist.ptr[0];
  }

  const MBEntityHandle *const list_end = list + list_size, 
                       *const input_end = pair_list + 2*num_pairs;
  MBEntityHandle* list_write = list;
  for (const MBEntityHandle* list_read = list; list_read != list_end; ++list_read) {
    const MBEntityHandle* ptr = std::lower_bound( pair_list, input_end, *list_read );
    if ((ptr != input_end && (*ptr == *list_read || (ptr - pair_list)%2)) && // if in delete list
        std::find(list_read+1, list_end, *list_read) == list_end) { // and is last occurance in list 
        // only remove adj if no previous occurance
      if (adj && std::find(list, list_write, *list_read) == list_write)
        adj->remove_adjacency( *list_read, my_handle );
    }
    else {
      *list_write = *list_read;
      ++list_write;
    }
  }

  resize_compact_list( count, clist, list_write - list );
  return MB_SUCCESS;
}

static MBErrorCode vector_remove_vector( MBMeshSet::Count& count, 
                                         MBMeshSet::CompactList& clist, 
                                         const MBEntityHandle* vect,
                                         size_t vect_size,
                                         MBEntityHandle my_handle, 
                                         AEntityFactory* adj )
{
  MBEntityHandle *list;
  size_t list_size;
  if (count < MBMeshSet::MANY) {
    list = clist.hnd;
    list_size = count;
  }
  else {
    list = clist.ptr[0];
    list_size = clist.ptr[1] - clist.ptr[0];
  }

  const MBEntityHandle *const list_end = list + list_size, 
                       *const input_end = vect + vect_size;
  MBEntityHandle* list_write = list;
  for (const MBEntityHandle* list_read = list; list_read != list_end; ++list_read) {
    if (std::find(vect, input_end, *list_read) != input_end && // if in delete list
        std::find(list_read+1, list_end, *list_read) == list_end) { // and is last occurance in list 
        // only remove adj if no previous occurance?
      if (adj ) // && std::find(list, list_write, *list_read) == list_write)
        adj->remove_adjacency( *list_read, my_handle );
    }
    else {
      *list_write = *list_read;
      ++list_write;
    }
  }

  resize_compact_list( count, clist, list_write - list );
  return MB_SUCCESS;
}

static MBErrorCode vector_insert_vector( MBMeshSet::Count& count, 
                                         MBMeshSet::CompactList& clist, 
                                         const MBEntityHandle* vect,
                                         size_t vect_size,
                                         MBEntityHandle my_handle, 
                                         AEntityFactory* adj )
{
  const size_t orig_size = count < MBMeshSet::MANY ? count : clist.ptr[1] - clist.ptr[0];
  MBEntityHandle* list = resize_compact_list( count, clist, orig_size + vect_size );
  if (adj) 
    for (size_t i = 0; i < vect_size; ++i)
      adj->add_adjacency( vect[i], my_handle, false );
  memcpy( list+orig_size, vect, sizeof(MBEntityHandle)*vect_size );
  return MB_SUCCESS;
}

MBErrorCode MBMeshSet::insert_entity_ranges( const MBEntityHandle* range_vect, size_t len, MBEntityHandle my_h, AEntityFactory* adj )
{
  typedef const std::pair<MBEntityHandle,MBEntityHandle>* pair_vect_t;
  pair_vect_t pair_vect = reinterpret_cast<pair_vect_t>(range_vect);
  MBMeshSet::Count count = static_cast<MBMeshSet::Count>(mContentCount);
  MBErrorCode rval;
  if (!vector_based())
    rval = range_tool<pair_vect_t>::ranged_insert_entities( count, contentList,  pair_vect, 
                                             pair_vect + len/2, my_h, tracking() ? adj : 0 );
  else
    rval = range_tool<pair_vect_t>::vector_insert_entities( count, contentList,  pair_vect, 
                                             pair_vect + len/2, my_h, tracking() ? adj : 0 );
  mContentCount = count;
  return rval;
}

MBErrorCode MBMeshSet::insert_entity_ranges( const MBRange& range, MBEntityHandle my_h, AEntityFactory* adj )
{
  MBErrorCode rval;
  MBMeshSet::Count count = static_cast<MBMeshSet::Count>(mContentCount);
  if (!vector_based())
    rval = range_tool<MBRange::const_pair_iterator>::ranged_insert_entities( count, 
                             contentList, range.const_pair_begin(), range.const_pair_end(), 
                             my_h, tracking() ? adj : 0 );
  else
    rval = range_tool<MBRange::const_pair_iterator>::vector_insert_entities( count, 
                             contentList, range.const_pair_begin(), range.const_pair_end(), 
                             my_h, tracking() ? adj : 0 );
  mContentCount = count;
  return rval;
}

MBErrorCode MBMeshSet::remove_entity_ranges( const MBEntityHandle* range_vect, size_t len, MBEntityHandle my_h, AEntityFactory* adj )
{
  MBErrorCode rval;
  MBMeshSet::Count count = static_cast<MBMeshSet::Count>(mContentCount);
  if (vector_based()) 
    rval = vector_remove_ranges( count, contentList, range_vect, len/2, my_h, 
                                 tracking() ? adj : 0 );
  else {
    typedef const std::pair<MBEntityHandle,MBEntityHandle>* pair_vect_t;
    pair_vect_t pair_vect = reinterpret_cast<pair_vect_t>(range_vect);
    rval = range_tool<pair_vect_t>::ranged_remove_entities( count, contentList, pair_vect, 
                                           pair_vect + len/2, my_h, tracking() ? adj : 0 );
  }
  mContentCount = count;
  return rval;
}

MBErrorCode MBMeshSet::remove_entity_ranges( const MBRange& range, MBEntityHandle my_h, AEntityFactory* adj )
{
  MBErrorCode rval;
  MBMeshSet::Count count = static_cast<MBMeshSet::Count>(mContentCount);
  if (vector_based()) 
    rval = vector_remove_range( count, contentList, range, my_h, tracking() ? adj : 0 );
  else 
    rval = range_tool<MBRange::const_pair_iterator>::ranged_remove_entities( count, 
                         contentList, range.const_pair_begin(), range.const_pair_end(), 
                         my_h, tracking() ? adj : 0 );
  mContentCount = count;
  return rval;
}


MBErrorCode MBMeshSet::intersect( const MBMeshSet* other, MBEntityHandle my_handle, AEntityFactory* adj )
{
  MBErrorCode rval;
  if (!vector_based() && !other->vector_based()) {
    size_t other_count = 0;
    const MBEntityHandle* other_vect = other->get_contents( other_count );
    if (!other_count)
      return clear( my_handle, adj );
    assert(0 == other_count%2);
    
    std::vector<MBEntityHandle> compliment;
    compliment.reserve( other_count + 4 );
    if (*other_vect > 0) {
      compliment.push_back( 0 );
      compliment.push_back( *other_vect - 1 );
    }
    ++other_vect;
    const MBEntityHandle *const other_end = other_vect + other_count - 2;
    for (; other_vect < other_end; other_vect += 2) {
      compliment.push_back( other_vect[0] + 1 );
      compliment.push_back( other_vect[1] - 1 );
    }
    if (*other_vect < ~(MBEntityHandle)0) {
      compliment.push_back( *other_vect + 1 );
      compliment.push_back( ~(MBEntityHandle)0 );
    }
    
    return remove_entity_ranges( &compliment[0], compliment.size(), my_handle, adj );
  }
  else {
    MBRange my_ents, other_ents;
    rval = get_entities(my_ents);
    if (MB_SUCCESS != rval)
      return rval;
    rval = other->get_entities(other_ents);
    return remove_entities( my_ents.subtract(other_ents), my_handle, adj );
  }
}

static void convert_to_ranges( const MBEntityHandle* vect_in, size_t vect_in_len,
                               std::vector<MBEntityHandle>& vect_out )
{
  vect_out.reserve( 2*vect_in_len );
  vect_out.resize( vect_in_len );
  std::copy( vect_in, vect_in+vect_in_len, vect_out.begin() );
  std::sort( vect_out.begin(), vect_out.end() );

    // duplicate all entries
  vect_out.resize( vect_out.size() * 2 );
  for (long i = vect_out.size() - 1; i >= 0; --i) 
    vect_out[i] = vect_out[i/2];
   
    // compact adjacent ranges
  std::vector<MBEntityHandle>::iterator r = vect_out.begin(), w = vect_out.begin();
  while (r != vect_out.end()) {
    *w = *r;
    ++w; 
    ++r;
    *w = *r;
    ++r;
    
    while (r != vect_out.end() && *w + 1 == *r) {
      ++r;
      *w = *r;
      ++r;
    }
    ++w;
  }
  
    // remove extra space
  vect_out.erase( w, vect_out.end() );
}

MBErrorCode MBMeshSet::insert_entity_vector( const MBEntityHandle* vect, size_t len, MBEntityHandle my_h, AEntityFactory* adj )
{
  MBMeshSet::Count count = static_cast<MBMeshSet::Count>(mContentCount);
  MBErrorCode rval;
  if (vector_based())
    rval = vector_insert_vector( count, contentList, vect, len, my_h, tracking() ? adj : 0 );
  else {
    std::vector<MBEntityHandle> rangevect;
    convert_to_ranges( vect, len, rangevect );
    typedef const std::pair<MBEntityHandle,MBEntityHandle>* pair_vect_t;
    pair_vect_t pair_vect = reinterpret_cast<pair_vect_t>(&rangevect[0]);
    rval = range_tool<pair_vect_t>::ranged_insert_entities( count, contentList, pair_vect, 
                                 pair_vect + rangevect.size()/2, my_h, tracking() ? adj : 0 );
  }
  mContentCount = count;
  return rval;
}

MBErrorCode MBMeshSet::remove_entity_vector( const MBEntityHandle* vect, size_t len, MBEntityHandle my_h, AEntityFactory* adj )
{
  MBMeshSet::Count count = static_cast<MBMeshSet::Count>(mContentCount);
  MBErrorCode rval;
  if (vector_based())
    rval = vector_remove_vector( count, contentList, vect, len, my_h, tracking() ? adj : 0 );
  else {
    std::vector<MBEntityHandle> rangevect;
    convert_to_ranges( vect, len, rangevect );
    typedef const std::pair<MBEntityHandle,MBEntityHandle>* pair_vect_t;
    pair_vect_t pair_vect = reinterpret_cast<pair_vect_t>(&rangevect[0]);
    rval = range_tool<pair_vect_t>::ranged_remove_entities( count, contentList, pair_vect, 
                                pair_vect + rangevect.size()/2, my_h, tracking() ? adj : 0 );
  }
  mContentCount = count;
  return rval;
}



bool MBMeshSet::replace_entities( MBEntityHandle my_handle,
                                  const MBEntityHandle* entities,
                                  size_t num_entities,
                                  AEntityFactory* adjfact )
{
  assert(0 == num_entities%2);
  if (vector_based()) {
    bool was_contained = false;
    size_t count;
    MBEntityHandle* vect = get_contents( count );
    MBEntityHandle* const vect_end = vect+count;
    for (size_t i = 0; i < num_entities; i+=2) {
      for (MBEntityHandle* p = vect; p != vect_end; ++p ) {
        if (*p == entities[i]) {
          if (tracking()) {
            adjfact->remove_adjacency( *p, my_handle );
            adjfact->add_adjacency( entities[i+1], my_handle, false );
          }
          *p = entities[i+1];
          was_contained = true;
        }
      }
    }
    return was_contained;
  }
  else {
    std::vector<MBEntityHandle> swap_list;
    swap_list.reserve( num_entities / 2 );

      // get list of handles to remove
    size_t count;
    MBEntityHandle* vect = get_contents( count );
    MBEntityHandle* const vect_end = vect+count;
    for (size_t i = 0; i < num_entities; i+=2) {
      MBEntityHandle* p = std::lower_bound(vect, vect_end, entities[i]);
      if (p != vect_end && (*p == entities[i] || (p-vect)%2 == 1))
        swap_list.push_back(entities[i]);
    }
    if (swap_list.empty())
      return false;
    
      // remove entities
    remove_entities( &swap_list[0], swap_list.size(), my_handle, adjfact );
    
      // get list of handles to add
    std::vector<MBEntityHandle>::iterator si = swap_list.begin();
    for (size_t i = 0; i < num_entities; ++i) {
      if (entities[i] == *si) {
        *si = entities[i+1];
        if (++si == swap_list.end())
          break;
      }
    }
    
      // add entities
    add_entities( &swap_list[0], swap_list.size(), my_handle, adjfact );
    return true;
  }
}


/*****************************************************************************************
 *                                  Misc. Methods                                        *
 *****************************************************************************************/

unsigned long MBMeshSet::get_memory_use() const
{
  unsigned long result = 0;
  if (mParentCount == MANY)
    result += parentMeshSets.ptr[1] - parentMeshSets.ptr[0];
  if (mChildCount == MANY)
    result += childMeshSets.ptr[1] - childMeshSets.ptr[0];
  if (mContentCount == MANY)
    result += contentList.ptr[1] - contentList.ptr[0];
  return sizeof(MBEntityHandle)*result;
}
