#ifndef TAG_COMPARE_HPP
#define TAG_COMPARE_HPP

#include "TagInfo.hpp"
#include "VarLenTag.hpp"

/* OPAQUE FUNCTORS */

/** Test fixed-length opaque tags for equality */
class TagBytesEqual {
  private:
    const void* value;
    int size;
  public:
    TagBytesEqual( const void* v, int s ) : value(v), size(s) {}
    bool operator()( const void* data ) const
      { return !memcmp(value, data, size); }
};
/** Test if fixed-length opaque tag values are less than a value */
class TagBytesLess {
  private:
    const void* value;
    int size;
  public:
    TagBytesLess( const void* v, int s ) : value(v), size(s) {}
    bool operator()( const void* data ) const
      { return 0 < memcmp(value, data, size); }
};
/** Test variable-length opaque tags for equality */
class TagVarBytesEqual {
  private:
    const void* value;
    int size;
  public:
    TagVarBytesEqual( const void* v, int s ) : value(v), size(s) {}
    bool operator()( const void* data ) const {
      const VarLenTag* vdata = reinterpret_cast<const VarLenTag*>(data);
      return (int)vdata->size() == size && !memcmp(value, vdata->data(), size); 
    }
};
/** Test if variable-length opaque tag values are less than a value */
class TagVarBytesLess {
  private:
    const void* value;
    int size;
  public:
    TagVarBytesLess( const void* v, int s ) : value(v), size(s) {}
    bool operator()( const void* data ) const {
      const VarLenTag* vdata = reinterpret_cast<const VarLenTag*>(data);
      if ((int)vdata->size() < size) 
        return 0 <= memcmp( vdata->data(), value, vdata->size() );
      else
        return 0 < memcmp( vdata->data(), value, size );
    }
};


/* TEMPLATE FUNCTORS */


/** Compare fixed-length tags containing a known data type */
template <typename T>
class TagTypeEqual {
  private:
    const T* value;
    int size;
  public:
    TagTypeEqual( const void* v, int s ) 
      : value(reinterpret_cast<const T*>(v)), 
        size(s/sizeof(T)) 
        {}
        
    bool operator()( const void* data ) const { 
      const T* ddata = reinterpret_cast<const T*>(data);
      for (int i = 0; i < size; ++i)
        if (value[i] != ddata[i])
          return false;
      return true;
    }
};

/** Compare fixed-length tags containing a known data type */
template <typename T>
class TagTypeLess {
  private:
    const T* value;
    int size;
  public:
    TagTypeLess( const void* v, int s ) 
      : value(reinterpret_cast<const T*>(v)), 
        size(s/sizeof(T)) 
        {}
    
    bool operator()( const void* data ) const {
      const T* ddata = reinterpret_cast<const T*>(data);
      for (int i = 0; i < size; ++i)
        if (value[i] <= ddata[i])
          return false;
      return true;
    }
};

/** Compare single-value tags containing a known data type
 * Optimization of TagTypeEqual for 1-value case. 
 */
template <typename T>
class TagOneTypeEqual {
  private:
    T value;
    int size;
  public:
    TagOneTypeEqual( const void* v ) 
      : value(*reinterpret_cast<const T*>(v))
        {}
        
    bool operator()( const void* data ) const { 
      const T* ddata = reinterpret_cast<const T*>(data);
      return *ddata == value;
    }
};

/** Compare single-value tags containing a known data type
 * Optimization of TagTypeLess for 1-value case. 
 */
template <typename T>
class TagOneTypeLess {
  private:
    T value;
    int size;
  public:
    TagOneTypeLess( const void* v ) 
      : value(*reinterpret_cast<const T*>(v))
        {}
    
    bool operator()( const void* data ) const {
      const T* ddata = reinterpret_cast<const T*>(data);
      return *ddata < value;
    }
};

/** Compare variable-length tags containing a known data type */
template <typename T>
class TagVarTypeEqual
{
  private:
    const T* value;
    int size;
  public:
    TagVarTypeEqual( const void* v, int s ) 
      : value(reinterpret_cast<const T*>(v)), 
        size(s/sizeof(T)) 
        {}
        
    bool operator()( const void* data ) const {
      const VarLenTag* vdata = reinterpret_cast<const VarLenTag*>(data);
      if (vdata->size() != size * sizeof(T))
        return false;
      const T* ddata = reinterpret_cast<const T*>(vdata->data());
      for (int i = 0; i < size; ++i)
        if (value[i] != ddata[i])
          return false;
      return true;
    }
};

/** Compare variable-length tags containing a known data type */
template <typename T>
class TagVarTypeLess
{
  private:
    const T* value;
    int size;
  public:
    TagVarTypeLess( const void* v, int s ) 
      : value(reinterpret_cast<const T*>(v)), 
        size(s/sizeof(T)) 
        {}
    bool operator()( const void* data ) const {
      const VarLenTag* vdata = reinterpret_cast<const VarLenTag*>(data);
      const T* ddata = reinterpret_cast<const T*>(vdata->data());
      if ((int)vdata->size() < sizeof(T)*size) {
        for (int i = 0; i < vdata->size()/sizeof(T); ++i)
          if (value[i] < ddata[i])
            return false;
      }
      else {
        for (int i = 0; i < vdata->size()/sizeof(T); ++i)
          if (value[i] <= ddata[i])
            return false;
      }
      return true;
    }
};

/* TYPE FUNCTORS */

typedef TagBytesEqual        TagIntsEqual;
typedef TagVarBytesEqual     TagVarIntsEqual;
typedef TagTypeLess    <int> TagIntsLess;
typedef TagVarTypeLess <int> TagVarIntsLess;
typedef TagOneTypeEqual<int> TagOneIntEqual;
typedef TagOneTypeLess <int> TagOneIntLess;

typedef TagBytesEqual                   TagHandlesEqual;
typedef TagVarBytesEqual                TagVarHandlesEqual;
typedef TagTypeLess    <MBEntityHandle> TagHandlesLess;
typedef TagVarTypeLess <MBEntityHandle> TagVarHandlesLess;
typedef TagOneTypeEqual<MBEntityHandle> TagOneHandleEqual;
typedef TagOneTypeLess <MBEntityHandle> TagOneHandleLess;

typedef TagTypeEqual   <double> TagDoublesEqual;
typedef TagVarTypeEqual<double> TagVarDoublesEqual;
typedef TagTypeLess    <double> TagDoublesLess;
typedef TagVarTypeLess <double> TagVarDoublesLess;
typedef TagOneTypeEqual<double> TagOneDoubleEqual;
typedef TagOneTypeLess <double> TagOneDoubleLess;

/* SEARCHING */

template <class Functor,
          class IteratorType>
void find_tag_values( Functor compare,
                      IteratorType begin,
                      IteratorType end,
                      MBRange& results )
{
  MBRange::iterator insert = results.begin();
  for (IteratorType i = begin; i != end; ++i) 
    if (compare( i->second ))
      insert = results.insert( insert, i->first, i->first );
}

template <class Functor,
          class IteratorType>
void find_tag_values( Functor compare,
                      IteratorType begin,
                      IteratorType end,
                      std::vector<MBEntityHandle>& results )
{
  MBRange::iterator insert = results.begin();
  for (IteratorType i = begin; i != end; ++i) 
    if (compare( i->second ))
      results.push_back( i->first );
}

/** Find all entities for which a tag has a specific value
 *\param IteratorType : an iterator that has map behavior:
 *                      the value of 'first' is the entity handle.
 *                      the value of 'second' is a pointer to the tag data.
 *\param ContainerType : std::vector<MBEntityHandle> or MBRange
 */
template <class IteratorType, class ContainerType>
void find_tag_values_equal( const TagInfo& tag_info,
                            const void* value,
                            int size,
                            IteratorType begin,
                            IteratorType end,
                            ContainerType& results )
{
  switch (tag_info.get_data_type()) {
    case MB_TYPE_INTEGER:
      switch (tag_info.get_size()) {
        case MB_VARIABLE_LENGTH:
          find_tag_values( TagVarIntsEqual( value, size ), begin, end, results );
          break;
        case sizeof(int):
          find_tag_values( TagOneIntEqual( value ), begin, end, results );
          break;
        default:
          find_tag_values( TagIntsEqual( value, size ), begin, end, results );
          break;
      }
      break;
        
    case MB_TYPE_DOUBLE:
      switch (tag_info.get_size()) {
        case MB_VARIABLE_LENGTH:
          find_tag_values( TagVarDoublesEqual( value, size ), begin, end, results );
          break;
        case sizeof(double):
          find_tag_values( TagOneDoubleEqual( value ), begin, end, results );
          break;
        default:
          find_tag_values( TagDoublesEqual( value, size ), begin, end, results );
          break;
      }
      break;
        
    case MB_TYPE_HANDLE:
      switch (tag_info.get_size()) {
        case MB_VARIABLE_LENGTH:
          find_tag_values( TagVarHandlesEqual( value, size ), begin, end, results );
          break;
        case sizeof(MBEntityHandle):
          find_tag_values( TagOneHandleEqual( value ), begin, end, results );
          break;
        default:
          find_tag_values( TagHandlesEqual( value, size ), begin, end, results );
          break;
      }
      break;
        
    default:
      if (tag_info.get_size() == MB_VARIABLE_LENGTH) 
        find_tag_values( TagVarBytesEqual( value, size ), begin, end, results );
      else
        find_tag_values( TagBytesEqual( value, size ), begin, end, results );
      break;
  }
}

/** Iterator to use in find_tag_values_equal for arrays of data */
class ByteArrayIterator 
{
  private:
    size_t step;
    typedef std::pair<MBEntityHandle, const char*> data_type;
    data_type data;
  public:
    ByteArrayIterator( MBEntityHandle start_handle,
                       const void* data_array,
                       size_t tag_size )
      : step(tag_size),
        data(start_handle, reinterpret_cast<const char*>(data_array))
        
      {}
    ByteArrayIterator( MBEntityHandle start_handle,
                       const void* data_array,
                       const TagInfo& tag_info )
      : step(tag_info.get_size() == MB_VARIABLE_LENGTH ? sizeof(VarLenTag) : tag_info.get_size()),
        data(start_handle, reinterpret_cast<const char*>(data_array))
        {}
    bool operator==( const ByteArrayIterator& other ) const
      { return data.first == other.data.first; }
    bool operator!=( const ByteArrayIterator& other ) const
      { return data.first != other.data.first; }
    ByteArrayIterator& operator++()
      { ++data.first; data.second += step; return *this; }
    ByteArrayIterator operator++(int)
      { ByteArrayIterator result(*this); operator++(); return result; }
    ByteArrayIterator& operator--()
      { --data.first; data.second -= step; return *this; }
    ByteArrayIterator operator--(int)
      { ByteArrayIterator result(*this); operator--(); return result; }
    ByteArrayIterator& operator+=(size_t amt)
      { data.first += amt; data.second += amt*step; return *this; }
    ByteArrayIterator& operator-=(size_t amt)
      { data.first -= amt; data.second -= amt*step; return *this; }
    MBEntityHandle operator-( const ByteArrayIterator& other ) const
      { return data.first - other.data.first; }
    const data_type& operator*() const 
      { return data; }
    const data_type* operator->() const 
      { return &data; }
};

#endif

