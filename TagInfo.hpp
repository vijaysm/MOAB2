#ifndef TAG_INFO_HPP
#define TAG_INFO_HPP

#include "MBTypes.h"

#include <string>
#include <string.h>  /* memcpy */
#include <stdlib.h>  /* realloc & free */

// ! stores information about a tag
class TagInfo
{
public:

  //! constructor
  TagInfo() : mDefaultValue(0),
              mMeshValue(0),
              mDefaultValueSize(0),
              mMeshValueSize(0),
              mDataSize(0), 
              dataType(MB_TYPE_OPAQUE),
              mTagName(""), 
              isValid(false)
              {}

  //! constructor that takes all parameters
  inline TagInfo( const char * name, 
                  int size, 
                  MBDataType type, 
                  const void * default_value,
                  int default_value_size);
  
  TagInfo( const TagInfo& copy );
  
  inline ~TagInfo();
  
  TagInfo& operator=( const TagInfo& copy );
  
  //! set the name of the tag
  void set_name( const std::string& name) { mTagName = name; }

  //! set the size of the data, might be MB_VARIABLE_LENGTH 
  void set_size( const int size ) { mDataSize = size; }

  //! get the name of the tag
  const std::string& get_name() const { return mTagName; }

  //! get the size of the data
  int get_size() const { return mDataSize; }

    //! get length of default value
  int default_value_size() const { return mDefaultValueSize; }

    //! get the default data
  const void *default_value() const  
    { return mDefaultValue; }
  
    //! set mesh value
  void set_mesh_value( const void* data, int size );
  
    //! get mesh value
  int get_mesh_value_size() const { return mMeshValueSize; }
  
    //! get mesh value
  const void* get_mesh_value() const 
    { return mMeshValue; }
  
    //! remove mesh value
  void remove_mesh_value();
  
  inline MBDataType get_data_type() const     { return dataType; }
  
  inline void set_data_type( MBDataType t )   { dataType = t; }

  static int size_from_data_type( MBDataType t );
  
  bool is_valid() const { return isValid; }
  void invalidate();
  
    // Check that all lengths are valid multiples of the type size.
    // Returns true if all lengths are valid, false othersize.
  bool check_valid_sizes( const int* sizes, int num_sizes ) const;

private:    

  MBErrorCode reserve_mesh_tag_id( int& id_out );

  //! stores the default data, if any
  void* mDefaultValue;
  
  //! store the mesh value, if any
  void* mMeshValue;
  
  //! Size of mDefaultValue and mMeshValue, in bytes
  //! NOTE: These sizes differ from mDataSize in two cases:
  //!    a) Variable-length tags
  //!    b) Bit tags (where mDataSize is bits, not bytes.)
  int mDefaultValueSize, mMeshValueSize;

  //! stores the size of the data for this tag
  int mDataSize;
  
  //! type of tag data
  MBDataType dataType;

  //! stores the tag name
  std::string mTagName;
  
  //! flag to mark unused entries
  bool isValid;

};

inline TagInfo::TagInfo( const char* name, 
                         int size, 
                         MBDataType type,
                         const void* default_value,
                         int default_value_size)
 : mDefaultValue(0),
   mMeshValue(0),
   mDefaultValueSize(default_value_size),
   mMeshValueSize(0),
   mDataSize(size),
   dataType(type),
   mTagName(name),
   isValid(true)
{
  if (default_value) {
    mDefaultValue = malloc( mDefaultValueSize );
    memcpy( mDefaultValue, default_value, mDefaultValueSize );
  }
}

inline TagInfo::~TagInfo() 
{
  free( mDefaultValue );
  free( mMeshValue );
  mDefaultValue = mMeshValue = 0;
  mDefaultValueSize = mMeshValueSize = 0;
}

#endif
