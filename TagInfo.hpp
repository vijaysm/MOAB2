#ifndef TAG_INFO_HPP
#define TAG_INFO_HPP

#include "MBTypes.h"
#include "VarLenTag.hpp"

#include <string>
#include <string.h>
#include <assert.h>

// ! stores information about a tag
class TagInfo
{
public:

  //! constructor
  TagInfo() : mTagName(""), 
              mDataSize(0), 
              isValid(false),
              dataType(MB_TYPE_OPAQUE)
              {}

  //! constructor that takes all parameters
  inline TagInfo( const char * name, 
                  int size, 
                  MBDataType type, 
                  const void * default_value,
                  int default_value_size);
  
  //! set the name of the tag
  void set_name( const std::string& name) { mTagName = name; }

  //! set the size of the data
  void set_size( const int size ) { mDataSize = size; }

  //! get the name of the tag
  const std::string& get_name() const { return mTagName; }

  //! get the size of the data
  int get_size() const { return mDataSize; }

    //! get length of default value
  int default_value_size() const { return mDefaultValue.size(); }

    //! get the default data
  const void *default_value() const  
    { return mDefaultValue.size() ? mDefaultValue.data() : 0;}
  
    //! set mesh value
  void set_mesh_value( const void* data, int size );
  
    //! get mesh value
  int get_mesh_value_size() const { return mMeshValue.size(); }
  
    //! get mesh value
  const void* get_mesh_value() const 
    { return mMeshValue.size() ? mMeshValue.data() : 0; }
  
    //! remove mesh value
  void remove_mesh_value() { mMeshValue.clear(); }
  
  inline MBDataType get_data_type() const     { return dataType; }
  
  inline void set_data_type( MBDataType t )   { dataType = t; }

  static int size_from_data_type( MBDataType t );
  
  bool is_valid() const { return isValid; }
  void invalidate();

private:    

  MBErrorCode reserve_mesh_tag_id( int& id_out );

  //! stores the tag name
  std::string mTagName;

  //! stores the size of the data for this tag
  unsigned short mDataSize;
  
  //! flag to mark unused entries
  bool isValid;

  //! stores the default data, if any
  VarLenTag mDefaultValue;
  
  //! store the mesh value, if any
  VarLenTag mMeshValue;
  
  //! type of tag data
  MBDataType dataType;

};

inline TagInfo::TagInfo( const char* name, 
                         int size, 
                         MBDataType type,
                         const void* default_value,
                         int default_value_size)
 : mTagName( name ),
   mDataSize( size ),
   isValid( true ),
   mDefaultValue( default_value_size, default_value ),
   dataType( type )
{
    // if tag is not variable-length and default_value_size is not zero,
    // then size and default_value_size must be the same.
  assert( size == MB_VARIABLE_LENGTH || default_value_size == 0 || default_value_size == size );
}


inline void TagInfo::set_mesh_value( const void* data, int size )
{
    // if tag is not variable-length, then size must be tag size
  assert( get_size() == MB_VARIABLE_LENGTH || get_size() == size );
  mMeshValue.set( data, size );
}

#endif
