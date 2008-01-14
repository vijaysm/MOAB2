#ifndef TAG_INFO_HPP
#define TAG_INFO_HPP

#include "MBTypes.h"

#include <string>
#include <string.h>

// ! stores information about a tag
class TagInfo
{
public:

  //! constructor
  TagInfo() : mTagName(""), 
              mDataSize(0), 
              isValid(false),
              mDefaultValue(NULL), 
              mMeshValue(NULL),
              dataType(MB_TYPE_OPAQUE)
              {}
  
  //! destructor
  inline ~TagInfo();
  
  //! copy constructor
  inline TagInfo(const TagInfo&);

  //! constructor that takes all parameters
  inline TagInfo(const char *, int, MBDataType type, const void *);

  //! assignment operator
  inline TagInfo &operator=(const TagInfo &rhs);
  
  //! set the name of the tag
  void set_name( const std::string& name) { mTagName = name; }

  //! set the size of the data
  void set_size( const int size ) { mDataSize = size; }

  //! get the name of the tag
  const std::string& get_name() const { return mTagName; }

  //! get the size of the data
  int get_size() const { return mDataSize; }

    //! get the default data
  const void *default_value() const  { return mDefaultValue;}
  
    //! set mesh value
  void set_mesh_value( const void* data );
  
    //! get mesh value
  const void* get_mesh_value() const { return mMeshValue; }
  
    //! remove mesh value
  void remove_mesh_value();
  
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
  unsigned char *mDefaultValue;
  
  //! store the mesh value, if any
  unsigned char *mMeshValue;
  
  //! type of tag data
  MBDataType dataType;

};


inline TagInfo::TagInfo(const TagInfo& copy)
  : mTagName( copy.mTagName ),
    mDataSize( copy.mDataSize ),
    isValid( copy.isValid ),
    mDefaultValue( 0 ),
    mMeshValue( 0 ),
    dataType( copy.dataType )
{
  if (copy.mDefaultValue) {
    mDefaultValue = new unsigned char[mDataSize];
    memcpy(mDefaultValue, copy.mDefaultValue, mDataSize);
  }
  
  if (copy.mMeshValue) {
    mMeshValue = new unsigned char[mDataSize];
    memcpy(mMeshValue, copy.mMeshValue, mDataSize);
  }
}

inline TagInfo::TagInfo( const char* name, 
                         int size, 
                         MBDataType type,
                         const void* default_value)
 : mTagName( name ),
   mDataSize( size ),
   isValid( true ),
   mDefaultValue( 0 ),
   mMeshValue( 0 ),
   dataType( type )
{
  if (default_value) {
    mDefaultValue = new unsigned char[size];
    memcpy(mDefaultValue, default_value, size);
  }
}

inline TagInfo &TagInfo::operator=(const TagInfo &rhs)
{
  mTagName = rhs.mTagName;
  mDataSize = rhs.mDataSize;
  isValid = rhs.isValid;
  
  delete [] mDefaultValue;
  delete [] mMeshValue;
  mDefaultValue = 0;
  mMeshValue = 0;
  
  if (rhs.mDefaultValue) {
    mDefaultValue = new unsigned char[mDataSize];
    memcpy( mDefaultValue, rhs.mDefaultValue, mDataSize );
  }
  
  if (rhs.mMeshValue) {
    mMeshValue = new unsigned char[mDataSize];
    memcpy( mMeshValue, rhs.mMeshValue, mDataSize );
  }
  
  return *this;
}

inline TagInfo::~TagInfo() 
{
  // clean up default value
  delete [] mDefaultValue;
  delete [] mMeshValue;
}


inline void TagInfo::set_mesh_value( const void* data )
{
  if (!mMeshValue)
    mMeshValue = new unsigned char[mDataSize];
  memcpy( mMeshValue, data, mDataSize );
}

inline void TagInfo::remove_mesh_value()
{ 
  delete [] mMeshValue;
  mMeshValue = 0;
}

#endif
