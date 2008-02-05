#include "TagInfo.hpp"


int TagInfo::size_from_data_type( MBDataType t )
{
  static const int sizes[] = { 1, 
                               sizeof(int), 
                               sizeof(double), 
                               1, 
                               sizeof(MBEntityHandle),
                               0 };  
   return sizes[t];
}


void TagInfo::invalidate()
{
  mTagName.clear();
  isValid = false;
  free( mMeshValue );
  free( mDefaultValue );
  mDefaultValue = mMeshValue = 0;
  mDefaultValueSize = mMeshValueSize = 0;
}

  
TagInfo::TagInfo( const TagInfo& copy )
  : mDefaultValue(0),
    mMeshValue(0),
    mDefaultValueSize(copy.mDefaultValueSize),
    mMeshValueSize(copy.mMeshValueSize),
    mDataSize(copy.mDataSize),
    dataType(copy.dataType),
    mTagName(copy.mTagName),
    isValid(copy.isValid)
{
  if (mDefaultValueSize) {
    mDefaultValue = malloc( mDefaultValueSize );
    memcpy( mDefaultValue, copy.mDefaultValue, mDefaultValueSize );
  }
  if (mMeshValueSize) {
    mMeshValue = malloc( mMeshValueSize );
    memcpy( mMeshValue, copy.mMeshValue, mMeshValueSize );
  }
}

TagInfo& TagInfo::operator=( const TagInfo& copy )
{
  if (copy.mDefaultValue) {
    if (mDefaultValueSize != copy.mDefaultValueSize)
      mDefaultValue = realloc( mDefaultValue,copy.mDefaultValueSize);
    mDefaultValueSize = copy.mDefaultValueSize;
    memcpy( mDefaultValue, copy.mDefaultValue, copy.mDefaultValueSize );
  }
  else if (mDefaultValue) {
    free( mDefaultValue );
    mDefaultValue = 0;
    mDefaultValueSize = 0;
  }

  if (copy.mMeshValue) {
    if (mMeshValueSize != copy.mMeshValueSize)
      mMeshValue = realloc( mMeshValue,copy.mMeshValueSize);
    mMeshValueSize = copy.mMeshValueSize;
    memcpy( mMeshValue, copy.mMeshValue, copy.mMeshValueSize );
  }
  else if (mMeshValue) {
    free( mMeshValue );
    mMeshValue = 0;
    mMeshValueSize = 0;
  }
  
  mDataSize = copy.mDataSize;
  dataType = copy.dataType;
  mTagName = copy.mTagName;
  isValid = copy.isValid;
  return *this;
}

void TagInfo::set_mesh_value( const void* data, int size )
{
  if (mMeshValueSize != size) {
    mMeshValueSize = size;
    mMeshValue = realloc( mMeshValue, size );
  }
  memcpy( mMeshValue, data, size );
}

    //! remove mesh value
void TagInfo::remove_mesh_value() 
{
  free( mMeshValue );
  mMeshValueSize = 0;
}

  
    // Check that all lengths are valid multiples of the type size.
    // Returns true if all lengths are valid, false othersize.
bool TagInfo::check_valid_sizes( const int* sizes, int num_sizes ) const
{
  unsigned sum = 0;
  const unsigned size = size_from_data_type( get_data_type() );
  if (size == 1)
    return true;
  for (int i = 0; i < num_sizes; ++i)
    sum |= ((unsigned)sizes[i]) % size;
  return (sum == 0);
}

