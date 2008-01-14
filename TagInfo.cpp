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
  delete [] mDefaultValue;
  mDefaultValue = 0;
  delete [] mMeshValue;
  mMeshValue = 0;
}
