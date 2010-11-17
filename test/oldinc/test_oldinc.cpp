  // test forrwards
#include "MBForward.hpp"  
typedef MBRange Foo1;
typedef MBInterface Foo2;
typedef MBProcConfig Foo3;

  // test types
#include "MBEntityHandle.h"
MBEntityHandle handle = 0;
#include "MBTypes.h"
MBEntityType type = MBVERTEX;
MBTag tag = 0;
MBErrorCode code = MB_SUCCESS;
int i = MB_VARIABLE_LENGTH;
MBTagType tagtype = MB_TAG_DENSE;
MBDataType datatype = MB_TYPE_HANDLE;
MBEntitySetProperty prop = MESHSET_SET;

  // test version info
#include "MBVersion.h"

#ifndef MB_VERSION
# error "MB_VERSION not defined"
#endif

#ifndef MB_VERSION_MAJOR
# error "MB_VERSION_MAJOR not defined"
#endif

#ifndef MB_VERSION_MINOR
# error "MB_VERSION_MINOR not defined"
#endif

#ifndef MB_VERSION_STRING
# error "MB_VERSION_STRING not defined"
#endif

#include "MBCN.hpp"

#include "MBUnknownInterface.hpp"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBUtil.hpp"

#ifdef USE_MPI
#  include "MBmpi.h"
#  include "MBParallelData.hpp"
#  include "MBParallelComm.hpp"
#  include "MBProcConfig.hpp"
#endif

#include "MBCore.hpp"
#include "MBReaderIface.hpp"
#include "MBWriterIface.hpp"
#include "MBReadUtilIface.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBReaderWriterSet.hpp"

#include "MBGeomUtil.hpp"
#include "MBAdaptiveKDTree.hpp"
#include "MBBSPTree.hpp"
#include "MBOrientedBoxTreeTool.hpp"
#include "MBSkinner.hpp"

#include "MBCartVect.hpp"
#include "MBBSPTreePoly.hpp"

int main()
{
    // check that the expected types are defined
  MBCore mb_core;
  MBInterface& mb = mb_core;
  MBRange range;
    // initialize and test the following to NULL to eliminate compiler warnings
  MBReaderIface* read_ptr = NULL;
  MBWriterIface* write_ptr = NULL;
  MBReadUtilIface* read_util = NULL;
  MBWriteUtilIface* write_util = NULL;
  MBReaderWriterSet* io_set_ptr = NULL;
  if (read_ptr || write_ptr || read_util || write_util || io_set_ptr) ;
  
  MBAdaptiveKDTree kdtree_tool(&mb);
  MBBSPTree bsptree_tool(&mb);
  MBOrientedBoxTreeTool obbtree_tool(&mb);
  MBSkinner skin_tool(&mb);

#ifdef WITH_MPI
  MBParallelComm* pcomm_ptr;
  MBProcConfig* pconf_ptr;
#endif

  MBCartVect vect(0.0);
  MBBSPTreePoly poly;

  return 0;
}

  
  
