#ifndef IMESH_MOAB_HPP
#define IMESH_MOAB_HPP

#include "iMesh.h"
#include "moab/Forward.hpp"
#include <cstring>
#include <cstdlib>

using namespace moab;

/* map from MB's entity type to TSTT's entity topology */
extern const iMesh_EntityTopology tstt_topology_table[MBMAXTYPE+1];

/* map from MB's entity type to TSTT's entity type */
extern const iBase_EntityType tstt_type_table[MBMAXTYPE+1];

/* map to MB's entity type from TSTT's entity topology */
extern const EntityType mb_topology_table[MBMAXTYPE+1];

/* map from TSTT's tag types to MOAB's */
extern const DataType mb_data_type_table[4];

/* map from MOAB's tag types to tstt's */
extern const iBase_TagValueType tstt_data_type_table[MB_MAX_DATA_TYPE+1];

/* map from MOAB's ErrorCode to tstt's */
extern "C" const iBase_ErrorType iBase_ERROR_MAP[MB_FAILURE+1];

/* Create ITAPS iterator */
iMesh_EntityIterator create_itaps_iterator( Range& swap_range,
                                            int array_size = 1 ); 

/* Define macro for quick reference to Interface instance */
static inline Interface* MBI_cast( iMesh_Instance i )  
  { return reinterpret_cast<Interface*>(i); }         
#define MBI MBI_cast(instance)

/* Most recently returned error code */
extern "C" iBase_Error iMesh_LAST_ERROR;

#define RETURN(a) do {iMesh_LAST_ERROR.error_type = *err = (a); \
                      iMesh_LAST_ERROR.description[0] = '\0'; \
                      return;} while(false)

#include "moab/Core.hpp"

class MBiMesh : public Core
{
private:
  bool haveDeletedEntities;
public:
  MBiMesh();

  virtual ~MBiMesh();
  bool have_deleted_ents( bool reset ) {
    bool result = haveDeletedEntities;
    if (reset)
      haveDeletedEntities = false;
    return result;
  }

  virtual ErrorCode delete_mesh();
  virtual ErrorCode delete_entities( const EntityHandle*, const int );
  virtual ErrorCode delete_entities( const Range& );
  int AdjTable[16];
};

#define MBimesh reinterpret_cast<MBiMesh*>(MBI)


static inline int
iMesh_processError( int code, const char* desc ) 
{
  std::strncpy( iMesh_LAST_ERROR.description, desc,
                sizeof(iMesh_LAST_ERROR.description) );
  iMesh_LAST_ERROR.description[sizeof(iMesh_LAST_ERROR.description)-1] = '\0';
  return (iMesh_LAST_ERROR.error_type = (iBase_ErrorType)code);
}

#define ERROR(CODE,MSG) do { *err = iMesh_setLastError( MBI, (CODE), (MSG) ); return; } while(false)
#define IBASE_ERROR(CODE,MSG) do { *err = iMesh_processError( (CODE), (MSG) ); return; } while(false)

static inline int iMesh_setLastError( Interface*, int code, const char* msg )
  { return iMesh_processError( code, msg ); }  
static inline int iMesh_setLastError( Interface* mbi, ErrorCode code, const char* msg )
  { 
    std::string message(msg);
    message += "  (MOAB Error Code: ";
    message += mbi->get_error_string(code);
    message += ")";
    return iMesh_processError( iBase_ERROR_MAP[code], message.c_str() ); 
  }

#define CHKERR(CODE,MSG) \
  if (iMesh_isError((CODE))) ERROR((CODE),(MSG))

static inline bool iMesh_isError(int code)
  { return (iBase_SUCCESS != code); }
static inline bool iMesh_isError(ErrorCode code)
  { return (MB_SUCCESS != code); }


// Check the array size, and allocate the array if necessary.
// Free the array upon leaving scope unless KEEP_ARRAY
// is invoked.
#define ALLOC_CHECK_ARRAY(array, this_size) \
  iMeshArrayManager array ## _manager ( reinterpret_cast<void**>(array), *(array ## _allocated), *(array ## _size), this_size, sizeof(**array), err ); \
  if (iBase_SUCCESS != *err) return

#define ALLOC_CHECK_TAG_ARRAY(array, this_size) \
  iMeshArrayManager array ## _manager ( reinterpret_cast<void**>(array), *(array ## _allocated), *(array ## _size), this_size, 1, err ); \
  if (iBase_SUCCESS != *err) return

#define KEEP_ARRAY(array) \
  array ## _manager .keep_array()

// Check the array size, and allocate the array if necessary.
// Do NOT free the array upon leaving scope.
#define ALLOC_CHECK_ARRAY_NOFAIL(array, this_size) \
  ALLOC_CHECK_ARRAY(array, this_size); KEEP_ARRAY(array)


// Implement RAII pattern for allocated arrays
class iMeshArrayManager
{
  void** arrayPtr;

public:


  iMeshArrayManager( void** array_ptr,
                     int& array_allocated_space,
                     int& array_size,
                     int count,
                     int val_size,
                     int* err ) : arrayPtr(0)
  {
    if (!*array_ptr) {
      *array_ptr = std::malloc(val_size * count);
      array_allocated_space = array_size = count;
      if (!*array_ptr) {
        IBASE_ERROR(iBase_MEMORY_ALLOCATION_FAILED, "Couldn't allocate array.");
      }
      arrayPtr = array_ptr;
    }
    else {
      array_size = count;
      if (array_allocated_space < count) {
        IBASE_ERROR(iBase_BAD_ARRAY_DIMENSION, 
          "Allocated array not large enough to hold returned contents.");
      }
    }
    RETURN(iBase_SUCCESS);
  }
  
  ~iMeshArrayManager() 
  {
    if (arrayPtr) {
      std::free(*arrayPtr);
      *arrayPtr = 0;
    }
  }
  
  void keep_array()
    { arrayPtr = 0; }
};

#endif // IMESH_MOAB_HPP
