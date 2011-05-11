#ifndef IGEOM_MOAB_HPP
#define IGEOM_MOAB_HPP

#include "iGeom.h"
#include "moab/Forward.hpp"
#include "iMesh.h"
#include "MBiMesh.hpp"

namespace moab { class GeomTopoTool; }

/* map from MB's entity type to TSTT's entity type */
extern const iBase_EntityType tstt_type_table[moab::MBMAXTYPE+1];

/* map to MB's entity type from TSTT's entity topology */
extern const moab::EntityType mb_topology_table[moab::MBMAXTYPE+1];

/* map from TSTT's tag types to MOAB's */
extern const moab::DataType mb_data_type_table[4];

/* map from MOAB's tag types to tstt's */
extern const iBase_TagValueType tstt_data_type_table[moab::MB_MAX_DATA_TYPE+1];

/* map from MOAB's MBErrorCode to tstt's */
extern "C" const iBase_ErrorType iBase_ERROR_MAP[moab::MB_FAILURE+1];

/* Create ITAPS iterator */
iBase_EntityIterator create_itaps_iterator( moab::Range& swap_range,
                                            int array_size = 1 ); 

/* Define macro for quick reference to MBInterface instance */
static inline moab::Interface* MBI_cast( iGeom_Instance i )  
  { return reinterpret_cast<MBiMesh*>(i)->mbImpl; }
#define MBI MBI_cast(instance)

/* Define macro for quick reference to moab::Interface instance */
static inline moab::EntityHandle MBH_cast( iBase_EntityHandle h )  
  { return reinterpret_cast<moab::EntityHandle>(h); }         

#define MIN(a,b) (a > b ? b : a)

#define IMESH_INSTANCE(a) reinterpret_cast<iMesh_Instance>(a)

#define GETGTT(a) {if (_my_geomTopoTool == NULL) _my_geomTopoTool =	\
	new GeomTopoTool(reinterpret_cast<MBiMesh*>(a)->mbImpl);}

static inline bool iGeom_isError(int code)
  { return (iBase_SUCCESS != code); }
static inline bool iGeom_isError(moab::ErrorCode code)
  { return (MB_SUCCESS != code); }

#define MBIGEOMI mbimeshi_instance(IMESH_INSTANCE(instance))

#define RETURN(CODE)                                                   \
  do {                                                                 \
    *err = MBIGEOMI->set_last_error((CODE), "");                       \
    return;                                                            \
  } while(false)

#define ERROR(CODE,MSG)                                                \
  do {                                                                 \
    *err = MBIGEOMI->set_last_error((CODE), (MSG));                    \
    return;                                                            \
  } while(false)

#define CHKERR(CODE,MSG)                                               \
  do {                                                                 \
    if (iGeom_isError((CODE)))                                         \
      ERROR((CODE),(MSG));                                             \
  } while(false)

#define FWDERR()                                                       \
  do {                                                                 \
    if (iGeom_isError(*err))                                           \
      return;                                                          \
  } while(false)

#define CHECK_SIZE(array, allocated, size, type, retval)               \
  do {                                                                 \
    if (0 != allocated && NULL != array && allocated < (size)) {       \
      ERROR(iBase_MEMORY_ALLOCATION_FAILED, "Allocated array not "     \
            "enough to hold returned contents.");                      \
    }                                                                  \
    if ((size) && ((allocated) == 0 || NULL == (array))) {             \
      array = (type*)malloc((size)*sizeof(type));                      \
      allocated=(size);                                                \
      if (NULL == array) {                                             \
        ERROR(iBase_MEMORY_ALLOCATION_FAILED,                          \
              "Couldn't allocate array.");                             \
      }                                                                \
    }                                                                  \
  } while(false)

#endif // IGEOM_MOAB_HPP
