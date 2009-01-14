#ifndef IMESH_MOAB_HPP
#define IMESH_MOAB_HPP

#include "iMesh.h"
#include "MBForward.hpp"

/* map from MB's entity type to TSTT's entity topology */
extern const iMesh_EntityTopology tstt_topology_table[MBMAXTYPE+1];

/* map from MB's entity type to TSTT's entity type */
extern const iBase_EntityType tstt_type_table[MBMAXTYPE+1];

/* map to MB's entity type from TSTT's entity topology */
extern const MBEntityType mb_topology_table[MBMAXTYPE+1];

/* map from TSTT's tag types to MOAB's */
extern const MBDataType mb_data_type_table[4];

/* map from MOAB's tag types to tstt's */
extern const iBase_TagValueType tstt_data_type_table[MB_MAX_DATA_TYPE+1];

/* map from MOAB's MBErrorCode to tstt's */
extern const iBase_ErrorType iBase_ERROR_MAP[MB_FAILURE+1];

/* Create ITAPS iterator */
iMesh_EntityIterator create_itaps_iterator( MBRange& swap_range,
                                            int array_size = 1 ); 

/* Define macro for quick reference to MBInterface instance */
static inline MBInterface* MBI_cast( iMesh_Instance i )  
  { return reinterpret_cast<MBInterface*>(i); }         
#define MBI MBI_cast(instance)

/* Most recently returned error code */
extern "C" iBase_Error iMesh_LAST_ERROR;

#define RETURN(a) do {iMesh_LAST_ERROR.error_type = *err = (a); return;} while(false)

#define MBRTN(a) RETURN(iBase_ERROR_MAP[(a)])

#define CHKERR(err) if (MB_SUCCESS != (err)) MBRTN(err)

#include "MBCore.hpp"

class MBiMesh : public MBCore
{
private:
  bool haveDeletedEntities;
  bool fullConnectivity;
public:
  MBiMesh(bool adj_includes_ho = false);

  virtual ~MBiMesh();
  bool have_deleted_ents( bool reset ) {
    bool result = haveDeletedEntities;
    if (reset)
      haveDeletedEntities = false;
    return result;
  }

  virtual MBErrorCode delete_mesh();
  virtual MBErrorCode delete_entities( const MBEntityHandle*, const int );
  virtual MBErrorCode delete_entities( const MBRange& );
  int AdjTable[16];
  bool adj_includes_ho_nodes() const { return fullConnectivity; }
};

#define MBimesh reinterpret_cast<MBiMesh*>(MBI)

#endif // IMESH_MOAB_HPP
