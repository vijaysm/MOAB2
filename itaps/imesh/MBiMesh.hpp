#ifndef MBIMESH_HPP
#define MBIMESH_HPP

#include "moab/Core.hpp"
using namespace moab;

class MBiMesh
{
private:
  bool haveDeletedEntities;
  bool iCreatedInterface;
public:
  MBiMesh(moab::Interface *mbImpl = NULL);

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
  moab::Interface *mbImpl;
};

static inline MBiMesh *mbimeshi_instance(iMesh_Instance instance) {return reinterpret_cast<MBiMesh*>(instance);}
#define MBIMESHI mbimeshi_instance(instance)
#define MOABI MBIMESHI->mbImpl

inline MBiMesh::MBiMesh(Interface *impl)
        : haveDeletedEntities(false), iCreatedInterface(false), mbImpl(impl)
{
  int tmp_table[] = {
      1, 1, 1, 1,
      1, 0, 2, 2,
      1, 2, 0, 2,
      1, 2, 2, 1
  };
  memcpy(AdjTable, tmp_table, 16*sizeof(int));

  if (!mbImpl) {
    mbImpl = new Core();
    iCreatedInterface = true;
  }
}

inline MBiMesh::~MBiMesh() 
{
  if (iCreatedInterface) delete mbImpl;
}

inline ErrorCode MBiMesh::delete_mesh() {
  haveDeletedEntities = true;
  return mbImpl->delete_mesh();
}

inline ErrorCode MBiMesh::delete_entities( const EntityHandle* a, const int n )
{
  if (n > 0)
    haveDeletedEntities = true;
  return mbImpl->delete_entities( a, n );
}

inline ErrorCode MBiMesh::delete_entities( const Range& r )
{
  if (!r.empty())
    haveDeletedEntities = true;
  return mbImpl->delete_entities( r );
}

#endif
