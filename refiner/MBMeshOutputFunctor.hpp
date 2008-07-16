#ifndef MB_MESHOUTPUTFUNCTOR_HPP
#define MB_MESHOUTPUTFUNCTOR_HPP

#include "MBTypes.h"
#include "MBEntityRefiner.hpp"
#include "MBProcessSet.hpp"

#include <vector>
#include <map>

#include <string.h>

class MBSplitVerticesBase;
class MBParallelComm;

class MBMeshOutputFunctor : public MBEntityRefinerOutputFunctor
{
public:
  MBMeshOutputFunctor( MBRefinerTagManager* tag_mgr );
  ~MBMeshOutputFunctor();

  void print_vert_crud( MBEntityHandle vout, int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags );
  void assign_global_ids( MBParallelComm* comm );

  void assign_tags( MBEntityHandle vhandle, const void* vtags );

  virtual MBEntityHandle operator () ( MBEntityHandle vhash, const double* vcoords, const void* vtags );
  virtual MBEntityHandle operator () ( int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags );
  virtual void operator () ( MBEntityHandle h );
  virtual void operator () ( MBEntityType etyp );

  MBInterface* mesh_in;
  MBInterface* mesh_out;
  bool input_is_output;
  std::vector<MBSplitVerticesBase*> split_vertices;
  std::vector<MBEntityHandle> elem_vert;
  MBRefinerTagManager* tag_manager;
  MBEntityHandle destination_set;
  std::map<MBProcessSet,int> proc_partition_counts;
};

#endif // MB_MESHOUTPUTFUNCTOR_HPP
