#include <iostream>
#include <fstream>
#include <vector>

#include "MBReaderIface.hpp"
#include "MBInterface.hpp"

#define VERTEX_LIST       2411 
#define MAKE_TETRAHEDRA   2412

#define MAT_PROP_TABLE_TAG "mat_prop_table"
#define PHYS_PROP_TABLE_TAG "phys_prop_table"

class MBReadUtilIface;

class ReadIDEAS : public MBReaderIface
{

public:

  static MBReaderIface* factory( MBInterface* );

  MBErrorCode load_file(const char* fname, 
			MBEntityHandle& meshset, 
			const FileOptions&,
			const int* material_set_list,
			int num_material_sets );
  
  //! Constructor
  ReadIDEAS(MBInterface* impl = NULL);
  
  //! Destructor
  virtual ~ReadIDEAS() {}

protected:
  
  void skip_header();
  void create_vertices();
  void create_tetrahedral_elements();
  
private:
  
  std::ifstream file;
  
  // Read mesh interface
  MBReadUtilIface* readMeshIface;
  
  // MOAB Interface
  MBInterface* MBI;
  
  // Handle for the mesh
  MBEntityHandle mesh_handle;

};
