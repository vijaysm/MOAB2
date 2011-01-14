#ifndef __vtkMOABMesh_h
#define __vtkMOABMesh_h

#include "vtkUnstructuredGrid.h"

#include "moab/Interface.hpp"
#include "moab/WriteUtilIface.hpp"
#include "moab/Range.hpp"

class vtkIntArray;

using namespace moab;

#include <vector>
#include <string>

class VTK_EXPORT vtkMOABMesh
{
public:
  static vtkMOABMesh *instance(Interface *iface = NULL, bool create_new = true);

  ErrorCode delete_instance();
  
  ErrorCode Update(EntityHandle file_set = 0);

  ErrorCode load_file(const char *file_name, const char *options, EntityHandle &file_set);

  vtkUnstructuredGrid *GetOutput();
  void SetOutput(vtkUnstructuredGrid *ug);

  void file_sets(Range &files);
  
  void file_names(std::vector<std::string> &fnames);

  bool file_loaded(const char *filename);
  
  void Execute();

protected:
  vtkMOABMesh(Interface *impl);
  ~vtkMOABMesh();

private:
  vtkMOABMesh(const vtkMOABMesh&);  // Not implemented.
  void operator=(const vtkMOABMesh&);  // Not implemented.

  static vtkMOABMesh *instance_;
  
  vtkUnstructuredGrid *myUG;

  Interface *mbImpl;
  
  WriteUtilIface *iFace;

  Range fileSets;
                               
  static const int vtk_cell_types[];

  int maxPointId;
  int maxCellId;

  Tag vtkIdTag;
  
  std::vector<std::string> fileNames;

  bool outOfDate;
  
  ErrorCode construct_mesh(EntityHandle file_set);

  ErrorCode create_points_vertices(EntityHandle file_set, Range &verts);
  
  ErrorCode create_elements(EntityHandle file_set);
  
  ErrorCode construct_filters();

  ErrorCode read_tags(EntityHandle file_set);
  
  ErrorCode read_sparse_tags(EntityHandle file_set);

  ErrorCode read_dense_tags(EntityHandle file_set);

  void add_name(vtkUnstructuredGrid *output, const char *prefix,
                const int id);
  
};

inline vtkUnstructuredGrid *vtkMOABMesh::GetOutput() 
{
  return myUG;
}

inline void vtkMOABMesh::SetOutput(vtkUnstructuredGrid *ug) 
{
  myUG = ug;
}

inline void vtkMOABMesh::file_sets(Range &fsets) 
{
  fsets.merge(fileSets);
}

inline void vtkMOABMesh::file_names(std::vector<std::string> &fnames) 
{
  std::copy(fileNames.begin(), fileNames.end(), std::back_inserter(fnames));
}

inline bool vtkMOABMesh::file_loaded(const char *filename) 
{
  return (std::find(fileNames.begin(), fileNames.end(), std::string(filename)) != fileNames.end());
}

#endif


