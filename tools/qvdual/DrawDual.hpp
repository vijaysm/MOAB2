#ifndef DRAWDUAL_HPP
#define DRAWDUAL_HPP

class Agnode_t;
class Agedge_t;
class Agraph_t;
class QVTKWidget;

#include "MBInterface.hpp"
#include "MBRange.hpp"
#include <map>
#include <vector>

class vtkPolyData;
class vtkRenderer;
class vtkCellPicker;
class vtkUnstructuredGrid;
class vtkExtractCells;
class vtkObject;
class vtkFloatArray;
class vtkCellArray;
class vtkActor;
class Agsym_t;

class DrawDual
{
public:
  DrawDual();
  ~DrawDual();

  bool draw_dual_surfs(MBRange &dual_surfs);
  MBErrorCode draw_dual_surf(MBEntityHandle dual_surf);
  
private:

  static DrawDual *gDrawDual;
  
  class GVEntity
  {
  public:
    int numGvizEntities;
    MBEntityHandle dualSurfs[3];
    MBEntityHandle moabEntity;
    int pointPos[3][2];
    int vtkEntityIds[3];
    vtkActor *myActors[3];
    Agnode_t *gvizPoints[3];
    Agedge_t *gvizEdges[3];

    GVEntity() 
      {
        numGvizEntities = 0;
        dualSurfs[0] = dualSurfs[1] = dualSurfs[2] = 0;
        moabEntity = 0;
        pointPos[0][0] = pointPos[0][1] = pointPos[0][2] = 
          pointPos[1][0] = pointPos[1][1] = pointPos[1][2] = 0;
        vtkEntityIds[0] = vtkEntityIds[1] = vtkEntityIds[2] = -1;
        myActors[0] = myActors[1] = myActors[2] = NULL;
        gvizPoints[0] = gvizPoints[1] = gvizPoints[2] = NULL;
        gvizEdges[0] = gvizEdges[1] = gvizEdges[2] = NULL;
      }
    int get_index(const MBEntityHandle dual_surf) 
      {
        if (dual_surf == dualSurfs[0]) return 0;
        else if (dual_surf == dualSurfs[1]) return 1;
        else if (dual_surf == dualSurfs[2]) return 2;
        else return -1;
      }
    int add_gvpoint(const MBEntityHandle dual_surf, Agnode_t *this_pt) 
      {
        int surf_num = get_index(0);
        if (surf_num != -1) {
          dualSurfs[surf_num] = dual_surf;
          gvizPoints[surf_num] = this_pt;
        }
        return surf_num;
      }
    int add_gvedge(const MBEntityHandle dual_surf, Agedge_t *this_ed) 
      {
        int surf_num = get_index(0);
        if (surf_num != -1) {
          dualSurfs[surf_num] = dual_surf;
          gvizEdges[surf_num] = this_ed;
        }
        return surf_num;
      }
  };

  class GraphWindows 
  {
  public:
    Agraph_t *gvizGraph;
    QVTKWidget *qvtkWidget;
    vtkExtractCells *pickExtractor;
    vtkPolyData *highPoly;

    GraphWindows() : gvizGraph(NULL), qvtkWidget(NULL), pickExtractor(NULL),
      highPoly(NULL) {}
  };
  
    //! make sure all dual vertices and edges have graphviz nodes and edges
  MBErrorCode construct_graphviz_data(MBEntityHandle dual_surf);
  
    //! given the loop vertices, compute and fix their points
  MBErrorCode compute_fixed_points(MBEntityHandle dual_surf, MBRange &dverts,
                            MBRange &face_verts);

    //! compute the position on the loop, accounting for multiple loops
  void get_loop_vertex_pos(unsigned int vert_num, 
                           unsigned int loop_num, 
                           unsigned int num_loops, 
                           double angle, int &xpos_pts, int &ypos_pts);
  
    //! construct the points & cells for the vtkPolyData from the MOAB/graphviz data
  MBErrorCode make_vtk_data(MBEntityHandle dual_surf,
                            vtkPolyData *pd,
                            vtkRenderer *this_ren);

    //! construct dim-dimensional cells
  MBErrorCode make_vtk_cells(const MBRange &cell_range, const int dim,
                             const float color_index,
                             const MBEntityHandle dual_surf,
                             std::map<MBEntityHandle, GVEntity *> &vert_gv_map,
                             vtkPolyData *pd,
                             vtkFloatArray *color_ids);
  
    //! given a qvtk widget, return the first polydata supplying data to it
  vtkPolyData *get_polydata(QVTKWidget *this_wid);
  
    //! get a clean polydata for this widget
  void get_clean_pd(MBEntityHandle dual_surf,
                    QVTKWidget *&this_wid, vtkPolyData *&pd);

    //! draw various labels with the sheet
  MBErrorCode draw_labels(MBEntityHandle dual_surf,
                          vtkPolyData *pd,
                          vtkPolyData *new_pd);

  MBErrorCode draw_chords(MBEntityHandle dual_surf,
                          vtkPolyData *pd,
                          vtkPolyData *&new_pd);
  
  void label_interior_verts(MBEntityHandle dual_surf,
                            vtkPolyData *pd,
                            vtkRenderer *ren);
  
  MBEntityHandle other_sheet(const MBEntityHandle this_chord,
                             const MBEntityHandle dual_surf);
  
  MBErrorCode get_primal_ids(const MBRange &ents, std::vector<int> &ids);
  
  MBErrorCode get_dual_entities(const MBEntityHandle dual_ent,
                                MBRange *dcells,
                                MBRange *dedges = NULL,
                                MBRange *dverts = NULL,
                                MBRange *dverts_loop = NULL);
  
  MBErrorCode allocate_points(MBEntityHandle dual_surf,
                              vtkPolyData *pd,
                              MBRange &verts,
                              std::map<MBEntityHandle, GVEntity*> &vert_gv_map);
  
  static void add_picker(vtkRenderer *this_ren);
  
  static void process_events(vtkObject *caller, 
                             unsigned long event,
                             void* clientdata, 
                             void* /*vtkNotUsed(calldata)*/);
  
  static void process_pick();
  
    //! map of dual surfaces and windows they're drawn in
  std::map<MBEntityHandle, GraphWindows> surfDrawrings;
  
    //! cache some of the tags we use
  MBTag gvEntityHandle, dualEntityTagHandle;

  static MBTag dualSurfaceTagHandle, dualCurveTagHandle;

    //! gviz graphics context, seems to be needed for layout
  void *gvizGvc;

    //! information about the sheet drawing window
  int xSize, ySize, xOrigin, yOrigin;

    //! picker for dual data
  static vtkCellPicker *dualPicker;

    //! entities which are currently picked
  static MBRange pickRange;

  MBEntityHandle get_picked_cell(MBEntityHandle cell_set,
                                 const int dim,
                                 const int picked_cell);

  void update_high_polydatas();
  
  MBErrorCode get_xform(MBEntityHandle dual_surf, Agsym_t *asym_pos, 
                        double &x_xform, double &y_xform);
  
};


#endif
