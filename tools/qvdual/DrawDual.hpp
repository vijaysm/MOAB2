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

class DualTool;
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
class QLineEdit;

class DrawDual
{
public:
  DrawDual(QLineEdit *pickline1, QLineEdit *pickline2);
  ~DrawDual();

  bool draw_dual_surfs(MBRange &dual_surfs,
                       const bool use_offsets = false);
  bool draw_dual_surfs(std::vector<MBEntityHandle> &dual_surfs,
                       const bool use_offsets = false);
  MBErrorCode draw_dual_surf(MBEntityHandle dual_surf,
                             int offset_num = 0);
  
  MBEntityHandle lastPickedEnt; // last picked entity
  MBEntityHandle secondLastPickedEnt; // second last picked entity

    //! reset the drawing data for a sheet
  MBErrorCode reset_drawing_data(MBEntityHandle dual_surf);

  MBErrorCode reset_drawn_sheets(MBRange &drawn_sheets);
  
private:

  static DrawDual *gDrawDual;
  DualTool *dualTool;
  QLineEdit *pickLine1, *pickLine2;

  class GVEntity
  {
  public:
    int numGvizEntities;
    MBEntityHandle dualSurfs[3];
    MBEntityHandle moabEntity;
    int pointPos[3][2];
    int vtkEntityIds[4]; // extra pt for edge mid-pts
    vtkActor *myActors[3];
    Agnode_t *gvizPoints[5]; // extra 2 for edge mid-pts
    Agedge_t *gvizEdges[4]; // extra 2 for extra edges
    
    GVEntity() 
      {
        numGvizEntities = 0;
        dualSurfs[0] = dualSurfs[1] = dualSurfs[2] = 0;
        moabEntity = 0;
        pointPos[0][0] = pointPos[0][1] = pointPos[0][2] = 
          pointPos[1][0] = pointPos[1][1] = pointPos[1][2] = 0;
        vtkEntityIds[0] = vtkEntityIds[1] = vtkEntityIds[2] = vtkEntityIds[3] = -1;
        myActors[0] = myActors[1] = myActors[2] = NULL;
        gvizPoints[0] = gvizPoints[1] = gvizPoints[2] = gvizPoints[3] = 
          gvizPoints[4] = NULL;
        gvizEdges[0] = gvizEdges[1] = gvizEdges[2] = gvizEdges[3] = NULL;
      }
    void reset(const int index);
    int get_index(const MBEntityHandle dual_surf) 
      {
        if (dual_surf == dualSurfs[0]) return 0;
        else if (dual_surf == dualSurfs[1]) return 1;
        else if (dual_surf == dualSurfs[2]) return 2;
        else if (dualSurfs[0] == 0) return -1;
        else if (dualSurfs[1] == 0) return -2;
        else if (dualSurfs[2] == 0) return -3;
        else return -10;
      }
  };

  class GraphWindows 
  {
  public:
    Agraph_t *gvizGraph;
    QVTKWidget *qvtkWidget;
    vtkActor *pickActor;

    GraphWindows() : gvizGraph(NULL), qvtkWidget(NULL), pickActor(NULL) {}
  };
  
    //! make sure all dual vertices and edges have graphviz nodes and edges
  MBErrorCode construct_graphviz_data(MBEntityHandle dual_surf);
  
    //! given the loop vertices, compute and fix their points
  MBErrorCode compute_fixed_points(MBEntityHandle dual_surf, MBRange &dverts,
                                   MBRange &face_verts, MBRange &loop_edges);

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

  MBErrorCode label_other_sheets(MBEntityHandle dual_surf,
                                 vtkPolyData *pd,
                                 vtkPolyData *&new_pd);
  
  void label_interior_verts(MBEntityHandle dual_surf,
                            vtkPolyData *pd,
                            vtkRenderer *ren);
  
  MBEntityHandle other_sheet(const MBEntityHandle this_chord,
                             const MBEntityHandle dual_surf);
  
  MBErrorCode get_primal_ids(const MBRange &ents, std::vector<int> &ids);
  
  MBErrorCode allocate_points(MBEntityHandle dual_surf,
                              vtkPolyData *pd,
                              MBRange &verts,
                              MBRange &loop_edges,
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
  
  MBErrorCode construct_graphviz_points(MBEntityHandle dual_surf, 
                                        MBRange &dverts, 
                                        Agsym_t *asym_pos,
                                        GVEntity **dvert_gv);
  
  MBErrorCode construct_graphviz_edges(MBEntityHandle dual_surf, 
                                       MBRange &dedges, 
                                       MBRange &loop_verts, 
                                       Agsym_t *asym_pos, 
                                       GVEntity **dvert_gv, 
                                       GVEntity **dege_gv);
  
  Agsym_t *get_asym(MBEntityHandle dual_surf, const int dim,
                    const char *name, const char *def_val = NULL);
  
  MBErrorCode fixup_degen_bchords(MBEntityHandle dual_surf);

  void print_picked_ent(MBEntityHandle picked_ent);

    //! given some entities, get the corresponding gviz points on the sheet
  void get_points(const MBEntityHandle *ents, const int num_ents, 
                  const bool extra,
                  MBEntityHandle dual_surf, Agnode_t **points);
  
};


#endif
