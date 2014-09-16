#ifdef WIN32
#pragma warning (disable : 4786)
#endif 

#include "NestedRefine.hpp"
#include "moab/HalfFacetRep.hpp"
#include "moab/ReadUtilIface.hpp"
#include <iostream>
#include <assert.h>
#include <vector>


namespace moab{

  NestedRefine::NestedRefine(Core *mesh_in)
  {
    assert(NULL != mesh_in);
    mb = mesh_in;
    HalfFacetRep ahf_mesh(*&mb);
    ahf = ahf_mesh;

    ErrorCode error;
    error = ahf->initialize();
    assert(error==MB_SUCCESS);

    for (int i = 0; i<20; i++)
      {
        access_lmem[i].num_verts = access_lmem[i].num_edges = access_lmem[i].num_faces = access_lmem[i].num_cells = 0;
        access_lmem[i].start_vertex = access_lmem[i].start_edge = access_lmem[i].start_face = access_lmem[i].start_cell = 0;
        access_lmem[i].coordinates = access_lmem[i].edge_conn = access_lmem[i].face_conn = access_lmem[i].cell_conn = NULL;
      }
  }

  ErrorCode NestedRefine::generate_mesh_hierarchy(int *level_seq, int num_level)
  {
    ErrorCode error;
    int hmest[num_level][4];
    error = estimate_hm_storage(level_seq, num_level, hmest);
    if (error != MB_SUCCESS) return error;

    error = generate_hm(level_seq, num_level, hmest);
    if (error != MB_SUCCESS) return error;

  }

  ErrorCode NestedRefine::estimate_hm_storage(int *level_seq, int num_level, int hmest[][4])
  {
    int nverts_in = ahf->_verts.size();
    int nedges_in = ahf->_edges.size();
    int nfaces_in = ahf->_faces.size();
    int ncells_in = ahf->_cells.size();

    hmest[0][0] = nverts_in + new_verts;
    hmest[0][1] = nedges_in*refTemplates[index][level_seq[0]].total_new_ents;
    hmest[0][2] = nfaces_in*refTemplates[index][level_seq[0]].total_new_ents;
    hmest[0][3] = ncells_in*refTemplates[index][level_seq[0]].total_new_ents;

    int nverts, nedges, nfaces, ncells;
    for (int l=1; l<num_level; l++)
      {
        nverts = hmest[l-1][0];
        nedges = hmest[l-1][1];
        nfaces = hmest[l-1][2];
        ncells = hmest[l-1][3];

        hmest[0][0] = nverts_in + new_verts;
        hmest[0][1] = nedges*refTemplates[index][level_seq[0]].total_new_ents;
        hmest[0][2] = nfaces*refTemplates[index][level_seq[0]].total_new_ents;
        hmest[0][3] = ncells*refTemplates[index][level_seq[0]].total_new_ents;

      }
  }

  ErrorCode NestedRefine::create_hm_storage(EntityHandle *set, int cur_level, int estL[4])
  {
    //Obtain chunks of memory for the current level. Add them to a particular meshset.
    //Set flags here
    ErrorCode error = mb->create_meshset(MESHSET_SET,*set);

    ReadUtilIface *read_iface;
    error = mb->query_interface(read_iface);

    //Vertices
    //std::vector<double *> coords;
   // EntityHandle start_nverts;
    error = read_iface->get_node_coords(3, estL[0], 0,  access_lmem[cur_level].start_vertex ,  access_lmem[cur_level].coordinates);
   // access_lmem[cur_level].start_vertex = start_nverts;
    access_lmem[cur_level].num_verts = estL[0];
   // access_lmem[cur_level].coordinates = coords;
    error = mb->add_entities(*set,*start_nverts, estL[0]);

    // Edges
    if (estL[1] != 0)
      {
       // EntityHandle start_nedges, *econnect;
        error = read_iface->get_element_connect(estL[1], 2, MBEDGE, 0, access_lmem[cur_level].start_edge, access_lmem[cur_level].edge_conn);
       // access_lmem[cur_level].start_edge = start_nedges;
        access_lmem[cur_level].num_edges = estL[1];
        //access_lmem[cur_level].edge_conn = econnect;
        error = mb->add_entities(*set,*start_nedges,estL[1]);
      }

    //Faces
    if (estL[2] != 0)
      {
        //EntityHandle start_nfaces, *fconnect;
        EntityHandle first_face = ahf->_faces.begin();
        EntityType type = mb->type_from_handle(first_face);
        int nvpf = ahf->local_maps_2d(first_face);
        error = read_iface->get_element_connect(estL[2], nvpf, type, 0, access_lmem[cur_level].start_face, access_lmem[cur_level].face_conn);
       // access_lmem[cur_level].start_face = start_nfaces;
        access_lmem[cur_level].num_faces = estL[2];
       // access_lmem[cur_level].face_conn = fconnect;
        error = mb->add_entities(*set,*start_nfaces, estL[2]);
      }

    //Cells
    if (estL[3] != 0)
      {
        //EntityHandle start_ncells, *cconnect;
        EntityHandle first_cell = ahf->_cells.begin();
        EntityType type = mb->type_from_handle(first_cell);
        int nvpc = ahf->local_maps_2d(first_cell);
        error = read_iface->get_element_connect(estL[3], nvpc, type, 0, access_lmem[cur_level].start_cell, access_lmem[cur_level].cell_conn);
       // access_lmem[cur_level].start_cell = start_ncells;
        access_lmem[cur_level].num_cells = estL[3];
       // access_lmem[cur_level].cell_conn = cconnect;
        error = mb->add_entities(*set,*start_ncells, estL[3]);
      }
  }

  ErrorCode NestedRefine::generate_hm(int *level_seq, int num_level)
  {
    ErrorCode error;

    EntityHandle set_label[num_level] = {0};

    for (int l = 0; l<num_level; l++)
      {
        //Create arrays for storing the current level
        error = create_hm_storage(&set_label[l], l, hmest[l][4]);

        //Copy the old vertices and add their map
        error = copy_old_vertices(l);

        //Create the new entities and new vertices
        error = construct_hm_entities(l, level_seq[l]);

        //Create ahf maps on the current level
        //error = construct_hm_ahf();
      }
  }

  ErrorCode NestedRefine::copy_old_vertices(int cur_level)
  {
    ErrorCode error;

    if (cur_level)
      {
        int nverts_prev = access_lmem[cur_level-1].num_verts;
        for (int i = 0; i < nverts_prev; i++)
          {
            access_lmem[cur_level].coordinates[0][i] = access_lmem[cur_level-1].coordinates[0][i];
            access_lmem[cur_level].coordinates[1][i] = access_lmem[cur_level-1].coordinates[1][i];
            access_lmem[cur_level].coordinates[2][i] = access_lmem[cur_level-1].coordinates[2][i] ;
          }
      }
    else // Copy the vertices from the input mesh
      {
        int nverts_in = ahf->_verts.size();
        double* x_coords, y_coords, z_coords;

        error = mb->get_coords(ahf->_verts, x_coords, y_coords, z_coords);

        for (int i = 0; i < nverts_in; i++)
          {
            access_lmem[cur_level].coordinates[0][i] = x_coords[i];
            access_lmem[cur_level].coordinates[1][i] = y_coords[i];
            access_lmem[cur_level].coordinates[2][i] = z_coords[i];
          }
      }

    //To add: Map from old vertices to new duplicates.
  }

  ErrorCode NestedRefine::construct_hm_entities(int cur_level, int deg)
  {
    ErrorCode error;

    //Step 1: Create a buffer for vertices. This buffer should be reused during mesh refinement.
    int nvpf, nvpc, vtotal;
   if (ahf->thismeshtype == CURVE)
     {
       int num_new_verts = refTemplates[index][cur_level].total_new_verts;
       vtotal = 2+num_new_verts;
     }
   else if(ahf->thismeshtype == SURFACE || ahf->thismeshtype == SURFACE_MIXED)
     {
      int num_new_verts = refTemplates[index][cur_level].total_new_verts;
      nvpf = ahf->local_maps_2d(first_face);
      vtotal = nvpf+num_new_verts;
     }
   else
     {
       int num_new_verts = refTemplates[index][cur_level].total_new_verts;
       nvpc = ahf->lConnMap3D[index].num_verts_in_cell;
       vtotal = nvpc+num_new_verts;
     }
   EntityHandle *vbuffer = new EntityHandle(vtotal);

    //Step 2: Refine the previous mesh
   if (nvpc) //Refine 3D:VOLUME mesh
     error = construct_hm_3D(cur_level, vbuffer);

   else if (nvpf)  // Refine 2D:SURFACE mesh
     error = construct_hm_2D(cur_level, vbuffer);

   else //Refine 1D:CURVE mesh
     error = construct_hm_1D(cur_level, vbuffer);

  }

  ErrorCode NestedRefine::construct_hm_1D(int cur_level, int deg, EntityHandle *vbuffer)
  {
    ErrorCode error;
    int nverts_prev;
    if (cur_level)
      nverts_prev = access_lmem[cur_level-1].num_verts;
    else
      nverts_prev = ahf->_verts.size();

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int eid = 0; eid < access_lmem[cur_level-1].num_edges; eid++)
      {
        //Add the vertex handles to vbuffer for the current level for the working edge
        EntityHandle edge = access_lmem[cur_level-1].start_edge + eid;
        error = get_connectivity(cur_level-1, edge, conn);
        for (int i=0; i<conn.size(); i++)
          vbuffer[i] = access_lmem[cur_level].start_vertex+ (conn(i)-access_lmem[cur_level-1].start_vertex);
        for (int i=0; i<num_new_verts; i++)
          vbuffer[i+conn.size()] = access_lmem[cur_level].start_vertex+nverts_prev+1;

        //Use the template to obtain the subentities
        for (int i = 0; i < refPatterns[0][deg].total_new_ents; i++)
          {
            int id1 = refPatterns[0][deg].ents_conn[i][0];
            int id2 = refPatterns[0][deg].ents_conn[i][1];
            access_lmem[cur_level].edge_conn[2i] = vbuffer[id1];
            access_lmem[cur_level].edge_conn[2i+1] = vbuffer[id2];
          }
      }

    //Step 2: Obtain the coordinates of the new vertices introduced in the current level

    //Step 3: Construct the ahf maps for the new level

  }

  ErrorCode NestedRefine::construct_hm_2D(int cur_level, int deg, EntityHandle *vbuffer)
  {
    ErrorCode error;
    int nverts_prev, nents_prev;
    if (cur_level)
      {
        nverts_prev = access_lmem[cur_level-1].num_verts;
        nents_prev = access_lmem[cur_level-1].num_faces;
      }
    else
      {
      nverts_prev = ahf->_verts.size();
      nents_prev = ahf->_faces.size();
      }

    //Create some book-keeping arrays over the old mesh to avoid introducing duplicate vertices
    int nepf = ahf->local_maps_2d(first_face);
    EntityHandle *trackverts = new EntityHandle[nepf*nents_prev];
    EntityHandle *flag_ents = new EntityHandle[nents_prev];

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int fid = 0; fid < nents_prev; fid++)
      {
        EntityHandle face = access_lmem[cur_level-1].start_face + fid;

        //Obtain incident edges on this face
        std::vector<EntityHandle> eids;
        error = ahf->get_adjacencies(face,1,eids);

        error = get_connectivity(cur_level-1, face, conn);
        for (int i=0; i<conn.size(); i++)
          vbuffer[i] = access_lmem[cur_level].start_vertex+(conn(i)-access_lmem[cur_level-1].start_vertex);

        //If no incident edge exists, obtain the sibling faces
        if (!eids.empty())
          {

          }
        else
          {
            error = ahf->get_adjacencies()
          };


        //Add the vertex handles to vbuffer for the current level for the working edge

        error = get_connectivity(cur_level-1, edge, conn);
        for (int i=0; i<conn.size(); i++)
          vbuffer[i] = access_lmem[cur_level].start_vertex+ (conn(i)-access_lmem[cur_level-1].start_vertex);
        for (int i=0; i<num_new_verts; i++)
          vbuffer[i+conn.size()] = access_lmem[cur_level].start_vertex+nverts_prev+1;

        //Use the template to obtain the subentities
        for (int i = 0; i < refPatterns[0][deg].total_new_ents; i++)
          {
            int id1 = refPatterns[0][deg].ents_conn[i][0];
            int id2 = refPatterns[0][deg].ents_conn[i][1];
            access_lmem[cur_level].edge_conn[2i] = vbuffer[id1];
            access_lmem[cur_level].edge_conn[2i+1] = vbuffer[id2];
          }
      }

    //Step 2: Obtain the coordinates of the new vertices introduced in the current level

    //Step 3: Construct the ahf maps for the new level



  }

  ErrorCode NestedRefine::construct_hm_3D(int cur_level, int deg, EntityHandle *vbuffer)
  {

  }

  ErrorCode NestedRefine::construct_hm_ahf()
  {

  }

}//namesapce moab
  
