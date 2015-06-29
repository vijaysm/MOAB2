
#include "moab/NestedRefine.hpp"
#include "moab/Templates.hpp"
#include "moab/HalfFacetRep.hpp"
#include "moab/MeasureTime.hpp"
#include "moab/ReadUtilIface.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "moab/ParallelMergeMesh.hpp"
#include "moab/Skinner.hpp"
#endif
#include "Internals.hpp"
#include "MBTagConventions.hpp"

#include <iostream>
#include <assert.h>
#include <vector>
#include <limits>
#include <cmath>


namespace moab{

  NestedRefine::NestedRefine(Core *impl, ParallelComm *comm, EntityHandle rset)
    : mbImpl(impl), pcomm(comm), _rset(rset)
  {
    ErrorCode error;
    assert(NULL != impl);

#ifdef MOAB_HAVE_MPI
    // Get the Parallel Comm instance to prepare all new sets to work in parallel
    // in case the user did not provide any arguments
    if (!comm)
      pcomm = moab::ParallelComm::get_pcomm(mbImpl, 0);
#endif
    error = initialize();
    if (error != MB_SUCCESS)
    {
      std::cout<<"Error initializing NestedRefine\n"<<std::endl;
      exit(1);
    }
  }

  NestedRefine::~NestedRefine()
  {
#ifdef MOAB_HAVE_AHF
    ahf = NULL;
#else
    delete ahf;
#endif

  }

  ErrorCode NestedRefine::initialize()
  {
    ErrorCode error;

    tm = new MeasureTime();
    if (!tm)
      return MB_MEMORY_ALLOCATION_FAILED;

#ifdef MOAB_HAVE_AHF
    ahf = mbImpl->a_half_facet_rep();
#else
    ahf = new HalfFacetRep(mbImpl, pcomm, _rset);
    if (!ahf)
      return MB_MEMORY_ALLOCATION_FAILED;
#endif

    //Check for mixed entity type
    bool chk_mixed = ahf->check_mixed_entity_type();
    if (chk_mixed) MB_SET_ERR(  MB_NOT_IMPLEMENTED, "Encountered a mesh with mixed entity types");

    error = ahf->initialize(); MB_CHK_ERR(error);
    error = ahf->get_entity_ranges(_inverts, _inedges, _infaces, _incells);  MB_CHK_ERR(error);

    //Check for supported entity type
    if (!_incells.empty())
      {
        EntityType type = mbImpl->type_from_handle(_incells[0]);
        if(type != MBTET && type != MBHEX)
          MB_SET_ERR(MB_FAILURE, "Not supported 3D entity types: MBPRISM, MBPYRAMID, MBKNIFE, MBPOLYHEDRON");

        meshdim = 3;
        elementype = type;
      }
   else  if (!_infaces.empty())
      {
        EntityType type = mbImpl->type_from_handle(_infaces[0]);
        if(type == MBPOLYGON)
          MB_SET_ERR(MB_FAILURE, "Not supported 2D entity type: POLYGON");

        meshdim = 2;
        elementype = type;
      }
    else if (!_inedges.empty())
      {
        meshdim = 1;
        elementype = MBEDGE;
      }
    else MB_SET_ERR(MB_NOT_IMPLEMENTED, "Encountered a mixed-dimensional or invalid mesh");

    //Initialize std::map to get indices of degrees.
    deg_index[2] = 0;
    deg_index[3] = 1;
    deg_index[5] = 2;

    //Set ghost flag to false
      hasghost = false;

    return MB_SUCCESS;
  }

  /************************************************************
   *     Interface Functions                                  *
   ************************************************************/

  ErrorCode NestedRefine::generate_mesh_hierarchy(int num_level, int *level_degrees, std::vector<EntityHandle> &level_sets)
  { 
    assert(num_level > 0);

    ErrorCode error;
    moab::EntityHandle *hmsets = new moab::EntityHandle[num_level];

    if (meshdim <=2)
      {
        for (int i=0; i<num_level; i++)
          {
            assert((level_degrees[i]==2) || (level_degrees[i]==3)  || (level_degrees[i]==5 ));
            level_dsequence[i] = level_degrees[i];
          }
      }
    else
      {
        for (int i=0; i<num_level; i++)
          {
            assert((level_degrees[i]==2) || (level_degrees[i]==3) );
            level_dsequence[i] = level_degrees[i];
          }
      }

    error = generate_hm(level_degrees, num_level, hmsets); MB_CHK_ERR(error);

    // copy the entity handles
    level_sets.resize(num_level + 1);
    level_sets[0] = _rset;
    for (int i=0; i<num_level; i++)
      level_sets[i+1]=hmsets[i];

    delete [] hmsets;
    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::get_connectivity(EntityHandle ent, int level, std::vector<EntityHandle> &conn)
  {
    ErrorCode error;
    EntityType type = mbImpl->type_from_handle(ent);
    EntityHandle start_ent ;
    if (!conn.empty())
      conn.clear();
    if (level> 0)
      {
        if (type == MBEDGE)
          {
            conn.reserve(2);
            start_ent = level_mesh[level-1].start_edge;
            EntityID offset = ID_FROM_HANDLE(ent) - ID_FROM_HANDLE(start_ent);
            conn.push_back( level_mesh[level - 1].edge_conn[2 * offset]);
            conn.push_back( level_mesh[level - 1].edge_conn[2 * offset + 1]);

          }
        else if (type == MBTRI || type == MBQUAD)
          {
            int num_corners = ahf->lConnMap2D[type-2].num_verts_in_face;
            conn.reserve(num_corners);
            start_ent = level_mesh[level-1].start_face;
            EntityID offset = ID_FROM_HANDLE(ent) - ID_FROM_HANDLE(start_ent);

            for (int i=0; i<num_corners; i++)
              conn.push_back(level_mesh[level-1].face_conn[num_corners* offset+i]);
          }
        else if (type == MBTET || type == MBHEX)
          {
            int index = ahf->get_index_in_lmap(*_incells.begin());
            int num_corners = ahf->lConnMap3D[index].num_verts_in_cell;
            conn.reserve(num_corners);
            start_ent = level_mesh[level-1].start_cell;
            EntityID offset = ID_FROM_HANDLE(ent) - ID_FROM_HANDLE(start_ent);
            for (int i=0; i<num_corners; i++)
              conn.push_back(level_mesh[level-1].cell_conn[num_corners*offset+i]);
          }
        else
          MB_SET_ERR(MB_FAILURE, "Requesting connectivity for an unsupported entity type");

      }
      else
      {
        error = mbImpl->get_connectivity(&ent, 1, conn); MB_CHK_ERR(error);
      }

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::get_coordinates(EntityHandle *verts, int num_verts, int level, double *coords)
  {
    if (level >0){
        EntityID vstart = ID_FROM_HANDLE(level_mesh[level - 1].start_vertex);
        for (int i = 0; i < num_verts; i++)
          {
            const EntityHandle &vid = verts[i];
            EntityID offset = ID_FROM_HANDLE(vid) - vstart;
            coords[3 * i]   = level_mesh[level - 1].coordinates[0][offset];
            coords[3 * i + 1] = level_mesh[level - 1].coordinates[1][offset];
            coords[3 * i + 2] = level_mesh[level - 1].coordinates[2][offset];
          }
      }
    else
      {
        ErrorCode error;
        error = mbImpl->get_coords(verts, num_verts, coords); MB_CHK_ERR(error);
      }

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::get_adjacencies(const EntityHandle source_entity,
                            const unsigned int target_dimension,
                            std::vector<EntityHandle> &target_entities)

  {
    ErrorCode error;
    error = ahf->get_adjacencies(source_entity, target_dimension, target_entities); MB_CHK_ERR(error);

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::child_to_parent(EntityHandle child, int child_level, int parent_level, EntityHandle *parent)
  {
    assert((child_level>0) &&(child_level>parent_level));
    EntityType type = mbImpl->type_from_handle(child);
    assert(type != MBVERTEX);

    int child_index;
    if (type == MBEDGE)
      child_index = child - level_mesh[child_level-1].start_edge;
    else if (type == MBTRI || type == MBQUAD)
      child_index = child - level_mesh[child_level-1].start_face;
    else if (type == MBTET || type == MBHEX)
      child_index = child - level_mesh[child_level-1].start_cell;
    else
      MB_SET_ERR(MB_FAILURE, "Requesting parent for unsupported entity type");


    int parent_index;
    int l = child_level - parent_level;
    for (int i=0; i< l; i++)
      {
        int d = get_index_from_degree(level_dsequence[child_level-i-1]);
        int nch = refTemplates[type-1][d].total_new_ents;
        child_index = child_index/nch;
      }
     parent_index = child_index;

    if (type == MBEDGE)
      {
        if (parent_level>0)
          *parent = level_mesh[parent_level-1].start_edge + parent_index;
        else
          *parent = _inedges[parent_index];
      }
    else if (type == MBTRI || type == MBQUAD)
      {
        if (parent_level>0)
          *parent = level_mesh[parent_level-1].start_face + parent_index;
        else
          *parent = _infaces[parent_index];
      }
    else if (type == MBTET || type == MBHEX)
      {
        if (parent_level>0)
          *parent = level_mesh[parent_level-1].start_cell + parent_index;
        else
          *parent = _incells[parent_index];
      }

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::parent_to_child(EntityHandle parent, int parent_level, int child_level,  std::vector<EntityHandle> &children)
  {
    assert ((child_level>0) && (child_level>parent_level));
    EntityType type = mbImpl->type_from_handle(parent);
    assert(type != MBVERTEX);

    int parent_index;
    if (type == MBEDGE)
      {
        if(parent_level>0)
          parent_index = parent - level_mesh[parent_level-1].start_edge;
        else
          parent_index = _inedges.index(parent);
      }
    else if (type == MBTRI || type == MBQUAD)
      {
        if (parent_level>0)
          parent_index = parent - level_mesh[parent_level-1].start_face;
        else
          parent_index = _infaces.index(parent);
      }
    else if (type == MBTET || type == MBHEX)
      {
        if (parent_level>0)
          parent_index = parent - level_mesh[parent_level-1].start_cell;
        else
          parent_index = _incells.index(parent);
      }
    else
      MB_SET_ERR(MB_FAILURE, "Requesting children for unsupported entity type");


    int start, end;
    start = end = parent_index;
    for (int i=parent_level; i< child_level; i++)
      {
        int d = get_index_from_degree(level_dsequence[i]);
        int nch = refTemplates[type-1][d].total_new_ents;
        start = start*nch;
        end = end*nch + nch-1;
      }

    int num_child = end-start;
    children.reserve(num_child);

    for (int i=start; i<=end; i++)
      {
        EntityHandle child;
        if (type == MBEDGE)
          child = level_mesh[child_level-1].start_edge + i;
        else if (type == MBTRI || type == MBQUAD)
          child = level_mesh[child_level-1].start_face + i;
        else if (type == MBTET || type == MBHEX)
          child = level_mesh[child_level-1].start_cell + i;

        children.push_back(child);
      }

    return MB_SUCCESS;

  }

  ErrorCode NestedRefine::vertex_to_entities(EntityHandle vertex, int level, std::vector<EntityHandle> &incident_entities)
  {
    assert(level>=0);
    ErrorCode error;

    //Step 1: Get the incident entities at the current level
    std::vector<EntityHandle> inents;
    if (meshdim == 1)
      {
        error = ahf->get_up_adjacencies_1d(vertex, inents); MB_CHK_ERR(error);
      }
    else if (meshdim == 2)
      {
        error = ahf->get_up_adjacencies_vert_2d(vertex, inents); MB_CHK_ERR(error);
      }
    else if (meshdim == 3)
      {
        error = ahf->get_up_adjacencies_vert_3d(vertex, inents); MB_CHK_ERR(error);
      }

    //Step 2: Loop over all the incident entities at the current level and gather their parents
    for (int i=0; i< (int)inents.size(); i++ )
      {
        EntityHandle ent = inents[i];
        EntityHandle parent;
        error = child_to_parent(ent, level, level-1, &parent); MB_CHK_ERR(error);
        incident_entities.push_back(parent);
      }

    //Step 3: Sort and remove duplicates
    std::sort(incident_entities.begin(), incident_entities.end());
    incident_entities.erase(std::unique(incident_entities.begin(), incident_entities.end()), incident_entities.end());

    return MB_SUCCESS;
  }

  bool NestedRefine::is_entity_on_boundary(const EntityHandle &entity)
  {
    bool is_border = false;
    EntityType type = mbImpl->type_from_handle(entity);

    if (type == MBVERTEX)
      is_border = is_vertex_on_boundary(entity);
    else if (type == MBEDGE )
      is_border = is_edge_on_boundary(entity);
    else if (type == MBTRI || type == MBQUAD )
      is_border = is_face_on_boundary(entity);
    else if (type == MBTET || type == MBHEX)
      is_border = is_cell_on_boundary(entity);
    else
      MB_SET_ERR(MB_FAILURE, "Requesting boundary information for unsupported entity type");

    return is_border;
  }

  ErrorCode NestedRefine::exchange_ghosts(std::vector<EntityHandle> &lsets, int num_glayers)
  {
    ErrorCode error;

    if (hasghost)
      return MB_SUCCESS;

    hasghost = true;
#ifdef MOAB_HAVE_MPI
    error = pcomm->exchange_ghost_cells(meshdim, 0, num_glayers, 0, true, false);MB_CHK_ERR(error);
#else
    MB_SET_ERR(MB_FAILURE,"Requesting ghost layers for a serial mesh");
#endif

    Range  * lverts = new Range[lsets.size()];
    Range *  lents   = new Range[lsets.size()];
    for (size_t i=0; i<lsets.size(); i++)
    {
      error = mbImpl->get_entities_by_dimension(lsets[i], meshdim, lents[i]);MB_CHK_ERR(error);
      error = mbImpl->get_connectivity(lents[i], lverts[i]);MB_CHK_ERR(error);

      for (int gl =0; gl< num_glayers; gl++)
      {
        error = mbImpl->get_adjacencies(lverts[i],meshdim,false,lents[i],Interface::UNION);MB_CHK_ERR(error);
        error = mbImpl->get_connectivity(lents[i], lverts[i]);MB_CHK_ERR(error);
      }
    }
    for (size_t i=0; i<lsets.size(); i++)
    {
      error = mbImpl->add_entities(lsets[i], lverts[i]);MB_CHK_ERR(error);
      error = mbImpl->add_entities(lsets[i], lents[i]);MB_CHK_ERR(error);
    }

    delete [] lverts;
    delete [] lents;
    return MB_SUCCESS;
  }



  /***********************************************
   *  Basic functionalities: generate HM         *
   ***********************************************/

  ErrorCode NestedRefine::estimate_hm_storage(EntityHandle set, int level_degree, int cur_level, int hmest[4])
  {
    ErrorCode error;

    //Obtain the size of input mesh.
    int nverts_prev, nedges_prev, nfaces_prev, ncells_prev;
    if (cur_level)
      {
        nverts_prev = level_mesh[cur_level-1].num_verts;
        nedges_prev = level_mesh[cur_level-1].num_edges;
        nfaces_prev = level_mesh[cur_level-1].num_faces;
        ncells_prev = level_mesh[cur_level-1].num_cells;
      }
      else
      {
        nverts_prev = _inverts.size();
        nedges_prev = _inedges.size();
        nfaces_prev = _infaces.size();
        ncells_prev = _incells.size();
      }

    //Estimate mesh size of current level mesh.
    int nedges=0, nfaces=0;
    error = count_subentities(set, cur_level-1, &nedges, &nfaces); MB_CHK_ERR(error);

    int d = get_index_from_degree(level_degree);
    int nverts = refTemplates[MBEDGE-1][d].nv_edge*nedges;
    hmest[0] = nverts_prev + nverts;
    hmest[1] = nedges_prev*refTemplates[MBEDGE-1][d].total_new_ents;
    hmest[2] = 0;
    hmest[3] = 0;

    int findex, cindex;
    if (nfaces_prev != 0)
      {
        EntityHandle start_face;
        if (cur_level)
          start_face = level_mesh[cur_level - 1].start_face;
        else
          start_face = *_infaces.begin();
        findex = mbImpl->type_from_handle(start_face) - 1;
        hmest[2] = nfaces_prev * refTemplates[findex][d].total_new_ents;

        if (meshdim == 2)
          hmest[0] += refTemplates[findex][d].nv_face * nfaces_prev;

        if (meshdim == 3)
          hmest[1] += nfaces_prev * intFacEdg[findex - 1][d].nie;
      }

      if (ncells_prev != 0)
      {
        cindex = mbImpl->type_from_handle(*(_incells.begin())) - 1;
        hmest[3] = ncells_prev * refTemplates[cindex][d].total_new_ents;

        hmest[0] += refTemplates[cindex][d].nv_face * nfaces;
        hmest[0] += refTemplates[cindex][d].nv_cell * ncells_prev;
      }

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::create_hm_storage_single_level(EntityHandle *set, int cur_level, int estL[4])
  {
    //Obtain chunks of memory for the current level. Add them to a particular meshset.
    EntityHandle set_handle;
    ErrorCode error = mbImpl->create_meshset(MESHSET_SET, set_handle); MB_CHK_SET_ERR(error, "Cannot create mesh for the current level");
    *set = set_handle;

    ReadUtilIface *read_iface;
    error = mbImpl->query_interface(read_iface); MB_CHK_ERR(error);

    //Vertices
    error = read_iface->get_node_coords(3, estL[0], 0,  level_mesh[cur_level].start_vertex , level_mesh[cur_level].coordinates); MB_CHK_ERR(error);
    level_mesh[cur_level].num_verts = estL[0];

    Range newverts(level_mesh[cur_level].start_vertex, level_mesh[cur_level].start_vertex+estL[0] - 1);
    error = mbImpl->add_entities(*set, newverts); MB_CHK_ERR(error);
    level_mesh[cur_level].verts = newverts;

    Tag gidtag;
    error = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, gidtag);MB_CHK_ERR(error);
    error = read_iface->assign_ids(gidtag, newverts, level_mesh[cur_level].start_vertex);MB_CHK_ERR(error);

    // Edges
    if (estL[1])
      {
        error = read_iface->get_element_connect(estL[1], 2, MBEDGE, 0, level_mesh[cur_level].start_edge, level_mesh[cur_level].edge_conn); MB_CHK_ERR(error);
        level_mesh[cur_level].num_edges = estL[1];

        Range newedges(level_mesh[cur_level].start_edge, level_mesh[cur_level].start_edge+estL[1] - 1);
        error = mbImpl->add_entities(*set, newedges); MB_CHK_ERR(error);
        level_mesh[cur_level].edges = newedges;
      }
    else
      level_mesh[cur_level].num_edges = 0;

    //Faces
    if (estL[2])
      {
        EntityType type = mbImpl->type_from_handle(*(_infaces.begin()));
        int nvpf = ahf->lConnMap2D[type-2].num_verts_in_face;
        error = read_iface->get_element_connect(estL[2], nvpf, type, 0, level_mesh[cur_level].start_face, level_mesh[cur_level].face_conn); MB_CHK_ERR(error);
        level_mesh[cur_level].num_faces = estL[2];

        Range newfaces(level_mesh[cur_level].start_face, level_mesh[cur_level].start_face+estL[2] - 1);
        error = mbImpl->add_entities(*set, newfaces); MB_CHK_ERR(error);
        level_mesh[cur_level].faces = newfaces;
      }
    else
      level_mesh[cur_level].num_faces = 0;

    //Cells
    if (estL[3])
      {
        EntityType type = mbImpl->type_from_handle(*(_incells.begin()));
        int index = ahf->get_index_in_lmap(*_incells.begin());
        int nvpc = ahf->lConnMap3D[index].num_verts_in_cell;
        error = read_iface->get_element_connect(estL[3], nvpc, type, 0, level_mesh[cur_level].start_cell, level_mesh[cur_level].cell_conn); MB_CHK_ERR(error);
        level_mesh[cur_level].num_cells = estL[3];

        Range newcells(level_mesh[cur_level].start_cell, level_mesh[cur_level].start_cell+estL[3] - 1);
        error = mbImpl->add_entities(*set, newcells); MB_CHK_ERR(error);
        level_mesh[cur_level].cells = newcells;
      }
    else
      level_mesh[cur_level].num_cells = 0;

    //Resize the ahf maps
    error = ahf->resize_hf_maps(level_mesh[cur_level].start_vertex, level_mesh[cur_level].num_verts, level_mesh[cur_level].start_edge, level_mesh[cur_level].num_edges, level_mesh[cur_level].start_face, level_mesh[cur_level].num_faces, level_mesh[cur_level].start_cell, level_mesh[cur_level].num_cells); MB_CHK_ERR(error);


    error = ahf->update_entity_ranges(*set); MB_CHK_ERR(error);

    //If the mesh type changes, then update the member variable in ahf to use the applicable adjacency matrix
    MESHTYPE nwmesh = ahf->get_mesh_type(level_mesh[cur_level].num_verts, level_mesh[cur_level].num_edges, level_mesh[cur_level].num_faces,level_mesh[cur_level].num_cells);MB_CHK_ERR(error);
    if (ahf->thismeshtype != nwmesh)
      ahf->thismeshtype = nwmesh;

    return MB_SUCCESS;
  }



  /**********************************
   *   Hierarchical Mesh Generation  *
   * *********************************/

  ErrorCode NestedRefine::generate_hm(int *level_degrees, int num_level, EntityHandle *hm_set)
  {
    ErrorCode error;

    Tag gidtag;
    error = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, gidtag);MB_CHK_ERR(error);

    timeall.tm_total = 0;
    timeall.tm_refine = 0;
    timeall.tm_presolve = 0;

    for (int l = 0; l<num_level; l++)
      {
        std::cout<<"Starting level = "<<l<<std::endl;
        double tstart;

#ifdef MOAB_HAVE_MPI
        tstart = MPI_Wtime();
#else
        tstart = tm->wtime();
#endif

        // Estimate storage
        int hmest[4] = {0,0,0,0};
        EntityHandle set;
        if (l)
          set = hm_set[l-1];
        else
          set =  _rset;
        error = estimate_hm_storage(set, level_degrees[l], l, hmest); MB_CHK_ERR(error);

        //Create arrays for storing the current level
        error = create_hm_storage_single_level(&hm_set[l], l, hmest); MB_CHK_ERR(error);

        //Copy the old vertices along with their coordinates
        error = copy_vertices_from_prev_level(l); MB_CHK_ERR(error);

        //Create the new entities and new vertices
        error = construct_hm_entities(l, level_degrees[l]); MB_CHK_ERR(error);

#ifdef MOAB_HAVE_MPI
        timeall.tm_refine += MPI_Wtime() - tstart;
#else
        timeall.tm_refine += tm->wtime() - tstart;
#endif
        //timeall.tm_refine += tm->wtime() - tstart;

        // Go into parallel communication
#ifdef MOAB_HAVE_MPI
        if (pcomm)
        {
          // TEMP: Add the adjacencies for MOAB-native DS
          // NOTE (VSM): This is expensive since it creates a doubly
          // redundant copy of the adjacency data in both MOAB-native
          // and AHF. Need to fix this with AHF optimized branch.
          ReadUtilIface *read_iface;
          error = mbImpl->query_interface(read_iface);MB_CHK_ERR(error);
          if (level_mesh[l].num_edges != 0)
          {
            error = read_iface->update_adjacencies(level_mesh[l].start_edge, level_mesh[l].num_edges, 2, level_mesh[l].edge_conn);MB_CHK_ERR(error);
          }
          if (level_mesh[l].num_faces != 0)
          {
              EntityType type = mbImpl->type_from_handle(*(_infaces.begin()));
              int nvpf = ahf->lConnMap2D[type - 2].num_verts_in_face;
              error = read_iface->update_adjacencies(level_mesh[l].start_face, level_mesh[l].num_faces, nvpf, level_mesh[l].face_conn);MB_CHK_ERR(error);
          }
          if (level_mesh[l].num_cells != 0)
          {
              int index = ahf->get_index_in_lmap(*_incells.begin());
              int nvpc = ahf->lConnMap3D[index].num_verts_in_cell;
              error = read_iface->update_adjacencies(level_mesh[l].start_cell, level_mesh[l].num_cells, nvpc, level_mesh[l].cell_conn);MB_CHK_ERR(error);
          }

          if (pcomm->size() > 1)
            {
             // double tpstart = tm->wtime();
              double tpstart = MPI_Wtime();

              // get all entities on the rootset
              moab::Range vtxs, edgs, facs, elms;
              error = mbImpl->get_entities_by_dimension(hm_set[l], 0, vtxs, false);MB_CHK_ERR(error);
              error = mbImpl->get_entities_by_dimension(hm_set[l], 1, edgs, false);MB_CHK_ERR(error);
              error = mbImpl->get_entities_by_dimension(hm_set[l], 2, facs, false);MB_CHK_ERR(error);
              error = mbImpl->get_entities_by_dimension(hm_set[l], 3, elms, false);MB_CHK_ERR(error);

              // set the parallel partition tag data
              moab::Tag part_tag;
              moab::EntityHandle part_set;
              int partid = pcomm->rank(), dum_id = -1;
              error = mbImpl->tag_get_handle("PARALLEL_PARTITION", 1, moab::MB_TYPE_INTEGER,
                                             part_tag, moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE, &dum_id);MB_CHK_ERR(error);

              error = mbImpl->create_meshset(moab::MESHSET_SET, part_set);MB_CHK_ERR(error);
              error = mbImpl->add_entities(part_set, vtxs);MB_CHK_ERR(error);
              error = mbImpl->add_entities(part_set, edgs);MB_CHK_ERR(error);
              error = mbImpl->add_entities(part_set, facs);MB_CHK_ERR(error);
              error = mbImpl->add_entities(part_set, elms);MB_CHK_ERR(error);
              error = mbImpl->add_entities(hm_set[l],&part_set,1);MB_CHK_ERR(error);
              error = mbImpl->tag_set_data(part_tag, &part_set, 1, &partid);MB_CHK_ERR(error);

              error = mbImpl->tag_set_data(part_tag, &hm_set[l], 1, &partid);MB_CHK_ERR(error);
              //mbImpl->list_entities(vtxs);
              //mbImpl->list_entity(hm_set[l]);

              //
              // Now that we have the local piece of the mesh refined consistently,
              // call parallel merge instead of resolved_shared to stitch together the meshes
              // and to handle parallel communication of remote proc/entity-handle pairs for
              // shared + non-owned interfaces entities.
              //
              // TODO: This needs to be replaced by the following scheme in the future as
              // an optimization step.
              //   > Assign global IDs consistently for shared entities so that parallel
              //   > resolve shared ents can happen out of the box. This is the fastest option.

              ParallelMergeMesh pm(pcomm, 1e-08);
              error = pm.merge(hm_set[l], true);MB_CHK_ERR(error);

              // std::cout<<"Writing level set"<<std::endl;
              // mbImpl->write_file("test.h5m", 0, ";;PARALLEL=WRITE_PART;DEBUG_IO=3", &hm_set[l], 1);

              timeall.tm_presolve += MPI_Wtime() - tpstart;
              //
              // Parallel Communication complete - all entities resolved
              //

              {
                // Assign new global IDs for all the entities we just generated to maintain contiguity
                // Range pents[4] = {vtxs, edgs, facs, elms};
                // error = pcomm->assign_global_ids(pents, 3, 1, true, true);MB_CHK_ERR(error);
                // error = pcomm->assign_global_ids(hm_set[l], 0, 1, false, true, false);MB_CHK_ERR(error);

                // Now that we have resolved all shared entities in parallel,
                // exchange GLOBAL_ID data for all local entities so that we are
                // up to date on remote changes.
                error = pcomm->exchange_tags(gidtag,vtxs);MB_CHK_ERR(error);
                error = pcomm->exchange_tags(gidtag,edgs);MB_CHK_ERR(error);
                error = pcomm->exchange_tags(gidtag,facs);MB_CHK_ERR(error);
                error = pcomm->exchange_tags(gidtag,elms);MB_CHK_ERR(error);
              }
             // timeall.tm_presolve += tm->wtime() - tpstart;
            }
          }
#endif
      }
    timeall.tm_total = timeall.tm_refine + timeall.tm_presolve;
  //  mbImpl->write_file("test.h5m", 0, ";;PARALLEL=WRITE_PART");

    return MB_SUCCESS;
  }


  ErrorCode NestedRefine::construct_hm_entities(int cur_level, int deg)
  {
    ErrorCode error;

    //Generate mesh for current level by refining previous level.
    if (ahf->thismeshtype == CURVE)
      {
        error = construct_hm_1D(cur_level, deg); MB_CHK_ERR(error);
      }
     else if(ahf->thismeshtype == SURFACE || ahf->thismeshtype == SURFACE_MIXED)
      {
        error = construct_hm_2D(cur_level, deg); MB_CHK_ERR(error);
      }
    else
      {
        error = construct_hm_3D(cur_level, deg); MB_CHK_ERR(error);
      }

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::construct_hm_1D(int cur_level, int deg)
  {
    ErrorCode error;
    int nverts_prev, nents_prev;
    if (cur_level)
      {
        nverts_prev = level_mesh[cur_level-1].num_verts;
        nents_prev = level_mesh[cur_level-1].num_edges;
      }
    else
      {
        nverts_prev = _inverts.size();
        nents_prev = _inedges.size();
      }

    int d = get_index_from_degree(deg);
    int vtotal = 2+ refTemplates[0][d].total_new_verts;
    EntityHandle *vbuffer = new EntityHandle[vtotal];

    std::vector<EntityHandle> conn;
    int count_nents = 0;
    int count_verts = nverts_prev;

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int eid = 0; eid < nents_prev; eid++)
      {
        conn.clear();

        // EntityHandle of the working edge
        EntityHandle edge;
        if (cur_level)
          edge = level_mesh[cur_level-1].start_edge + eid;
        else
          edge = _inedges[eid]; // Makes the assumption initial mesh is contiguous in memory

        error = get_connectivity(edge, cur_level, conn); MB_CHK_ERR(error);

        //Add the vertex handles to vbuffer for the current level for the working edge

        // Since the old vertices are copied first, their local indices do not change as new levels are added.
        // Clearly the local indices of the new vertices introduced in the current level is still the same
        // when the old vertices are copied. Thus, there is no need to explicitly store another map between
        // the old and duplicates in the subsequent levels. The second part in the following sum is the local
        // index in the previous level.

        // Add the corners to the vbuffer first.

        for (int i=0; i<(int)conn.size(); i++)
          {
            if (cur_level)
              vbuffer[i] = level_mesh[cur_level].start_vertex+ (conn[i] - level_mesh[cur_level-1].start_vertex);
            else
              vbuffer[i] = level_mesh[cur_level].start_vertex + (conn[i] - *_inverts.begin());
          }

        //Adding rest of the entityhandles to working buffer for vertices.
        int num_new_verts = refTemplates[0][d].total_new_verts;

        for (int i=0; i<num_new_verts; i++)
          {
            vbuffer[i+2] = level_mesh[cur_level].start_vertex+count_verts;
            count_verts += 1;
          }

        //Use the template to obtain the subentities
        int id1, id2;
        int etotal = refTemplates[0][d].total_new_ents;
        EntityHandle *ent_buffer = new EntityHandle[etotal];

        for (int i = 0; i < etotal; i++)
          {
            id1 = refTemplates[0][d].ents_conn[i][0];
            id2 = refTemplates[0][d].ents_conn[i][1];
            level_mesh[cur_level].edge_conn[2*(count_nents)] = vbuffer[id1];
            level_mesh[cur_level].edge_conn[2*(count_nents)+1] = vbuffer[id2];
            ent_buffer[i] = level_mesh[cur_level].start_edge+count_nents;
            count_nents += 1;
          };

        error = update_local_ahf(deg, MBEDGE,  vbuffer, ent_buffer, etotal); MB_CHK_ERR(error);

        // Compute the coordinates of the new vertices: Linear interpolation
        int idx;
        double xi;
        for (int i=0; i< num_new_verts; i++)
          {
            xi = refTemplates[0][d].vert_nat_coord[i][0];
            idx =  vbuffer[i+2] - level_mesh[cur_level].start_vertex; // index of new vertex in current level
            if (cur_level){
                id1 = conn[0]-level_mesh[cur_level-1].start_vertex; //index of old end vertices in current level
                id2 = conn[1]-level_mesh[cur_level-1].start_vertex;
              }
            else
              {
                id1 = _inverts.index(conn[0]);
                id2 = _inverts.index(conn[1]);
              }

            level_mesh[cur_level].coordinates[0][idx] = (1-xi)*level_mesh[cur_level].coordinates[0][id1] + xi*level_mesh[cur_level].coordinates[0][id2];
            level_mesh[cur_level].coordinates[1][idx] = (1-xi)*level_mesh[cur_level].coordinates[1][id1] + xi*level_mesh[cur_level].coordinates[1][id2];
            level_mesh[cur_level].coordinates[2][idx] = (1-xi)*level_mesh[cur_level].coordinates[2][id1] + xi*level_mesh[cur_level].coordinates[2][id2];
          }

        delete [] ent_buffer;
      }

    error = update_global_ahf(MBEDGE, cur_level, deg); MB_CHK_ERR(error);

    delete [] vbuffer;

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::construct_hm_1D(int cur_level, int deg, EntityType type, std::vector<EntityHandle> &trackverts)
  {
    ErrorCode error;

    int  nedges_prev;
    if (cur_level)
      nedges_prev = level_mesh[cur_level-1].num_edges;
    else
      nedges_prev = _inedges.size();

    int d = get_index_from_degree(deg);
    int nve = refTemplates[0][d].nv_edge;
    int vtotal =  2+ refTemplates[0][d].total_new_verts;
    int etotal = refTemplates[0][d].total_new_ents;
    int ne=0, dim=0, index=0;
    if (type == MBTRI || type == MBQUAD)
      {
        index = type-2;
        ne = ahf->lConnMap2D[index].num_verts_in_face;
        dim = 2;
      }
    else if (type == MBTET || type == MBHEX)
      {
        index = ahf->get_index_in_lmap(*(_incells.begin()));
        ne =  ahf->lConnMap3D[index].num_edges_in_cell;
        dim = 3;
      }

    EntityHandle *vbuffer = new EntityHandle[vtotal];
    EntityHandle *ent_buffer = new EntityHandle[etotal];

    std::vector<EntityHandle> adjents, econn, fconn;
    std::vector<int> leids;
    int count_nents = 0;

    //Loop over all the edges and gather the vertices to be used for refinement
    for (int eid=0; eid< nedges_prev; eid++)
      {
        adjents.clear();
        leids.clear();
        econn.clear();
        fconn.clear();
        for (int i=0; i<vtotal; i++)
          vbuffer[i] = 0;
        for (int i=0; i<etotal; i++)
          ent_buffer[i] = 0;

        EntityHandle edge;
        if (cur_level)
          edge = level_mesh[cur_level-1].start_edge + eid;
        else
          edge = _inedges[eid];

        error = get_connectivity(edge, cur_level, econn); MB_CHK_ERR(error);

        for (int i=0; i<(int)econn.size(); i++)
          {
            if (cur_level)
              vbuffer[i] = level_mesh[cur_level].start_vertex+ (econn[i] - level_mesh[cur_level-1].start_vertex);
            else
              vbuffer[i] = level_mesh[cur_level].start_vertex + (econn[i] - *_inverts.begin());
          }

        int fid=-1, lid=-1, idx1=-1, idx2=-1;

        if (dim==2)
          {
            error = ahf->get_up_adjacencies_2d(edge, adjents, &leids); MB_CHK_ERR(error);
            if (cur_level)
              fid = adjents[0] - level_mesh[cur_level-1].start_face;
            else
              fid = _infaces.index(adjents[0]);

            lid = leids[0];
            idx1 = lid;
            idx2 = ahf->lConnMap2D[index].next[lid];
          }
        else if (dim==3)
          {
            error = ahf->get_up_adjacencies_edg_3d(edge, adjents, &leids); MB_CHK_ERR(error);
            if (cur_level)
              fid = adjents[0] - level_mesh[cur_level-1].start_cell;
            else
              fid = _incells.index(adjents[0]);

            lid = leids[0];
            idx1 = ahf->lConnMap3D[index].e2v[lid][0];
            idx2 = ahf->lConnMap3D[index].e2v[lid][1];
          }

        error = get_connectivity(adjents[0], cur_level, fconn); MB_CHK_ERR(error);

        bool orient = false;
        if ((fconn[idx1] == econn[0])&&(fconn[idx2] == econn[1]))
          orient = true;

           if (orient)
             {
               for (int j=0; j<nve; j++)
                 vbuffer[j+2] = trackverts[fid*ne*nve+nve*lid+j];
             }
           else
             {
               for (int j=0; j<nve; j++)
                 vbuffer[(nve-j-1)+2] = trackverts[fid*ne*nve+nve*lid+j];
             }

           //Use the template to obtain the subentities
           int id1, id2;

           for (int i = 0; i < etotal; i++)
             {
               id1 = refTemplates[0][d].ents_conn[i][0];
               id2 = refTemplates[0][d].ents_conn[i][1];
               level_mesh[cur_level].edge_conn[2*(count_nents)] = vbuffer[id1];
               level_mesh[cur_level].edge_conn[2*(count_nents)+1] = vbuffer[id2];
               ent_buffer[i] = level_mesh[cur_level].start_edge+count_nents;

               count_nents += 1;
             };

           error = update_local_ahf(deg, MBEDGE,  vbuffer, ent_buffer, etotal); MB_CHK_ERR(error);
      }

    error = update_global_ahf_1D_sub(cur_level, deg); MB_CHK_ERR(error);

    delete [] vbuffer;
    delete [] ent_buffer;

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::construct_hm_2D(int cur_level, int deg)
  {
    ErrorCode error;
    int nverts_prev, nents_prev;
    if (cur_level)
      {
        nverts_prev = level_mesh[cur_level-1].num_verts;
        nents_prev = level_mesh[cur_level-1].num_faces;
      }
    else
      {
        nverts_prev = _inverts.size();
        nents_prev = _infaces.size();
      }

    //Create some book-keeping arrays over the old mesh to avoid introducing duplicate vertices and calculating vertices more than once.
    EntityType ftype = mbImpl->type_from_handle(*_infaces.begin());
    int nepf = ahf->lConnMap2D[ftype-2].num_verts_in_face;
    int findex = ftype-1;

    int d = get_index_from_degree(deg);
    int tnv = refTemplates[findex][d].total_new_verts;
    int vtotal = nepf + tnv;
    EntityHandle *vbuffer = new EntityHandle[vtotal];
    int etotal = refTemplates[findex][d].total_new_ents;
    EntityHandle *ent_buffer = new EntityHandle[etotal];

    int nve = refTemplates[findex][d].nv_edge;
    std::vector<EntityHandle> trackvertsF(nents_prev*nepf*nve, 0);
    int cur_nverts = level_mesh[cur_level].num_verts;
    std::vector<int> flag_verts(cur_nverts-nverts_prev, 0);

    int count_nverts = nverts_prev;
    int count_nents = 0;
    std::vector<EntityHandle> conn, cur_conn;

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int fid = 0; fid < nents_prev; fid++)
      {        
        conn.clear();
        cur_conn.clear();
        for (int i=0; i<vtotal; i++)
          vbuffer[i] = 0;
        for (int i=0; i<etotal; i++)
          ent_buffer[i] = 0;

        //EntityHandle of the working face
        EntityHandle face;
        if (cur_level)
          face = level_mesh[cur_level-1].start_face + fid;
        else
          face = _infaces[fid];

        error = get_connectivity(face, cur_level, conn); MB_CHK_ERR(error);

        //Step 1: Add vertices from the current level for the working face that will be used for subdivision.
        // Add the corners to vbuffer
        for (int i=0; i<(int)conn.size(); i++)
          {
            if (cur_level)
              vbuffer[i] = level_mesh[cur_level].start_vertex + (conn[i]-level_mesh[cur_level-1].start_vertex);
            else
              vbuffer[i] = level_mesh[cur_level].start_vertex + (conn[i] - *_inverts.begin());

            cur_conn.push_back(vbuffer[i]);
          }

        //Gather vertices already added to tracking array due to refinement of the sibling faces

        for (int i = 0; i < nepf; i++){
            for (int j = 0; j < nve; j++)
              {
                int id = refTemplates[findex][d].vert_on_edges[i][j];
                vbuffer[id] = trackvertsF[fid*nve*nepf+nve*i+j];
              }
          }

        //Add the remaining vertex handles to vbuffer for the current level for the working face
        for (int i=0; i<tnv; i++)
          {
            if (!vbuffer[i+nepf]){
                vbuffer[i+nepf] = level_mesh[cur_level].start_vertex+count_nverts;
                count_nverts += 1;
              }
          }

        //Step 2: Create the subentities using the template and the vbuffer
        int idx;       
        for (int i = 0; i < etotal; i++)
          {
            for (int k = 0; k < nepf; k++)
              {
                idx = refTemplates[findex][d].ents_conn[i][k];
                level_mesh[cur_level].face_conn[nepf*count_nents+k] = vbuffer[idx];
              }
            ent_buffer[i] = level_mesh[cur_level].start_face+count_nents;
            count_nents += 1;
          }

        // Step 3: Update the local AHF maps
        error = update_local_ahf(deg, ftype, vbuffer, ent_buffer, etotal); MB_CHK_ERR(error);

        //Step 4: Add the new vertices to the tracking array
        int id;

        for (int i = 0; i < nepf; i++)
          {
            // Add the vertices to trackvertsF for fid
            for (int j = 0; j < nve; j++)
              {
                id = refTemplates[findex][d].vert_on_edges[i][j];
                trackvertsF[fid*nepf*nve+nve*i+j] = vbuffer[id];
              }

            std::vector<EntityHandle> sibfids;
            std::vector<int> sibleids;
            std::vector<int> siborient;

            //Add the vertices to trackvertsF for siblings of fid, if any.
            error = ahf->get_up_adjacencies_2d(face, i, false, sibfids, &sibleids, &siborient); MB_CHK_ERR(error);

            if (!sibfids.size())
              continue;

            for (int s = 0; s < (int)sibfids.size(); s++)
              {
                int sibid;
                if (cur_level)
                  sibid = sibfids[s] - level_mesh[cur_level-1].start_face;
                else
                  sibid = sibfids[s] - *_infaces.begin();

                if (siborient[s]) // Same half-edge direction as the current half-edge
                  {
                    for (int j = 0; j < nve; j++)
                      {
                        id = refTemplates[findex][d].vert_on_edges[i][j];
                        trackvertsF[sibid*nepf*nve+nve*sibleids[s]+j] = vbuffer[id];
                      }
                  }
                else
                  {
                    for (int j = 0; j < nve; j++)
                      {
                        id = refTemplates[findex][d].vert_on_edges[i][nve-j-1];
                        trackvertsF[sibid*nepf*nve+nve*sibleids[s]+j] = vbuffer[id];
                      }
                  }
              }
          }

        //Step 5: Compute the coordinates of the new vertices, avoids computing more than once via the flag_verts array.
        double *corner_coords = new double[nepf*3];
        error = get_coordinates(&cur_conn[0], nepf, cur_level+1, corner_coords);  MB_CHK_ERR(error);

        error = compute_coordinates(cur_level, deg, ftype, vbuffer, vtotal, corner_coords, flag_verts, nverts_prev);  MB_CHK_ERR(error);

        delete [] corner_coords;
      }

    // Step 6: Update the global maps
    error = update_global_ahf(ftype, cur_level, deg);  MB_CHK_ERR(error);

    //Step 7: If edges exists, refine them.
    if (!_inedges.empty())
      {
        error = construct_hm_1D(cur_level, deg, ftype, trackvertsF); MB_CHK_ERR(error);
      }

    delete [] vbuffer;
    delete [] ent_buffer;

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::construct_hm_2D(int cur_level, int deg, EntityType type, std::vector<EntityHandle> &trackvertsE, std::vector<EntityHandle> &trackvertsF)
  {
    ErrorCode error;

    EntityType ftype=MBTRI;
    if (type == MBHEX)
      ftype = MBQUAD;

    int d = get_index_from_degree(deg);
    int findex = ftype-1;
    int cidx = ahf->get_index_in_lmap(*(_incells.begin()));

    int nepf = ahf->lConnMap2D[ftype-2].num_verts_in_face;
    int nepc = ahf->lConnMap3D[cidx].num_edges_in_cell;
    int nfpc = ahf->lConnMap3D[cidx].num_faces_in_cell;

    int tnv = refTemplates[findex][d].total_new_verts;
    int nve = refTemplates[findex][d].nv_edge;
    int nvf = refTemplates[findex][d].nv_face;
    int vtotal = nepf + tnv;
    int etotal = refTemplates[findex][d].total_new_ents;

    EntityHandle *vbuffer = new EntityHandle[vtotal];
    EntityHandle *ent_buffer = new EntityHandle[etotal];

    std::vector<EntityHandle> adjents, fconn, cconn;
    std::vector<int> leids;
    int count_nents = 0;

    int nents_prev, ecount;
     if (cur_level)
     {
       nents_prev = level_mesh[cur_level - 1].num_faces;
       ecount = level_mesh[cur_level - 1].num_edges * refTemplates[MBEDGE - 1][d].total_new_ents;;
     }
     else
     {
       nents_prev = _infaces.size();
       ecount  = _inedges.size() * refTemplates[MBEDGE - 1][d].total_new_ents;;
     }

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int it = 0; it < nents_prev; it++)
      {
        fconn.clear(); cconn.clear(); adjents.clear(); leids.clear();
        for (int i=0; i<vtotal; i++)
          vbuffer[i] = 0;
        for (int i=0; i<etotal; i++)
          ent_buffer[i] = 0;

        //EntityHandle of the working face
        EntityHandle face;
        if (cur_level)
          face = level_mesh[cur_level-1].start_face + it;
        else
          face = _infaces[it];

        error = get_connectivity(face, cur_level, fconn); MB_CHK_ERR(error);

        // Add the new handles for old connectivity in the buffer
        for (int i=0; i<(int)fconn.size(); i++)
          {
            if (cur_level)
              vbuffer[i] = level_mesh[cur_level].start_vertex + (fconn[i]-level_mesh[cur_level-1].start_vertex);
            else
              vbuffer[i] = level_mesh[cur_level].start_vertex + (fconn[i] - *_inverts.begin());
          }

        // Add handles for vertices on edges and faces from the already refined cell
        int fid, lid;
        error = ahf->get_up_adjacencies_face_3d(face, adjents, &leids); MB_CHK_ERR(error);

        if (cur_level)
          fid = adjents[0] - level_mesh[cur_level-1].start_cell;
        else
          fid = _incells.index(adjents[0]);

        lid = leids[0];

        error = get_connectivity(adjents[0], cur_level, cconn); MB_CHK_ERR(error);

        //Find the orientation w.r.t the half-face and then add vertices properly.
        EntityHandle *fac_conn = new EntityHandle[nepf];
        EntityHandle *lfac_conn = new EntityHandle[nepf];
        for (int j=0; j<nepf; j++)
          {
            fac_conn[j] = fconn[j];
            int id = ahf->lConnMap3D[cidx].hf2v[lid][j];
            lfac_conn[j] = cconn[id];
          }

        std::vector<int> le_idx, indices;

        error = reorder_indices(deg, fac_conn, lfac_conn, nepf, le_idx, indices); MB_CHK_ERR(error);

        delete [] fac_conn;
        delete [] lfac_conn;

        //Add the existing vertices on edges of the already refined cell to the vbuffer
        for (int j=0; j<nepf; j++)
          {
            int id = le_idx[j]; //Corresponding local edge
            int idx = ahf->lConnMap3D[cidx].f2leid[lid][id]; //Local edge in the cell

            //Get the orientation of the local edge of the face wrt the corresponding local edge in the cell
            bool eorient = false;
            int fnext = ahf->lConnMap2D[ftype-2].next[j];
            int idx1 = ahf->lConnMap3D[cidx].e2v[idx][0];
            int idx2 = ahf->lConnMap3D[cidx].e2v[idx][1];
            if ((fconn[j] == cconn[idx1] ) && (fconn[fnext] == cconn[idx2]))
              eorient = true;

            if (eorient){
                for (int k=0; k<nve; k++ )
                  {
                    int ind = refTemplates[findex][d].vert_on_edges[j][k];
                    vbuffer[ind] = trackvertsE[fid*nepc*nve+nve*idx+k];
                  }
              }
            else
              {
                for (int k=0; k<nve; k++ )
                  {
                    int ind = refTemplates[findex][d].vert_on_edges[j][nve-k-1];
                    vbuffer[ind] = trackvertsE[fid*nepc*nve+nve*idx+k];
                  }
              }
          }

        //Add the existing vertices on the face of the refine cell to vbuffer
        if (nvf)
          {
            for (int k=0; k<nvf; k++)
              {
                int ind = refTemplates[findex][d].vert_on_faces[0][k];
                vbuffer[ind] = trackvertsF[fid*nfpc*nvf+nvf*lid+indices[k]-1];
              }
          }

        // Create the subentities using the template and the vbuffer
        for (int i = 0; i < etotal; i++)
          {
            for (int k = 0; k < nepf; k++)
              {
                int idx = refTemplates[findex][d].ents_conn[i][k];
                level_mesh[cur_level].face_conn[nepf*count_nents+k] = vbuffer[idx];
              }
            ent_buffer[i] = level_mesh[cur_level].start_face+count_nents;
            count_nents += 1;
          }

        error = update_local_ahf(deg, ftype, vbuffer, ent_buffer, etotal); MB_CHK_ERR(error);

        //Create the interior edges
        int id1, id2;

        int ne = intFacEdg[ftype - 2][d].nie;
        for (int i = 0; i < ne; i++)
          {
            id1 = intFacEdg[ftype - 2][d].ieconn[i][0];
            id2 = intFacEdg[ftype - 2][d].ieconn[i][1];
            level_mesh[cur_level].edge_conn[2 * (ecount)] = vbuffer[id1];
            level_mesh[cur_level].edge_conn[2 * (ecount) + 1] = vbuffer[id2];      
            ecount += 1;
          }
      }


    // Step 6: Update the global maps
    error = update_global_ahf_2D_sub(cur_level, deg); MB_CHK_ERR(error);

    //Step 7: Update the hf-maps for the edges
    error = update_ahf_1D(cur_level);MB_CHK_ERR(error);

    delete [] vbuffer;
    delete [] ent_buffer;

    return MB_SUCCESS;

  }

  ErrorCode NestedRefine::construct_hm_3D(int cur_level, int deg)
  {
    ErrorCode error;
    EntityType type = mbImpl->type_from_handle(*(_incells.begin()));
    if (type == MBTET)
      {
        error = subdivide_tets(cur_level, deg);  MB_CHK_ERR(error);
      }
    else
      {
        error = subdivide_cells(type, cur_level, deg);  MB_CHK_ERR(error);
      }

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::subdivide_cells(EntityType type, int cur_level, int deg)
  {
    ErrorCode error;
    int nverts_prev, nents_prev;
    if (cur_level)
      {
        nverts_prev = level_mesh[cur_level-1].num_verts;
        nents_prev = level_mesh[cur_level-1].num_cells;
      }
    else
      {
        nverts_prev = _inverts.size();
        nents_prev = _incells.size();
      }

    //Create some book-keeping arrays over the parent mesh to avoid introducing duplicate vertices
    int cindex = type -1;
    int d = get_index_from_degree(deg);
    int ne = refTemplates[cindex][d].nv_edge;
    int nvf = refTemplates[cindex][d].nv_face;
    int nvtotal = refTemplates[cindex][d].total_new_verts;

    int index = ahf->get_index_in_lmap(*(_incells.begin()));
    int nvpc = ahf->lConnMap3D[index].num_verts_in_cell;
    int nepc = ahf->lConnMap3D[index].num_edges_in_cell;
    int nfpc = ahf->lConnMap3D[index].num_faces_in_cell;

    int vtotal = nvpc + nvtotal;
    EntityHandle *vbuffer = new EntityHandle[vtotal];

    std::vector<EntityHandle> trackvertsC_edg(nepc*ne*nents_prev, 0);
    std::vector<EntityHandle> trackvertsC_face(nfpc*nvf*nents_prev, 0);

    int cur_nverts = level_mesh[cur_level].num_verts;
    std::vector<int> flag_verts(cur_nverts-nverts_prev, 0);

    int count_nverts = nverts_prev;
    int count_ents = 0;
    std::vector<EntityHandle> conn, cur_conn;

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int cid = 0; cid < nents_prev; cid++)
      {
        conn.clear();
        cur_conn.clear();
        for (int i=0; i<vtotal; i++)
          vbuffer[i] = 0;

        //EntityHandle of the working cell
        EntityHandle cell;
        if (cur_level)
          cell = level_mesh[cur_level-1].start_cell + cid;
        else
          cell = _incells[cid];

        error = get_connectivity(cell, cur_level, conn);  MB_CHK_ERR(error);

        //Step 1: Add vertices from the current level for the working face that will be used for subdivision.
        // Add the corners to vbuffer
        for (int i=0; i<(int)conn.size(); i++)
          {
            if (cur_level)
              vbuffer[i] =  level_mesh[cur_level].start_vertex + (conn[i]-level_mesh[cur_level-1].start_vertex);
            else
              vbuffer[i] = level_mesh[cur_level].start_vertex +  (conn[i]-*_inverts.begin());

            cur_conn.push_back(vbuffer[i]);
          }

        //Gather vertices already added to tracking array due to refinement of the sibling cells
        for (int i=0; i<nepc; i++){
            for (int j=0; j<ne; j++)
              {
                int idx = refTemplates[cindex][d].vert_on_edges[i][j];
                vbuffer[idx] = trackvertsC_edg[cid*nepc*ne+ne*i+j];
              }
        }

        //Add remaining new vertex handles
        for (int i=0; i<nfpc; i++){
            for (int j=0; j<nvf; j++)
              {
                int idx = refTemplates[cindex][d].vert_on_faces[i][j];
                vbuffer[idx] = trackvertsC_face[cid*nfpc*nvf+nvf*i+j];
              }
          }

        //Add the remaining vertex handles to vbuffer for the current level for the working cell
        for (int i=0; i<nvtotal; i++){
            if (!vbuffer[i+nvpc]){
                vbuffer[i+nvpc] = level_mesh[cur_level].start_vertex+count_nverts;
                count_nverts += 1;
              }
          }

        //Step 2: Use the template to obtain the subentities. The coordinates and local ahf maps are also constructed.
        //Connectivity of the children
        int etotal = refTemplates[type-1][d].total_new_ents;
        EntityHandle  *ent_buffer = new EntityHandle[etotal];

        for (int i = 0; i < etotal; i++)
          {
            for (int k = 0; k < nvpc; k++)
              {
                int idx = refTemplates[type-1][d].ents_conn[i][k];
                level_mesh[cur_level].cell_conn[nvpc*count_ents+k] = vbuffer[idx];
              }
            ent_buffer[i] =  level_mesh[cur_level].start_cell + count_ents;
            count_ents += 1;
          }

        //Step 3: Update local ahf maps
        error = update_local_ahf(deg, type, vbuffer, ent_buffer, etotal);  MB_CHK_ERR(error);

        //Step 4: Update tracking information
        error = update_tracking_verts(cell, cur_level, deg, trackvertsC_edg, trackvertsC_face, vbuffer);  MB_CHK_ERR(error);

        //Step 5: Coordinates of the new vertices
        double *corner_coords = new double[nvpc*3];
        error = get_coordinates(&cur_conn[0], nvpc, cur_level+1, corner_coords); MB_CHK_ERR(error);

        error = compute_coordinates(cur_level, deg, type, vbuffer, vtotal, corner_coords, flag_verts, nverts_prev);  MB_CHK_ERR(error);

        delete [] ent_buffer;
        delete [] corner_coords;
      }

    //Step 6: Update the global maps
    error = update_global_ahf(type, cur_level, deg); MB_CHK_ERR(error);

    //Step 7: If edges exists, refine them as well.
    if (level_mesh[cur_level].num_edges != 0)
      {
        error = construct_hm_1D(cur_level,deg, type, trackvertsC_edg); MB_CHK_ERR(error);
      }

    //Step 8: If faces exists, refine them as well.
    if (!_infaces.empty())
      {
        error = construct_hm_2D(cur_level, deg, type, trackvertsC_edg, trackvertsC_face); MB_CHK_ERR(error);
      }

    delete [] vbuffer;

    //error = ahf->print_tags(3);

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::subdivide_tets(int cur_level, int deg)
  {
    ErrorCode error;
    int nverts_prev, nents_prev;
    if (cur_level)
      {
        nverts_prev = level_mesh[cur_level-1].num_verts;
        nents_prev = level_mesh[cur_level-1].num_cells;
      }
    else
      {
        nverts_prev = _inverts.size();
        nents_prev = _incells.size();
      }

    EntityType type = MBTET;
    int cindex = type -1;
    int d = get_index_from_degree(deg);
    int ne = refTemplates[cindex][d].nv_edge;
    int nvf = refTemplates[cindex][d].nv_face;
    int nvtotal = refTemplates[cindex][d].total_new_verts;

    int index = ahf->get_index_in_lmap(*(_incells.begin()));
    int nvpc = ahf->lConnMap3D[index].num_verts_in_cell;
    int nepc = ahf->lConnMap3D[index].num_edges_in_cell;
    int nfpc = ahf->lConnMap3D[index].num_faces_in_cell;

    // Create vertex buffer
    int vtotal = nvpc + nvtotal;
    EntityHandle *vbuffer = new EntityHandle[vtotal];

    //Create book-keeping arrays over the parent mesh to avoid introducing duplicate vertices
    std::vector<EntityHandle> trackvertsC_edg(nepc*ne*nents_prev, 0);
    std::vector<EntityHandle> trackvertsC_face(nfpc*nvf*nents_prev, 0);

    int cur_nverts = level_mesh[cur_level].num_verts;
    std::vector<int> flag_verts(cur_nverts-nverts_prev, 0);
    std::vector<int> cell_patterns(nents_prev,0);

    int count_nverts = nverts_prev;
    int count_ents = 0;
    std::vector<EntityHandle> conn, cur_conn;

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int cid = 0; cid < nents_prev; cid++)
      {
        conn.clear();
        cur_conn.clear();
        for (int i=0; i<vtotal; i++)
          vbuffer[i] = 0;

        //EntityHandle of the working cell
        EntityHandle cell;
        if (cur_level)
          cell = level_mesh[cur_level-1].start_cell + cid;
        else
          cell = _incells[cid];

        error = get_connectivity(cell, cur_level, conn);  MB_CHK_ERR(error);

        //Step 1: Add vertices from the current level for the working face that will be used for subdivision.
        // Add the corners to vbuffer
        for (int i=0; i<(int)conn.size(); i++)
          {
            if (cur_level)
              vbuffer[i] =  level_mesh[cur_level].start_vertex + (conn[i]-level_mesh[cur_level-1].start_vertex);
            else
              vbuffer[i] = level_mesh[cur_level].start_vertex +  (conn[i]-*_inverts.begin()) ;

             cur_conn.push_back(vbuffer[i]);
          }

        //Gather vertices already added to tracking array due to refinement of the sibling cells
        for (int i=0; i<nepc; i++){
            for (int j=0; j<ne; j++)
              {
                int idx = refTemplates[cindex][d].vert_on_edges[i][j];
                vbuffer[idx] = trackvertsC_edg[cid*nepc*ne+ne*i+j];
              }
        }

        //Add remaining new vertex handles
        for (int i=0; i<nfpc; i++){
            for (int j=0; j<nvf; j++)
              {
                int idx = refTemplates[cindex][d].vert_on_faces[i][j];
                vbuffer[idx] = trackvertsC_face[cid*nfpc*nvf+nvf*i+j];
              }
          }

        //Add the remaining vertex handles to vbuffer for the current level for the working cell
        for (int i=0; i<nvtotal; i++){
            if (!vbuffer[i+nvpc]){
                vbuffer[i+nvpc] = level_mesh[cur_level].start_vertex+count_nverts;
                count_nverts += 1;
              }
          }

        //Step 2: Coordinates of the new vertices
        double *corner_coords = new double[nvpc*3];
        error = get_coordinates(&cur_conn[0], nvpc, cur_level+1, corner_coords);  MB_CHK_ERR(error);

        error = compute_coordinates(cur_level, deg, type, vbuffer, vtotal, corner_coords, flag_verts, nverts_prev);  MB_CHK_ERR(error);

        //Step 3: Choose the tet refine pattern to be used for this tet
        int diag = find_shortest_diagonal_octahedron(cur_level, deg, vbuffer);
        int pat_id = diag + 2;
        cell_patterns[cid] = pat_id;

        //Step 4: Use the template to obtain the subentities. The coordinates and local ahf maps are also constructed.
        //Connectivity of the children
        int etotal = refTemplates[pat_id][d].total_new_ents;
        EntityHandle  *ent_buffer = new EntityHandle[etotal];

        for (int i = 0; i < etotal; i++)
          {
            for (int k = 0; k < nvpc; k++)
              {
                int idx = refTemplates[pat_id][d].ents_conn[i][k];
                level_mesh[cur_level].cell_conn[nvpc*count_ents+k] = vbuffer[idx];
              }
            ent_buffer[i] =  level_mesh[cur_level].start_cell + count_ents;
            count_ents += 1;
          }

        //Step 5: Update local ahf maps
        error = update_local_ahf(deg, MBTET, pat_id, vbuffer, ent_buffer, etotal); MB_CHK_ERR(error);

        //Step 6: Update tracking information
        error = update_tracking_verts(cell, cur_level, deg, trackvertsC_edg, trackvertsC_face, vbuffer);  MB_CHK_ERR(error);

        delete [] ent_buffer;
        delete [] corner_coords;
      }

    //Step 7: Update the global maps
    error = update_global_ahf(cur_level, deg, cell_patterns); MB_CHK_ERR(error);

    //Step 8: If edges exists, refine them as well.
    if (level_mesh[cur_level].num_edges != 0)
      {
        error = construct_hm_1D(cur_level,deg, type, trackvertsC_edg); MB_CHK_ERR(error);
      }

    //Step 9: If faces exists, refine them as well.
    if (!_infaces.empty())
      {
        error = construct_hm_2D(cur_level, deg, type, trackvertsC_edg, trackvertsC_face); MB_CHK_ERR(error);
      }

    delete [] vbuffer;

    return MB_SUCCESS;
  }


  ErrorCode NestedRefine::compute_coordinates(int cur_level, int deg, EntityType type, EntityHandle *vbuffer, int vtotal, double *corner_coords, std::vector<int> &vflag, int nverts_prev)
  {
    EntityHandle vstart = level_mesh[cur_level].start_vertex;
    int d = get_index_from_degree(deg);

    if (type == MBTRI)
      {
        double xi, eta, N[3];
        int findex = mbImpl->type_from_handle(*(_infaces.begin()))-1;

        for (int i=3; i<vtotal; i++)
          {
            if  (vflag[vbuffer[i]-vstart - nverts_prev])
              continue;

            xi = refTemplates[findex][d].vert_nat_coord[i-3][0];
            eta = refTemplates[findex][d].vert_nat_coord[i-3][1];
            N[0] = 1-xi-eta; N[1] = xi; N[2] = eta;

            double x =0, y = 0, z = 0;
            for (int j=0; j<3; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;
            vflag[vbuffer[i]-vstart - nverts_prev] = 1;
          }
      }
      else if (type == MBQUAD)
      {
        double xi, eta, N[4];
        int findex = mbImpl->type_from_handle(*(_infaces.begin()))-1;

        for (int i=4; i<vtotal; i++)
          {
            if (vflag[vbuffer[i]-vstart-nverts_prev])
              continue;

            xi = refTemplates[findex][d].vert_nat_coord[i-4][0];
            eta = refTemplates[findex][d].vert_nat_coord[i-4][1];
            N[0] = (1-xi)*(1-eta)/4; N[1] = (1+xi)*(1-eta)/4; N[2] = (1+xi)*(1+eta)/4, N[3] = (1-xi)*(1+eta)/4;

            double x =0, y = 0, z = 0;
            for (int j=0; j<4; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;
            vflag[vbuffer[i]-vstart-nverts_prev] = 1;
          }

      }
      else if (type == MBTET)
      {
        double xi, eta, mu, N[4];
        int cindex = mbImpl->type_from_handle(*(_incells.begin()))-1;

        for (int i=4; i<vtotal; i++)
          {
            if (vflag[vbuffer[i]-vstart-nverts_prev])
              continue;

            xi = refTemplates[cindex][d].vert_nat_coord[i-4][0];
            eta = refTemplates[cindex][d].vert_nat_coord[i-4][1];
            mu = refTemplates[cindex][d].vert_nat_coord[i-4][2];

            N[0] = 1-xi-eta-mu; N[1] = xi; N[2] = eta, N[3] = mu;

            double x =0, y = 0, z = 0;
            for (int j=0; j<4; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;
            vflag[vbuffer[i]-vstart-nverts_prev] = 1;
          }

      }
      else if (type == MBPRISM)
      {
        double xi, eta, mu, N[6];
        int cindex = mbImpl->type_from_handle(*(_incells.begin()))-1;

        for (int i=6; i<vtotal; i++)
          {
            if (vflag[vbuffer[i]-vstart-nverts_prev])
              continue;

            xi = refTemplates[cindex][d].vert_nat_coord[i-6][0];
            eta = refTemplates[cindex][d].vert_nat_coord[i-6][1];
            mu = refTemplates[cindex][d].vert_nat_coord[i-6][2];

            N[0] = (1-xi-eta)*(1-mu), N[1] = xi*(1-mu), N[2] = eta*(1-mu), N[3] = (1-xi-eta)*(1+mu), N[4] = xi*(1+mu), N[5] = eta*(1+mu);

            double x =0, y = 0, z = 0;
            for (int j=0; j<6; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;
            vflag[vbuffer[i]-vstart-nverts_prev] = 1;
          }

      }
      else if (type == MBHEX)
      {
        double xi, eta, mu, N[8];
        double d1,d2,d3, s1,s2,s3;
        int cindex = mbImpl->type_from_handle(*(_incells.begin()))-1;

        for (int i=8; i<vtotal; i++)
          {

            if (vflag[vbuffer[i] - vstart - nverts_prev])
              continue;

            xi = refTemplates[cindex][d].vert_nat_coord[i-8][0];
            eta = refTemplates[cindex][d].vert_nat_coord[i-8][1];
            mu = refTemplates[cindex][d].vert_nat_coord[i-8][2];

            d1 = 1-xi; d2 = 1- eta; d3 = 1-mu;
            s1 = 1+xi; s2 = 1+eta; s3 = 1+mu;
            N[0] = (d1*d2*d3)/8; N[1] = (s1*d2*d3)/8; N[2] = (s1*s2*d3)/8; N[3] = (d1*s2*d3)/8;
            N[4] = (d1*d2*s3)/8; N[5] = (s1*d2*s3)/8; N[6] =  (s1*s2*s3)/8; N[7] = (d1*s2*s3)/8;

            double x =0, y = 0, z = 0;
            for (int j=0; j<8; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;


            vflag[vbuffer[i]-vstart-nverts_prev] = 1;
          }
      }
    return MB_SUCCESS;
  }

  /**********************************
   *          Update AHF maps           *
   * ********************************/

ErrorCode NestedRefine::update_local_ahf(int deg, EntityType type, int pat_id, EntityHandle *vbuffer, EntityHandle *ent_buffer, int etotal)
{
  ErrorCode error;
  int nhf = 0, nv = 0, total_new_verts = 0;
  int d = get_index_from_degree(deg);

  //Get the number of half-facets
  if (type == MBEDGE)
    {
      nhf = 2;
      nv = 2;
      total_new_verts = refTemplates[0][d].total_new_verts;
    }
  else if (type == MBTRI || type == MBQUAD)
    {
      nhf = ahf->lConnMap2D[type-2].num_verts_in_face;
      nv = nhf;
      total_new_verts = refTemplates[pat_id][d].total_new_verts;
    }
  else if (type == MBTET || type == MBHEX)
    {
      int index = ahf->get_index_in_lmap(*_incells.begin());
      nhf = ahf->lConnMap3D[index].num_faces_in_cell;
      nv =  ahf->lConnMap3D[index].num_verts_in_cell;
      total_new_verts = refTemplates[pat_id][d].total_new_verts;
    }

  std::vector<EntityHandle> ent;
  std::vector<int> lid;

  //Update the vertex to half-facet map
  for (int i=0; i<total_new_verts; i++)
    {
      ent.clear();lid.clear();
      EntityHandle vid = vbuffer[i+nv];
      error = ahf->get_incident_map(type, vid, ent, lid);  MB_CHK_ERR(error);

      if (ent[0])
        continue;

      int id = refTemplates[pat_id][d].v2hf[i+nv][0]-1;
      ent[0] = ent_buffer[id];
      lid[0] = refTemplates[pat_id][d].v2hf[i+nv][1];

      error = ahf->set_incident_map(type, vid, ent, lid);  MB_CHK_ERR(error);
    }

  //Update the sibling half-facet map
  for (int i=0; i< etotal; i++)
    {
      EntityHandle  *sib_entids = new EntityHandle[nhf];
      int *sib_lids = new int[nhf];

      error = ahf->get_sibling_map(type, ent_buffer[i], &sib_entids[0], &sib_lids[0], nhf);  MB_CHK_ERR(error);

      for (int l=0; l< nhf; l++)
        {
          if (sib_entids[l])
            continue;

          // Fill out the sibling values
          int id = refTemplates[pat_id][d].ents_opphfs[i][2*l];
          if (id)
            {
              sib_entids[l] = ent_buffer[id-1];
              sib_lids[l] = refTemplates[pat_id][d].ents_opphfs[i][2*l+1];
            }
          else
            {
              sib_entids[l] = 0;
              sib_lids[l] = 0;
            }
        }

      error = ahf->set_sibling_map(type, ent_buffer[i], &sib_entids[0], &sib_lids[0], nhf); MB_CHK_ERR(error);

      for (int l=0; l< nhf; l++)
        {
          if (sib_entids[l]){

              EntityHandle set_entid = ent_buffer[i];
              int set_lid = l;

              error = ahf->set_sibling_map(type, sib_entids[l], sib_lids[l], set_entid, set_lid);  MB_CHK_ERR(error);
            }
        }

      delete [] sib_entids;
      delete [] sib_lids;

    }
  return MB_SUCCESS;
}

ErrorCode NestedRefine::update_local_ahf(int deg, EntityType type, EntityHandle *vbuffer, EntityHandle *ent_buffer, int etotal)
{
  ErrorCode error;
  assert(type != MBTET);
  error = update_local_ahf(deg, type, type-1, vbuffer, ent_buffer, etotal); MB_CHK_ERR(error);

  return MB_SUCCESS;
}

ErrorCode NestedRefine::update_global_ahf(EntityType type, int cur_level, int deg)
{

  ErrorCode error;

  //Get the number of half-facets and number of children of each type
  if (type == MBEDGE)
    {
      error = update_global_ahf_1D(cur_level, deg); MB_CHK_ERR(error);
    }
  else if (type == MBTRI || type == MBQUAD)
    {
      error = update_global_ahf_2D(cur_level, deg); MB_CHK_ERR(error);
    }
  else if (type == MBHEX)
    {
      error = update_global_ahf_3D(cur_level, deg); MB_CHK_ERR(error);
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::update_global_ahf(int cur_level, int deg, std::vector<int> &pattern_ids)
{
  ErrorCode error;
  error = update_global_ahf_3D(cur_level, deg, pattern_ids); MB_CHK_ERR(error);

  return MB_SUCCESS;
}

ErrorCode NestedRefine::update_global_ahf_1D(int cur_level, int deg)
 {
   ErrorCode error;
   int d = get_index_from_degree(deg);
   int nhf, nchilds, nverts_prev, nents_prev;
   nhf = 2;
   nchilds = refTemplates[0][d].total_new_ents;
   if (cur_level)
     {
       nverts_prev = level_mesh[cur_level-1].num_verts;
       nents_prev = level_mesh[cur_level-1].num_edges;
     }
   else
     {
       nverts_prev = _inverts.size();
       nents_prev = _inedges.size();
     }

   std::vector<EntityHandle> inci_ent, child_ents;
   std::vector<int> inci_lid, child_lids;

   //Update the vertex to half-facet maps for duplicate vertices
   for (int i=0; i<nverts_prev; i++)
     {
       inci_ent.clear(); inci_lid.clear(); child_ents.clear(); child_lids.clear();

       //Vertex id in the previous mesh and the current one
       EntityHandle vid;
       if (cur_level)
         vid = level_mesh[cur_level-1].start_vertex + i;
       else
         vid = _inverts[i];
       EntityHandle cur_vid = level_mesh[cur_level].start_vertex + i;

       //Get the incident half-vert in the previous mesh     
       error = ahf->get_incident_map(MBEDGE, vid, inci_ent, inci_lid);  MB_CHK_ERR(error);

       // Obtain the corresponding incident child in the current mesh
       int lvid = get_local_vid(vid, inci_ent[0], cur_level-1);
       int chid = refTemplates[0][d].v2hf[lvid][0]-1;

       int pid;
       if (cur_level)
         pid = inci_ent[0] - level_mesh[cur_level-1].start_edge;
       else
         pid = inci_ent[0] - *_inedges.begin();

       int ind = nchilds*pid;

       child_ents.push_back(level_mesh[cur_level].start_edge + ind+chid);
       child_lids.push_back(refTemplates[0][d].v2hf[lvid][1]);

       error = ahf->set_incident_map(MBEDGE, cur_vid, child_ents, child_lids);  MB_CHK_ERR(error);
     }

   //Update the sibling half-facet maps across entities
   for (int i=0; i< nents_prev; i++)
     {
       EntityHandle ent;
       if (cur_level)
         ent = level_mesh[cur_level-1].start_edge + i;
       else
         ent = _inedges[i];

       EntityHandle *sib_entids = new EntityHandle[nhf];
       int *sib_lids = new int[nhf];

       error = ahf->get_sibling_map(MBEDGE, ent, &sib_entids[0], &sib_lids[0], nhf);  MB_CHK_ERR(error);

       int id, idx;

       for (int l=0; l < nhf; l++)
         {
           if (!sib_entids[l])
             continue;

           //Find the child incident on the half-facet
           id = refTemplates[0][d].ents_on_pent[l][1]-1;
           idx = nchilds*i;
           EntityHandle child_ent = level_mesh[cur_level].start_edge + idx+id ;
           int ch_lid = l;

           //Find the sibling of the child
           EntityHandle *sib_childs = new EntityHandle[nhf];
           int *sib_chlids = new int[nhf];

           error = ahf->get_sibling_map(MBEDGE, child_ent, &sib_childs[0], &sib_chlids[0], nhf);  MB_CHK_ERR(error);

           //If the sibling already exists, dont do anything
           if (sib_childs[ch_lid])
             continue;

           //Get the correponding child of the sibling of the current parent
           int psib;
           if (cur_level)
             psib = sib_entids[l] - level_mesh[cur_level-1].start_edge;
           else
             psib = sib_entids[l] - *_inedges.begin();

           int plid = sib_lids[l];

           id = refTemplates[0][d].ents_on_pent[plid][1]-1;
           idx = nchilds*psib;

           EntityHandle psib_child = level_mesh[cur_level].start_edge + idx+id ;
           int psib_chlid = plid;

           //Set the siblings
           sib_childs[ch_lid] = psib_child;
           sib_chlids[ch_lid] = psib_chlid;

           error = ahf->set_sibling_map(MBEDGE, child_ent,  &sib_childs[0], &sib_chlids[0], nhf);  MB_CHK_ERR(error);

           delete [] sib_childs;
           delete [] sib_chlids;
         }
       delete [] sib_entids;
       delete [] sib_lids;
     }

   return MB_SUCCESS;
 }

ErrorCode NestedRefine::update_global_ahf_1D_sub(int cur_level, int deg)
{
  ErrorCode error;
  int d = get_index_from_degree(deg);
  int nhf, nchilds, nents_prev;
  nhf = 2;
  nchilds = refTemplates[0][d].total_new_ents;
  if (cur_level)
    {
      nents_prev = level_mesh[cur_level-1].num_edges;
    }
  else
    {
      nents_prev = _inedges.size();
    }

  //Update the sibling half-facet maps across entities

  std::vector<EntityHandle> conn;
  for (int i=0; i< nents_prev; i++)
    {
      EntityHandle ent;
      if (cur_level)
        ent = level_mesh[cur_level-1].start_edge + i;
      else
        ent = _inedges[i];

      //Set incident hv maps
      conn.clear();
      error = get_connectivity(ent,cur_level,conn); MB_CHK_ERR(error);

      std::vector<EntityHandle> inci_ent, child_ents;
      std::vector<int> inci_lid, child_lids;
      for (int j=0; j<2; j++)
        {
          inci_ent.clear(); inci_lid.clear(); child_ents.clear(); child_lids.clear();

          // Get the entityhandle of the vertex from previous level in the current level
          EntityHandle cur_vid;
          if (cur_level)
            cur_vid = level_mesh[cur_level].start_vertex + (conn[j]-level_mesh[cur_level-1].start_vertex);
          else
            cur_vid = level_mesh[cur_level].start_vertex + (conn[j] - *_inverts.begin());

          //Obtain the incident half-facet. If exists, then no need to assign another
          error = ahf->get_incident_map(MBEDGE, cur_vid, inci_ent, inci_lid);  MB_CHK_ERR(error);
          if (inci_ent[0] != 0)
            continue;

          //Get the incident half-facet on the old vertex
          error = ahf->get_incident_map(MBEDGE, conn[j], inci_ent, inci_lid);  MB_CHK_ERR(error);

          // Obtain the corresponding incident child in the current mesh
          int lvid = get_local_vid(conn[j], inci_ent[0], cur_level-1);
          int chid = refTemplates[0][d].v2hf[lvid][0]-1;

          int pid;
          if (cur_level)
            pid = inci_ent[0] - level_mesh[cur_level-1].start_edge;
          else
            pid = inci_ent[0] - *_inedges.begin();

          int ind = nchilds*pid;

          child_ents.push_back(level_mesh[cur_level].start_edge + ind+chid);
          child_lids.push_back(refTemplates[0][d].v2hf[lvid][1]);

          error = ahf->set_incident_map(MBEDGE, cur_vid, child_ents, child_lids);  MB_CHK_ERR(error);
        }

      EntityHandle *sib_entids = new EntityHandle[nhf];
      int *sib_lids = new int[nhf];

      error = ahf->get_sibling_map(MBEDGE, ent, &sib_entids[0], &sib_lids[0], nhf);  MB_CHK_ERR(error);

      int id, idx;

      for (int l=0; l < nhf; l++)
        {
          if (!sib_entids[l])
            continue;

          //Find the child incident on the half-facet
          id = refTemplates[0][d].ents_on_pent[l][1]-1;
          idx = nchilds*i;
          EntityHandle child_ent = level_mesh[cur_level].start_edge + idx+id ;
          int ch_lid = l;

          //Find the sibling of the child
          EntityHandle *sib_childs = new EntityHandle[nhf];
          int *sib_chlids = new int[nhf];

          error = ahf->get_sibling_map(MBEDGE, child_ent, &sib_childs[0], &sib_chlids[0], nhf);  MB_CHK_ERR(error);

          //If the sibling already exists, dont do anything
          if (sib_childs[ch_lid])
            continue;

          //Get the correponding child of the sibling of the current parent
          int psib;
          if (cur_level)
            psib = sib_entids[l] - level_mesh[cur_level-1].start_edge;
          else
            psib = sib_entids[l] - *_inedges.begin();

          int plid = sib_lids[l];

          id = refTemplates[0][d].ents_on_pent[plid][1]-1;
          idx = nchilds*psib;

          EntityHandle psib_child = level_mesh[cur_level].start_edge + idx+id ;
          int psib_chlid = plid;

          //Set the siblings
          sib_childs[ch_lid] = psib_child;
          sib_chlids[ch_lid] = psib_chlid;

          error = ahf->set_sibling_map(MBEDGE, child_ent, &sib_childs[0], &sib_chlids[0], nhf);  MB_CHK_ERR(error);

          delete [] sib_childs;
          delete [] sib_chlids;
        }
      delete [] sib_entids;
      delete [] sib_lids;
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::update_ahf_1D(int cur_level)
{
  ErrorCode error;
  error = ahf->determine_sibling_halfverts(level_mesh[cur_level].verts, level_mesh[cur_level].edges);MB_CHK_ERR(error);

  error = ahf->determine_incident_halfverts( level_mesh[cur_level].edges);MB_CHK_ERR(error);

  return MB_SUCCESS;
}

 ErrorCode NestedRefine::update_global_ahf_2D(int cur_level, int deg)
 {
   ErrorCode error;

   EntityType type = mbImpl->type_from_handle(*_infaces.begin());
   int nhf, nchilds, nverts_prev, nents_prev;

   nhf = ahf->lConnMap2D[type-2].num_verts_in_face;
   int d = get_index_from_degree(deg);
   nchilds = refTemplates[type-1][d].total_new_ents;

   if (cur_level)
     {
       nverts_prev = level_mesh[cur_level-1].num_verts;
       nents_prev = level_mesh[cur_level-1].num_faces;
     }
   else
     {
       nverts_prev = _inverts.size();
       nents_prev = _infaces.size();
     }

   std::vector<EntityHandle> inci_ent, child_ents;
   std::vector<int> inci_lid, child_lids;

   //Update the vertex to half-edge maps for old/duplicate vertices
   for (int i=0; i<nverts_prev; i++)
     {
       inci_ent.clear(); inci_lid.clear(); child_ents.clear(); child_lids.clear();

       //Vertex id in the previous mesh
       EntityHandle vid;
       if (cur_level)
         vid = level_mesh[cur_level-1].start_vertex + i;
       else
         vid = _inverts[i];
       EntityHandle cur_vid = level_mesh[cur_level].start_vertex + i;

       //Get the incident half-vert in the previous mesh
       error = ahf->get_incident_map(type, vid, inci_ent, inci_lid); MB_CHK_ERR(error);

       // Obtain the corresponding incident child in the current mesh
       for (int j=0; j<(int) inci_ent.size(); j++){
           int lvid = get_local_vid(vid, inci_ent[j], cur_level-1);
           int chid = refTemplates[type-1][d].v2hf[lvid][0]-1;

           int pid;
           if (cur_level)
             pid = inci_ent[j] - level_mesh[cur_level-1].start_face;
           else
             pid = inci_ent[j] - *_infaces.begin();

           int ind = nchilds*pid;

           child_ents.push_back(level_mesh[cur_level].start_face + ind+chid) ;
           child_lids.push_back(refTemplates[type-1][d].v2hf[lvid][1]);
         }
       error = ahf->set_incident_map(type, cur_vid, child_ents, child_lids);  MB_CHK_ERR(error);
     }


   EntityHandle fedge[2];

   //Update the sibling half-facet maps across entities
   for (int i=0; i< nents_prev; i++)
     {
       EntityHandle ent;
       if (cur_level)
         ent = level_mesh[cur_level-1].start_face + i;
       else
         ent = _infaces[i];

       std::vector<EntityHandle> fid_conn;
       error = get_connectivity(ent, cur_level, fid_conn);
       if (MB_SUCCESS != error) return error;

       EntityHandle *sib_entids = new EntityHandle[nhf];
       int *sib_lids = new int[nhf];

       error = ahf->get_sibling_map(type, ent, &sib_entids[0], &sib_lids[0], nhf);  MB_CHK_ERR(error);

       int id, idx;

       for (int l=0; l < nhf; l++)
         {
           if (!sib_entids[l])
             continue;

           int nidx = ahf->lConnMap2D[type-2].next[l];
           fedge[0] = fid_conn[l];
           fedge[1] = fid_conn[nidx];

           EntityHandle sfid = sib_entids[l];
           int slid = sib_lids[l];

           std::vector<EntityHandle> conn;
           error = get_connectivity(sfid, cur_level, conn);
           if (MB_SUCCESS != error) return error;

           bool orient = true;
           nidx = ahf->lConnMap2D[type-2].next[slid];
           if ((fedge[1] == conn[slid])&&(fedge[0] == conn[nidx]))
             orient = false;

           if (orient)
             assert((fedge[0] == conn[slid])&&(fedge[1] == conn[nidx]));

           //Find the childrens incident on the half-facet
           int nch = refTemplates[type-1][d].ents_on_pent[l][0];
           idx = nchilds*i;

           //Loop over all the incident childrens
           for (int k=0; k<nch; k++)
             {
               id = refTemplates[type-1][d].ents_on_pent[l][k+1]-1;
               EntityHandle child_ent = level_mesh[cur_level].start_face + idx+id ;
               int child_lid = l;

               //Find the sibling of the child
               EntityHandle child_sibent;
               int child_siblid;
               error = ahf->get_sibling_map(type, child_ent, child_lid, child_sibent, child_siblid); MB_CHK_ERR(error);

               if (child_sibent != 0)
                 continue;

               //Get the correponding child of the sibling of the current parent
               int psib;
               if (cur_level)
                 psib = sfid- level_mesh[cur_level-1].start_face;
               else
                 psib = sfid - *_infaces.begin();

               int plid = slid;

               if (orient)
                 id = refTemplates[type-1][d].ents_on_pent[plid][k+1]-1;
               else
                 id = refTemplates[type-1][d].ents_on_pent[plid][nch-k]-1;

              int sidx = nchilds*psib;

               EntityHandle psib_child = level_mesh[cur_level].start_face + sidx+id ;
               int psib_chlid = plid;

               //Set the siblings
               error = ahf->set_sibling_map(type, child_ent, child_lid, psib_child, psib_chlid); MB_CHK_ERR(error);
             }
         }

       delete [] sib_entids;
       delete [] sib_lids;
     }

   return MB_SUCCESS;

 }

 ErrorCode NestedRefine::update_global_ahf_2D_sub(int cur_level, int deg)
 {
   ErrorCode error;
   int d = get_index_from_degree(deg);
   EntityType type = mbImpl->type_from_handle(*_infaces.begin());
   int nhf, nchilds, nents_prev;
   nhf = ahf->lConnMap2D[type-2].num_verts_in_face;
   nchilds = refTemplates[type-1][d].total_new_ents;

   if (cur_level)
       nents_prev = level_mesh[cur_level-1].num_faces;
   else
       nents_prev = _infaces.size();

   EntityHandle fedge[2];

   //Update the sibling half-facet maps across entities
   for (int i=0; i< nents_prev; i++)
     {
       EntityHandle ent;
       if (cur_level)
         ent = level_mesh[cur_level-1].start_face + i;
       else
         ent = _infaces[i];

       std::vector<EntityHandle> fid_conn;
       error = get_connectivity(ent, cur_level, fid_conn);
       if (MB_SUCCESS != error) return error;

       std::vector<EntityHandle> inci_ent, child_ents;
       std::vector<int> inci_lid, child_lids;

       //Set incident half-edges
       for (int j=0; j<nhf; j++)
         {
           inci_ent.clear(); inci_lid.clear(); child_ents.clear(); child_lids.clear();
           EntityHandle cur_vid;
           if (cur_level)
             cur_vid = level_mesh[cur_level].start_vertex + (fid_conn[j]-level_mesh[cur_level-1].start_vertex);
           else
             cur_vid = level_mesh[cur_level].start_vertex + (fid_conn[j] - *_inverts.begin());

           //Obtain the incident half-facet. If exists, then no need to assign another
           error = ahf->get_incident_map(type, cur_vid, inci_ent, inci_lid);  MB_CHK_ERR(error);
           if (inci_ent[0] != 0)
             continue;

           //Get the incident half-facet on the old vertex
           error = ahf->get_incident_map(type, fid_conn[j], inci_ent, inci_lid); MB_CHK_ERR(error);

           // Obtain the corresponding incident child in the current mesh
           for (int k=0; k<(int)inci_ent.size(); k++){
               int lvid = get_local_vid(fid_conn[j], inci_ent[k], cur_level-1);
               int chid = refTemplates[type-1][d].v2hf[lvid][0]-1;

               int pid;
               if (cur_level)
                 pid = inci_ent[k] - level_mesh[cur_level-1].start_face;
               else
                 pid = inci_ent[k] - *_infaces.begin();

               int ind = nchilds*pid;

               child_ents.push_back(level_mesh[cur_level].start_face + ind+chid);
               child_lids.push_back(refTemplates[type-1][d].v2hf[lvid][1]);
             }

           error = ahf->set_incident_map(type, cur_vid, child_ents, child_lids);  MB_CHK_ERR(error);
         }

       //Set sibling half-edges
       EntityHandle *sib_entids = new EntityHandle[nhf];
       int *sib_lids = new int[nhf];

       error = ahf->get_sibling_map(type, ent, &sib_entids[0], &sib_lids[0], nhf);  MB_CHK_ERR(error);

       int id, idx;

       for (int l=0; l < nhf; l++)
         {
           if (!sib_entids[l])
             continue;

           int nidx = ahf->lConnMap2D[type-2].next[l];
           fedge[0] = fid_conn[l];
           fedge[1] = fid_conn[nidx];

           EntityHandle sfid = sib_entids[l];
           int slid = sib_lids[l];

           std::vector<EntityHandle> conn;
           error = get_connectivity(sfid, cur_level, conn);MB_CHK_ERR(error);

           assert((int)conn.size() > nidx && (int)conn.size() > slid);


           bool orient = true;
           nidx = ahf->lConnMap2D[type-2].next[slid];
           if ((fedge[1] == conn[slid])&&(fedge[0] == conn[nidx]))
             orient = false;

           if (orient)
             assert((fedge[0] == conn[slid])&&(fedge[1] == conn[nidx]));

           //Find the childrens incident on the half-facet
           int nch = refTemplates[type-1][d].ents_on_pent[l][0];
           idx = nchilds*i;

           //Loop over all the incident childrens
           for (int k=0; k<nch; k++)
             {
               id = refTemplates[type-1][d].ents_on_pent[l][k+1]-1;
               EntityHandle child_ent = level_mesh[cur_level].start_face + idx+id ;
               int child_lid = l;

               //Find the sibling of the child
               EntityHandle child_sibent;
               int child_siblid;
               error = ahf->get_sibling_map(type, child_ent, child_lid, child_sibent, child_siblid); MB_CHK_ERR(error);

               if (child_sibent != 0)
                 continue;

               //Get the correponding child of the sibling of the current parent
               int psib;
               if (cur_level)
                 psib = sfid- level_mesh[cur_level-1].start_face;
               else
                 psib = sfid - *_infaces.begin();

               int plid = slid;

               if (orient)
                 id = refTemplates[type-1][d].ents_on_pent[plid][k+1]-1;
               else
                 id = refTemplates[type-1][d].ents_on_pent[plid][nch-k]-1;

              int sidx = nchilds*psib;

               EntityHandle psib_child = level_mesh[cur_level].start_face + sidx+id ;
               int psib_chlid = plid;

               //Set the siblings
               error = ahf->set_sibling_map(type, child_ent, child_lid, psib_child, psib_chlid); MB_CHK_ERR(error);
             }
         }

       delete [] sib_entids;
       delete [] sib_lids;
     }

   return MB_SUCCESS;

 }


ErrorCode NestedRefine::update_global_ahf_3D(int cur_level, int deg)
{
  ErrorCode error;
  int nvpc, ne, nhf, nchilds, nverts_prev, nents_prev;

  EntityType type = mbImpl->type_from_handle(*_incells.begin());
  int index = ahf->get_index_in_lmap(*_incells.begin());
  int d = get_index_from_degree(deg);

  nhf = ahf->lConnMap3D[index].num_faces_in_cell;
  ne = ahf->lConnMap3D[index].num_edges_in_cell;
  nvpc = ahf->lConnMap3D[index].num_verts_in_cell;
  nchilds = refTemplates[type-1][d].total_new_ents;

  if (cur_level)
    {
      nverts_prev = level_mesh[cur_level-1].num_verts;
      nents_prev = level_mesh[cur_level-1].num_cells;
    }
  else
    {
      nverts_prev = _inverts.size();
      nents_prev = _incells.size();
    }

  std::vector<EntityHandle> inci_ent, child_ents;
  std::vector<int> inci_lid, child_lids;

  //Step 1: Update the V2HF maps for old/duplicate vertices
  for (int i=0; i<nverts_prev; i++)
    {
      inci_ent.clear(); inci_lid.clear(); child_ents.clear(); child_lids.clear();

      //Vertex id in the previous mesh
      EntityHandle vid;
      if (cur_level)
        vid = level_mesh[cur_level-1].start_vertex + i;
      else
        vid = _inverts[i];
      EntityHandle cur_vid = level_mesh[cur_level].start_vertex + i;

      //Get the incident half-vert in the previous mesh
      error = ahf->get_incident_map(type, vid, inci_ent, inci_lid);  MB_CHK_ERR(error);

      // Obtain the corresponding incident child in the current mesh
      for (int j=0;j<(int)inci_ent.size(); j++){
          int lvid = get_local_vid(vid, inci_ent[j], cur_level-1);
          int chid = refTemplates[type-1][d].v2hf[lvid][0]-1;

          int pid;
          if (cur_level)
            pid = inci_ent[j] - level_mesh[cur_level-1].start_cell;
          else
            pid = inci_ent[j] - *_incells.begin();

          int ind = nchilds*pid;

        //  EntityHandle child_ent = level_mesh[cur_level].start_cell + ind+chid ;
         // int child_lid = refTemplates[type-1][d].v2hf[lvid][1];
          child_ents.push_back(level_mesh[cur_level].start_cell + ind+chid);
          child_lids.push_back(refTemplates[type-1][d].v2hf[lvid][1]);
        }

      error = ahf->set_incident_map(type, cur_vid, child_ents, child_lids);  MB_CHK_ERR(error);
    }

//  error = ahf->determine_incident_halffaces( level_mesh[cur_level].cells);MB_CHK_ERR(error);

  //Step 2: Update SIBHFS maps
  for (int i=0; i< nents_prev; i++)
    {
      EntityHandle ent;
      if (cur_level)
        ent = level_mesh[cur_level-1].start_cell + i;
      else
        ent = _incells[i];

      EntityHandle *sib_entids = new EntityHandle[nhf];
      int *sib_lids = new int[nhf];

      error = ahf->get_sibling_map(type, ent, &sib_entids[0], &sib_lids[0], nhf);  MB_CHK_ERR(error);

      int id, idx;

      for (int l=0; l < nhf; l++)
        {

          if (!sib_entids[l])
            continue;

          //Get the number of children incident on this half-face
          int nch = refTemplates[type-1][d].ents_on_pent[l][0];

          //Get the order of children indices incident on this half-face
          int *id_sib = new int[nch];
          for (int k=0; k<nch; k++)
            id_sib[k] = 0;

          error = reorder_indices(cur_level, deg, ent, l, sib_entids[l], sib_lids[l], 1, id_sib);  MB_CHK_ERR(error);

          //Get the parent index of the sibling cell
          int psib;
          if (cur_level)
            psib = sib_entids[l] - level_mesh[cur_level-1].start_cell;
          else
            psib = sib_entids[l] - *_incells.begin();

          int plid = sib_lids[l];
          int  sidx = nchilds*psib;

          //Loop over all the childs incident on the working half-face
          idx = nchilds*i;

          for (int k=0; k<nch; k++)
            {
              id = refTemplates[type-1][d].ents_on_pent[l][k+1]-1;
              EntityHandle child_ent = level_mesh[cur_level].start_cell + idx + id ;
              int child_lid = l;

              //Find the sibling of the working child
              EntityHandle child_sibent;
              int child_siblid;
              error = ahf->get_sibling_map(type, child_ent, child_lid, child_sibent, child_siblid); MB_CHK_ERR(error);

              if (child_sibent != 0)
                continue;

              //Get the correponding child of the sibling of the current parent
              // We have already computed the order the children on incident corresponding to the working half-face
              id = refTemplates[type-1][d].ents_on_pent[plid][id_sib[k]]-1;

              EntityHandle psib_child = level_mesh[cur_level].start_cell + sidx+id ;
              int psib_chlid = plid;

              //Set the siblings of children incident on current half-face
              error = ahf->set_sibling_map(type, child_ent, child_lid, psib_child, psib_chlid);  MB_CHK_ERR(error);

              //Set the sibling of the sibling of the children to the children
              error = ahf->set_sibling_map(type, psib_child, psib_chlid, child_ent, child_lid); MB_CHK_ERR(error);
            }

          delete [] id_sib;
        }

      delete [] sib_entids;
      delete [] sib_lids;


      //Loop over edges to check if there are any non-manifold edges. If there are then the v2hfs map should be updated for the new vertices on it.
      const EntityHandle* conn;
      error = mbImpl->get_connectivity(ent, conn, nvpc);MB_CHK_ERR(error);

      for (int l=0; l<ne; l++)
        {
          id = ahf->lConnMap3D[index].e2v[l][0];
          EntityHandle v_start = conn[id];
          id = ahf->lConnMap3D[index].e2v[l][1];
          EntityHandle v_end = conn[id];

          std::vector<EntityHandle> inci_ent1, inci_ent2, inci_vent;
          std::vector<int> inci_lid1, inci_lid2, inci_vlid;

          error = ahf->get_incident_map(type, v_start, inci_ent1, inci_lid1);MB_CHK_ERR(error);
          error = ahf->get_incident_map(type, v_end, inci_ent2, inci_lid2);MB_CHK_ERR(error);

          if (inci_ent1.size() >1 && inci_ent2.size()>1)
            {
              std::vector<EntityHandle> cell_comps;
              std::vector<int> leid_comps, lfid_comps;

              error = ahf->get_half_facet_in_comp(ent, l , cell_comps, leid_comps, lfid_comps);MB_CHK_ERR(error);

              for (int j=0; j<(int)cell_comps.size(); j++)
                {
                  std::cout<<"cell_comps["<<j<<"] = "<<cell_comps[j]<<std::endl;
                  std::cout<<"leid_comps["<<j<<"] = "<<leid_comps[j]<<std::endl;
                  std::cout<<"lfid_comps["<<j<<"] = "<<lfid_comps[j]<<std::endl;
                }

              std::cout<<"ncomps = "<<cell_comps.size()<<std::endl;

              int nv =  refTemplates[type-1][d].nv_edge;
              std::vector<EntityHandle> set_child, everts;
              std::vector<int> set_chlid;
              std::vector<int> child_ids, child_lvids;

              //Get all the vertices on the local edge of the first component
              error = get_lid_inci_child(type, deg, lfid_comps[0], leid_comps[0], child_ids, child_lvids);



              int ind;
              if (cur_level)
                ind = level_mesh[cur_level].cells.index(cell_comps[0]);
              else
                ind = _incells.index(cell_comps[0]);

              for (int j=0; j<(int)child_ids.size(); j++)
                {
                  EntityHandle cent = level_mesh[cur_level].start_cell + ind*nchilds + child_ids[j];
                  const EntityHandle* econn;
                  error = mbImpl->get_connectivity(cent, econn, nvpc);MB_CHK_ERR(error);

                  everts.push_back(econn[child_lvids[j]]);
                }

              std::sort(everts.begin(), everts.end());
              std::vector<EntityHandle>::iterator last = std::unique(everts.begin(), everts.end());
              everts.erase(last, everts.end());

              std::cout<<"everts.size = "<<everts.size()<<std::endl;

              for (int k=0; k< nv; k++)
                {
                  set_child.clear(); set_chlid.clear();

                  EntityHandle vid = everts[k];

                  error = ahf->get_incident_map(type, vid, inci_vent, inci_vlid);MB_CHK_ERR(error);

                  if (inci_vent.size() > 1)
                    continue;

                  for (int c=0; c<(int)cell_comps.size(); c++){
                      child_ids.clear(); child_lvids.clear();
                      error = get_lid_inci_child(type, deg, lfid_comps[c], leid_comps[c], child_ids, child_lvids);MB_CHK_ERR(error);

                      if (cur_level)
                        ind = level_mesh[cur_level].cells.index(cell_comps[c]);
                      else
                        ind = _incells.index(cell_comps[c]);

                      for (int j=0; j<(int)child_ids.size(); j++)
                        {
                          std::cout<<"child_id = "<<child_ids[j]<<std::endl;
                          EntityHandle cent = level_mesh[cur_level].start_cell + ind*nchilds + child_ids[j];
                          std::cout<<"child_entid = "<<cent<<std::endl;
                          const EntityHandle* econn;
                          error = mbImpl->get_connectivity(cent, econn, nvpc);MB_CHK_ERR(error);

                          if (econn[child_lvids[j]] == vid)
                            {
                              set_child.push_back(cent);
                              set_chlid.push_back(lfid_comps[c]);
                              break;
                            }
                        }
                    }

                  std::cout<<"set_size = "<<set_child.size()<<std::endl;

                  error = ahf->set_incident_map(type, vid, set_child, set_chlid);MB_CHK_ERR(error);
                }
            }
        }
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::get_lid_inci_child(EntityType type, int deg, int lfid, int leid, std::vector<int> &child_ids, std::vector<int> &child_lvids)
{
  int index = ahf->get_index_in_lmap(*_incells.begin());
  int d = get_index_from_degree(deg);

 // int lv0 = ahf->lConnMap3D[index].e2v[leid][0];
//  int lv1 = ahf->lConnMap3D[index].e2v[leid][1];
  int nvpc = ahf->lConnMap3D[index].num_verts_in_cell;

  int nv =  refTemplates[type-1][d].nv_edge;
  int nch = refTemplates[type-1][d].ents_on_pent[lfid][0];

  for (int i=0; i< nch; i++)
    {
      int id = refTemplates[type-1][d].ents_on_pent[lfid][i+1]-1;
      for (int j=0; j< nvpc; j++)
        {
          int lv = refTemplates[type-1][d].ents_conn[id][j];
          for (int k=0; k<nv; k++)
            {
              if (lv == refTemplates[type-1][d].vert_on_edges[leid][k])
                {
                  child_ids.push_back(id);
                  child_lvids.push_back(j);
                }
            }
        }
    }

  return MB_SUCCESS;
}


ErrorCode NestedRefine::update_global_ahf_3D(int cur_level, int deg, std::vector<int> &pattern_ids)
{
  ErrorCode error;
  int nhf, nchilds, nverts_prev, nents_prev;

  EntityType type = MBTET;
  int index = ahf->get_index_in_lmap(*_incells.begin());
  int d = get_index_from_degree(deg);

  nhf = ahf->lConnMap3D[index].num_faces_in_cell;
  nchilds = refTemplates[type-1][d].total_new_ents;

  if (cur_level)
    {
      nverts_prev = level_mesh[cur_level-1].num_verts;
      nents_prev = level_mesh[cur_level-1].num_cells;
    }
  else
    {
      nverts_prev = _inverts.size();
      nents_prev = _incells.size();
    }

  std::vector<EntityHandle> inci_ent, child_ents;
  std::vector<int> inci_lid, child_lids;

  //Step 1: Update the V2HF maps for old/duplicate vertices
  for (int i=0; i<nverts_prev; i++)
    {
      inci_ent.clear(); inci_lid.clear(); child_ents.clear(); child_lids.clear();

      //Vertex id in the previous mesh
      EntityHandle vid;
      if (cur_level)
        vid = level_mesh[cur_level-1].start_vertex + i;
      else
        vid = _inverts[i];
      EntityHandle cur_vid = level_mesh[cur_level].start_vertex + i;

      //Get the incident half-vert in the previous mesh
      error = ahf->get_incident_map(type, vid, inci_ent, inci_lid);  MB_CHK_ERR(error);

      // Obtain the corresponding incident child in the current mesh
      for (int j=0; j<(int)inci_ent.size();j++){
          int lvid = get_local_vid(vid, inci_ent[j], cur_level-1);
          int chid = refTemplates[type-1][d].v2hf[lvid][0]-1;

          int pid;
          if (cur_level)
            pid = inci_ent[j] - level_mesh[cur_level-1].start_cell;
          else
            pid = inci_ent[j] - *_incells.begin();

          int ind = nchilds*pid;

          //EntityHandle child_ent = level_mesh[cur_level].start_cell + ind+chid ;
          //int child_lid = refTemplates[type-1][d].v2hf[lvid][1];
          child_ents.push_back(level_mesh[cur_level].start_cell + ind+chid);
          child_lids.push_back(refTemplates[type-1][d].v2hf[lvid][1]);
        }

      error = ahf->set_incident_map(type, cur_vid, child_ents, child_lids);  MB_CHK_ERR(error);
    }

 //   error = ahf->determine_incident_halffaces( level_mesh[cur_level].cells);MB_CHK_ERR(error);

  //Step 2: Update SIBHFS maps
  for (int i=0; i< nents_prev; i++)
    {
      EntityHandle ent;
      if (cur_level)
        ent = level_mesh[cur_level-1].start_cell + i;
      else
        ent = _incells[i];

      EntityHandle *sib_entids = new EntityHandle[nhf];
      int *sib_lids = new int[nhf];

      error = ahf->get_sibling_map(type, ent, &sib_entids[0], &sib_lids[0], nhf);  MB_CHK_ERR(error);

      int id, idx;

      for (int l=0; l < nhf; l++)
        {
          if (!sib_entids[l])
            continue;

          //Get the number of children incident on this half-face
          int pat_id = pattern_ids[i];
          int nch = refTemplates[pat_id][d].ents_on_pent[l][0];

          //Get the order of children indices incident on this half-face
          int *id_sib = new int[nch];
          for (int k=0; k<nch; k++)
            id_sib[k] = 0;

          error = reorder_indices(cur_level, deg, ent, l, sib_entids[l], sib_lids[l], 1, id_sib);  MB_CHK_ERR(error);

          //Get the parent index of the sibling cell
          int psib;
          if (cur_level)
            psib = sib_entids[l] - level_mesh[cur_level-1].start_cell;
          else
            psib = sib_entids[l] - *_incells.begin();

          int plid = sib_lids[l];
          int  sidx = nchilds*psib;
          int sibpat_id = pattern_ids[psib];

          //Loop over all the childs incident on the working half-face
          idx = nchilds*i;

          for (int k=0; k<nch; k++)
            {
              id = refTemplates[pat_id][d].ents_on_pent[l][k+1]-1;
              EntityHandle child_ent = level_mesh[cur_level].start_cell + idx + id ;
              int child_lid = l;

              EntityHandle child_sibent;
              int child_siblid;
              error = ahf->get_sibling_map(type, child_ent, child_lid, child_sibent, child_siblid); MB_CHK_ERR(error);

              if (child_sibent != 0)
                continue;

              //Get the correponding child of the sibling of the current parent
              // We have already computed the order the children on incident corresponding to the working half-face
              id = refTemplates[sibpat_id][d].ents_on_pent[plid][id_sib[k]]-1;

              EntityHandle psib_child = level_mesh[cur_level].start_cell + sidx+id ;
              int psib_chlid = plid;

              //Set the siblings of children incident on current half-face
              error = ahf->set_sibling_map(type, child_ent, child_lid, psib_child, psib_chlid);  MB_CHK_ERR(error);

              //Set the sibling of the sibling of the children to the children
              error = ahf->set_sibling_map(type, psib_child, psib_chlid, child_ent, child_lid); MB_CHK_ERR(error);
            }

          delete [] id_sib;
        }

      delete [] sib_entids;
      delete [] sib_lids;

    }

  return MB_SUCCESS;
}


/* **********************************
 *  *          Boundary Functions      *
 ************************************/

bool NestedRefine::is_vertex_on_boundary(const EntityHandle &vertex)
{
  ErrorCode error;
  EntityHandle sibents[27];
  int siblids[27];
  std::vector<EntityHandle> ent;
  std::vector<int> lid;

  int nhf;
  if (elementype == MBEDGE)
    nhf = 2;
  else if ((elementype == MBTRI)||(elementype == MBQUAD))
    nhf = ahf->lConnMap2D[elementype-2].num_verts_in_face;
  else if ((elementype == MBTET) || (elementype == MBHEX))
    {
      int idx = ahf->get_index_in_lmap(*_incells.begin());
      nhf = ahf->lConnMap3D[idx].num_faces_in_cell;
    }
  else
     MB_SET_ERR(MB_FAILURE, "Requesting vertex boundary information for an unsupported entity type");

  error = ahf->get_incident_map(elementype, vertex, ent, lid);MB_CHK_ERR(error);
  error = ahf->get_sibling_map(elementype, ent[0], &sibents[0], &siblids[0], nhf);MB_CHK_ERR(error);

  return (sibents[lid[0]] == 0);
}

bool NestedRefine::is_edge_on_boundary(const EntityHandle &entity)
{
  ErrorCode error;
  bool is_border = false;
  if (meshdim == 1) //The edge has a vertex on the boundary in the curve mesh
  {
    EntityHandle sibents[2];
    int siblids[2];
    error = ahf->get_sibling_map(MBEDGE, entity, &sibents[0], &siblids[0], 2);MB_CHK_ERR(error);
    for (int i = 0; i < 2; i++)
    {
      if (sibents[i] == 0)
      {
        is_border = true;
        break;
      }
    }
  }
  else if (meshdim == 2) //The edge is on the boundary of the 2d mesh
  {
    std::vector<EntityHandle> adjents;
    error = ahf->get_up_adjacencies_2d(entity, adjents);MB_CHK_ERR(error);
    if (adjents.size() == 1)
      is_border = true;
  }
  else if (meshdim == 3) //The edge is part of a face on the boundary of the 3d mesh
  {
    std::vector<EntityHandle> adjents;
    std::vector<int> leids;
    error = ahf->get_up_adjacencies_edg_3d(entity, adjents, &leids);MB_CHK_ERR(error);
    assert(!adjents.empty());

    int index = ahf->get_index_in_lmap(adjents[0]);
    int nhf = ahf->lConnMap3D[index].num_faces_in_cell;

    for (int i = 0; i < (int)adjents.size(); i++)
    {
      EntityHandle sibents[6];
      int siblids[6];
      error = ahf->get_sibling_map(elementype, adjents[0], &sibents[0], &siblids[0], nhf);MB_CHK_ERR(error);
      for (int k = 0; k < 2; k++)
      {
        int hf = ahf->lConnMap3D[index].e2hf[leids[0]][k];
        if (sibents[hf] == 0)
        {
          is_border = true;
          break;
        }
      }
    }
  }
  return is_border;
}

bool NestedRefine::is_face_on_boundary(const EntityHandle &entity)
{
  ErrorCode error;
  bool is_border = false;

  if (meshdim == 1)
    MB_SET_ERR(MB_FAILURE, "Requesting boundary information for a face entity type on a curve mesh");
  else if (meshdim == 2) //The face has a local edge on the boundary of the 2d mesh
  {
    EntityHandle sibents[4];
    int siblids[4];
    int nepf = ahf->lConnMap2D[elementype - 2].num_verts_in_face;
    error = ahf->get_sibling_map(elementype, entity, &sibents[0], &siblids[0], nepf);MB_CHK_ERR(error);

    for (int i = 0; i < nepf; i++)
    {
      if (sibents[i] == 0)
      {
        is_border = true;
        break;
      }
    }
  }
  else if (meshdim == 3)//The face lies on the boundary of the 3d mesh
  {
    std::vector<EntityHandle> adjents;
    error = ahf->get_up_adjacencies_face_3d(entity, adjents);MB_CHK_ERR(error);
    if (adjents.size() == 1)
      is_border = true;
  }
  return is_border;
}

bool NestedRefine::is_cell_on_boundary(const EntityHandle &entity)
{
  if (meshdim != 3)
    MB_SET_ERR(MB_FAILURE, "Requesting boundary information for a cell entity type on a curve or surface mesh");

  bool is_border = false;
  int index = ahf->get_index_in_lmap(*_incells.begin());
  int nfpc = ahf->lConnMap3D[index].num_faces_in_cell;
  EntityHandle sibents[6];
  int siblids[6];

  ErrorCode error = ahf->get_sibling_map(elementype, entity, &sibents[0], &siblids[0], nfpc);MB_CHK_ERR(error);

  for (int i = 0; i < nfpc; i++)
  {
    if (sibents[i] == 0)
    {
      is_border = true;
      break;
    }
  }
  return is_border;
}

/* **********************************
 *          Helper Functions              *
************************************/
ErrorCode NestedRefine::copy_vertices_from_prev_level(int cur_level)
{
  ErrorCode error;

  if (cur_level)
    {
      int nverts_prev = level_mesh[cur_level-1].num_verts;
      for (int i = 0; i < nverts_prev; i++)
        {
          level_mesh[cur_level].coordinates[0][i] = level_mesh[cur_level-1].coordinates[0][i];
          level_mesh[cur_level].coordinates[1][i] = level_mesh[cur_level-1].coordinates[1][i];
          level_mesh[cur_level].coordinates[2][i] = level_mesh[cur_level-1].coordinates[2][i] ;
        }
    }
  else // Copy the vertices from the input mesh
    {
      int nverts_in = _inverts.size();
      double *vcoords = new double[3*nverts_in];
      error = mbImpl->get_coords(_inverts, vcoords); MB_CHK_ERR(error);

      for (int i = 0; i < nverts_in; i++)
        {
          level_mesh[cur_level].coordinates[0][i] = vcoords[3*i];
          level_mesh[cur_level].coordinates[1][i] = vcoords[3*i+1];
          level_mesh[cur_level].coordinates[2][i] = vcoords[3*i+2];
        }

      delete [] vcoords;

    }
  return MB_SUCCESS;
  //To add: Map from old vertices to new duplicates: NOT NEEDED
}

ErrorCode NestedRefine::update_tracking_verts(EntityHandle cid, int cur_level, int deg, std::vector<EntityHandle> &trackvertsC_edg, std::vector<EntityHandle> &trackvertsC_face, EntityHandle *vbuffer)
{
  //The vertices in the vbuffer are added to appropriate edges and faces of cells that are incident on the working cell.
  ErrorCode error;

  EntityHandle cstart_prev;
  if (cur_level)
    cstart_prev = level_mesh[cur_level-1].start_cell;
  else
    cstart_prev = *_incells.begin();

  EntityType cell_type = mbImpl->type_from_handle(cstart_prev);
  int cindex = cell_type -1;
  int d = get_index_from_degree(deg);

  int nve = refTemplates[cindex][d].nv_edge;
  int nvf = refTemplates[cindex][d].nv_face;

  int index = ahf->get_index_in_lmap(*(_incells.begin()));
  int nepc = ahf->lConnMap3D[index].num_edges_in_cell;
  int nfpc = ahf->lConnMap3D[index].num_faces_in_cell;

  //Step 1: Add the vertices on an edge of the working cell to tracking array of incident cells.
  for (int i=0; i<nepc; i++)
    {
      //Add the vertices to edges of the current cell
      for (int j=0; j< nve; j++)
        {
          int id = refTemplates[cindex][d].vert_on_edges[i][j];
          int idx = cid - cstart_prev;
          int aid = idx*nve*nepc+nve*i+j;

          if (!trackvertsC_edg[aid])
            trackvertsC_edg[aid] = vbuffer[id];
        }

      //Obtain all the incident cells
      std::vector<EntityHandle> inc_cids;
      std::vector<int> inc_leids, inc_orient;

      error = ahf->get_up_adjacencies_edg_3d(cid, i, inc_cids,  &inc_leids, &inc_orient);   MB_CHK_ERR(error);

      if (inc_cids.size() ==1)
        continue;

      //Add the vertices to the edges of the incident cells
      for (int k=0; k<(int)inc_cids.size(); k++)
        {
          if (inc_cids[k] == cid)
            continue;

          int idx = inc_cids[k] - cstart_prev;

          if (inc_orient[k]) // Same edge direction as the current edge
            {
              for (int j=0; j< nve; j++)
                {
                  int  id = refTemplates[cindex][d].vert_on_edges[i][j];
                  int aid = idx*nve*nepc+nve*inc_leids[k]+j;

                  if (!trackvertsC_edg[aid])
                    trackvertsC_edg[aid] = vbuffer[id];
                }
            }
          else
            {
              for (int j=0; j< nve; j++)
                {
                  int id = refTemplates[cindex][d].vert_on_edges[i][nve-j-1];
                  int aid = idx*nve*nepc+nve*inc_leids[k]+j;

                  if (!trackvertsC_edg[aid])
                    trackvertsC_edg[aid] = vbuffer[id];
                }
            }
        }
    }

  //Step 2: Add the vertices on a face of the working cell to tracking array of incident cells.
  if (nvf) {

      for (int i=0; i< nfpc; i++)
        {
          //Add vertices to the tracking array of vertices on faces for the current cell
          std::vector<EntityHandle> face_vbuf(nvf,0);
          for (int j=0; j< nvf; j++)
            {
              int  id = refTemplates[cindex][d].vert_on_faces[i][j];
              int idx = cid - cstart_prev;
              int aid = idx*nvf*nfpc+nvf*i+j;

              if (!trackvertsC_face[aid])
                trackvertsC_face[aid] = vbuffer[id];

              face_vbuf[j] = vbuffer[id];
            }

          //Obtain all the incident cells
          std::vector<EntityHandle> sib_cids;
          std::vector<int>  sib_lfids;
          error = ahf->get_up_adjacencies_face_3d(cid, i, sib_cids, &sib_lfids);   MB_CHK_ERR(error);

          if (sib_cids.size()==1)
            continue;

          //Reorder the vertex local ids incident on the half-face
          int *id_sib = new int[nvf];
          for (int k=0; k<nvf; k++)
            id_sib[k] = 0;

          error = reorder_indices(cur_level, deg, sib_cids[1], sib_lfids[1], cid, i, 0, id_sib);  MB_CHK_ERR(error);

          //Add vertices to the tracking array of vertices on faces for the sibling cell of the current cell
          for (int j=0; j< nvf; j++)
            {
              int idx = sib_cids[1] - cstart_prev;
              int aid = idx*nvf*nfpc+nvf*sib_lfids[1]+j;

              if (!trackvertsC_face[aid])
                trackvertsC_face[aid] = face_vbuf[id_sib[j]-1];
            }

          delete [] id_sib;
        }
    }
  return MB_SUCCESS;
}

ErrorCode NestedRefine::reorder_indices(int cur_level, int deg, EntityHandle cell, int lfid, EntityHandle sib_cell, int sib_lfid, int index, int *id_sib)
{
  // Reorders the indices of either vertices or children cell local ids to match with order of the given cell and a local face.
  // index = 0 : vertices,
  //           = 1 : face


  assert(deg ==2 || deg == 3);

  ErrorCode error;
  int idx = ahf->get_index_in_lmap(*_incells.begin());
  int nvF = ahf->lConnMap3D[idx].hf2v_num[lfid];
  int nco = permutation[nvF-3].num_comb;

  if (!index && ((nvF==3 && deg==3)||(nvF==4 && deg==2)))
    {
      id_sib[0] = 1;
    }
  else
    {
      //Get connectivity of the cell and its sibling cell
      std::vector<EntityHandle> conn, sib_conn;
      error = get_connectivity(cell, cur_level, conn);   MB_CHK_ERR(error);

      error = get_connectivity(sib_cell, cur_level, sib_conn);   MB_CHK_ERR(error);

      //Get the connectivity of the local face in the cell and its sibling
      EntityHandle *lface = new EntityHandle[nvF];
      EntityHandle *lface_sib = new EntityHandle[nvF];
      for (int i=0; i<nvF; i++)
        {
          int id = ahf->lConnMap3D[idx].hf2v[lfid][i];
          lface[i] = conn[id];

          id = ahf->lConnMap3D[idx].hf2v[sib_lfid][i];
          lface_sib[i] = sib_conn[id];
        }

      //Find the combination
      int c = 0;
      for (int i=0; i<nco; i++)
        {
          int count = 0;
          for (int j=0; j<nvF; j++)
            {
              int id = permutation[nvF-3].comb[i][j];
              if (lface[j] == lface_sib[id])
                count += 1;
            }

          if (count == nvF)
            {
              c=i;
              break;
            }
        }

      if (c>nco)
        MB_SET_ERR(MB_FAILURE, "Getting a combination number more than currently supported");

      //Get the ordered indices
      if (((!index)&&(nvF==4)&&(deg==3)) || (deg==2))
        {
          for (int i=0; i<4; i++)
            id_sib[i] = permutation[nvF-3].porder2[c][i];
        }
        else
        {
          for (int i=0; i<9; i++)
            id_sib[i] = permutation[nvF-3].porder3[c][i];
        }

      delete [] lface;
      delete [] lface_sib;
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::reorder_indices(int deg, EntityHandle *face1_conn, EntityHandle *face2_conn, int nvF, std::vector<int> &lemap, std::vector<int> &vidx, int *leorient)
{
  //Given the connectivities of two faces, get the permuted indices w.r.t first face.
  //Step 1: First find the orientation
  int nco = permutation[nvF-3].num_comb;
  int c = 0;
  for (int i=0; i<nco; i++)
    {
      int count = 0;
      for (int j=0; j<nvF; j++)
        {
          int id = permutation[nvF-3].comb[i][j];
          if (face1_conn[j] == face2_conn[id])
            count += 1;
        }

      if (count == nvF)
        {
          c=i;
          break;
        }
    }

  if (c>nco)
    MB_SET_ERR(MB_FAILURE, "Getting a combination number more than currently supported");

  //Add the corresponding local edges
  lemap.reserve(nvF);
  for (int i=0; i<nvF; i++)
    {
      lemap.push_back(permutation[nvF-3].lemap[c][i]);
    }
  if (leorient)
    leorient[0] = permutation[nvF-3].orient[c];

  if (nvF==3&&deg==2)
    return MB_SUCCESS;

  if ((nvF==3 && deg==3)||(nvF==4 && deg==2))
    {
      vidx.push_back(1);
    }
  else if (nvF==4 && deg==3)
    {
      for (int i=0; i<4; i++)
        vidx.push_back(permutation[nvF-3].porder2[c][i]);
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::count_subentities(EntityHandle set, int cur_level, int *nedges, int *nfaces)
{
  ErrorCode error;

  if (cur_level>=0){
      Range edges, faces, cells;

      error = mbImpl->get_entities_by_dimension( set, 1, edges); MB_CHK_ERR(error);

      error = mbImpl->get_entities_by_dimension( set, 2, faces); MB_CHK_ERR(error);

      error = mbImpl->get_entities_by_dimension( set, 3, cells); MB_CHK_ERR(error);

      error = ahf->count_subentities(edges, faces, cells, nedges, nfaces); MB_CHK_ERR(error);
    }
  else
    {
      error = ahf->count_subentities(_inedges, _infaces, _incells, nedges, nfaces); MB_CHK_ERR(error);
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::get_octahedron_corner_coords(int cur_level, int deg, EntityHandle *vbuffer, double *ocoords)
{
  int lid[6]={0,0,0,0,0,0};

  if (deg==2)
    {
      lid[0] = 5; lid[1] = 8; lid[2] = 9;
      lid[3] = 6; lid[4] = 4; lid[5] = 7;
    }
  else if (deg ==3)
    {
      lid[0] = 19; lid[1] = 16; lid[2] = 18;
      lid[3] = 9; lid[4] = 4; lid[5] = 10;
    }

  EntityHandle vstart = level_mesh[cur_level].start_vertex;

  for (int i=0; i<6; i++)
    {
      EntityHandle vid = vbuffer[lid[i]];
      ocoords[3*i] = level_mesh[cur_level].coordinates[0][vid-vstart];
      ocoords[3*i+1] =  level_mesh[cur_level].coordinates[1][vid-vstart];
      ocoords[3*i+2] =  level_mesh[cur_level].coordinates[2][vid-vstart];
    }

  return MB_SUCCESS;
}

int NestedRefine::find_shortest_diagonal_octahedron(int cur_level, int deg, EntityHandle *vbuffer)
{
  ErrorCode error;
  double coords[18];
  error = get_octahedron_corner_coords(cur_level, deg, vbuffer, coords);
  if (error != MB_SUCCESS)
    MB_SET_ERR(MB_FAILURE,"Error in obtaining octahedron corner coordinates");

  int diag_map[6] = {1,3,2,4,5,0};
  double length = std::numeric_limits<double>::max();

  int diag = 0;
  double x, y,z;
  x = y = z = 0;

  for (int d=0; d<3; d++)
    {
      int id1 = diag_map[2*d];
      int id2 = diag_map[2*d+1];
      x = coords[3*id1] - coords[3*id2];
      y = coords[3*id1+1] - coords[3*id2+1];
      z = coords[3*id1+2] - coords[3*id2+2];
      double  dlen =  sqrt(x*x+y*y+z*z);
      if ( dlen < length)
        {
          length = dlen;
          diag = d+1;
        }
    }

  return diag;
}

int NestedRefine::get_local_vid(EntityHandle vid, EntityHandle ent, int level)
{
  ErrorCode error;
  //Given a vertex, find its local id in the given entity
  std::vector<EntityHandle> conn;

  error = get_connectivity(ent, level+1, conn);
  if (error != MB_SUCCESS)
    MB_SET_ERR(MB_FAILURE, "Error in getting connectivity of the requested entity");

  int lid=-1;
  for (int i=0; i<(int)conn.size(); i++)
    {
      if (conn[i] == vid)
        {
          lid = i;
          break;
        }
    }
  return lid;
}

int NestedRefine::get_index_from_degree(int degree)
{
  int d = deg_index.find(degree)->second;
  return d;
}

/*
ErrorCode NestedRefine::print_maps_1D(int level)
{
  ErrorCode error;
  int nv, ne;
  nv = level_mesh[level].num_verts;
  ne = level_mesh[level].num_edges;

  EntityHandle start_edge = level_mesh[level].start_edge;

  //V2HV
    std::cout<<"<V2HV_EID, V2HV_LVID>"<<std::endl;
  for (int i=0; i<nv; i++)
    {
      EntityHandle eid=0;
      int lvid=0;
      EntityHandle vid = level_mesh[level].start_vertex+i;
      error = ahf->get_incident_map(MBEDGE, vid, eid, lvid); MB_CHK_ERR(error);

      std::cout<<"For vertex = "<<vid<<"::Incident halfvertex "<<eid<<"  "<<lvid<<std::endl;
    }

  //SIBHVS
  std::cout<<"start_edge = "<<start_edge<<std::endl;
  std::cout<<"<SIBHVS_EID,SIBHVS_LVID>"<<std::endl;
  for (int i=0; i<ne; i++)
    {
      EntityHandle ent = start_edge+i;

      EntityHandle eid[2];  int lvid[2];
      error = ahf->get_sibling_map(MBEDGE, ent, &eid[0], &lvid[0], 2); MB_CHK_ERR(error);
      std::cout<<"<"<<eid[0]<<","<<lvid[0]<<">"<<"      "<<"<"<<eid[1]<<","<<lvid[1]<<">"<<std::endl;
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::print_maps_2D(int level, EntityType type)
{
  ErrorCode error;
  int nv, nf;
  nv = level_mesh[level].num_verts;
  nf = level_mesh[level].num_faces;

  EntityHandle start_face = level_mesh[level].start_face;

  //V2HV
    std::cout<<"<V2HE_FID, V2HE_LEID>"<<std::endl;
  for (int i=0; i<nv; i++)
    {
      EntityHandle fid=0;
      int leid=0;
      EntityHandle vid = level_mesh[level].start_vertex+i;
      error = ahf->get_incident_map(type, vid, fid, leid); MB_CHK_ERR(error);

      std::cout<<"For vertex = "<<vid<<"::Incident halfedge "<<fid<<"  "<<leid<<std::endl;
    }

  //SIBHES
  std::cout<<"start_face = "<<start_face<<std::endl;
  std::cout<<"<SIBHES_FID,SIBHES_LEID>"<<std::endl;
  EntityType ftype = mbImpl->type_from_handle(*_infaces.begin());
  int nepf = ahf->lConnMap2D[ftype-2].num_verts_in_face;

  EntityHandle *fid = new EntityHandle[nepf];
  int *leid = new int[nepf];

  for (int i=0; i<nf; i++)
    {
      for (int j=0; j<nepf; j++)
        {
          fid[j] = 0;
          leid[j] = 0;
        }

      EntityHandle ent = start_face+i;
      error = ahf->get_sibling_map(type, ent, fid, leid, nepf); MB_CHK_ERR(error);

      for (int j=0; j<nepf; j++){
          std::cout<<"<"<<fid[j]<<","<<leid[j]<<">"<<"      ";
      }
      std::cout<<std::endl;
    }

  delete [] fid;
  delete [] leid;

  return MB_SUCCESS;
}

ErrorCode NestedRefine::print_maps_3D(int level, EntityType type)
{
  ErrorCode error;
  int nv, nc;
  nv = level_mesh[level].num_verts;
  nc = level_mesh[level].num_cells;
  EntityHandle start_cell = level_mesh[level].start_cell;

  //V2HF
    std::cout<<"<V2HF_CID, V2HF_LFID>"<<std::endl;
  for (int i=0; i<nv; i++)
    {
      EntityHandle cid=0;
      int lfid=0;
      EntityHandle vid = level_mesh[level].start_vertex+i;
      error = ahf->get_incident_map(type, vid, cid, lfid); MB_CHK_ERR(error);

      std::cout<<"For vertex = "<<vid<<"::Incident halfface "<<cid<<"  "<<lfid<<std::endl;
    }

  //SIBHFS
  std::cout<<"start_cell = "<<start_cell<<std::endl;
  std::cout<<"<SIBHFS_CID,SIBHFS_LFID>"<<std::endl;
  int index = ahf->get_index_in_lmap(start_cell);
  int nfpc = ahf->lConnMap3D[index].num_faces_in_cell;

  EntityHandle *cid = new EntityHandle[nfpc];
  int *lfid = new int[nfpc];
  for (int i=0; i<nc; i++)
    {
      for (int k=0; k<nfpc; k++)
        {
          cid[k] = 0;
          lfid[k] = 0;
        }

      EntityHandle ent = start_cell+i;
      error = ahf->get_sibling_map(type, ent, cid, lfid, nfpc); MB_CHK_ERR(error);

      for (int j=0; j<nfpc; j++){
          std::cout<<"<"<<cid[j]<<","<<lfid[j]<<">"<<"      ";
      }
      std::cout<<std::endl;
    }

  delete [] cid;
  delete [] lfid;

  return MB_SUCCESS;
}

*/
}//namesapce moab

