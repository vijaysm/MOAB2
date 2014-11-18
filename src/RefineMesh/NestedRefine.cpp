#ifdef WIN32
#pragma warning (disable : 4786)
#endif 

#include "moab/NestedRefine.hpp"
#include "moab/HalfFacetRep.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Templates.hpp"
#include <iostream>
#include <assert.h>
#include <vector>
#include <limits>
#include <cmath>


namespace moab{

  NestedRefine::NestedRefine(Core *impl)
  {
    assert(NULL != impl);
    mbImpl = impl;

    ahf = new HalfFacetRep(mbImpl);
    assert(NULL != ahf);

    ErrorCode error;
    error = ahf->initialize();
    assert(error==MB_SUCCESS);
    error = ahf->get_entity_ranges(_inverts, _inedges, _infaces, _incells);
    assert(error==MB_SUCCESS);
  }

  NestedRefine::~NestedRefine()
  {
    delete ahf;
  }

  /************************************************************
   *     Interface Functions                                  *
   ************************************************************/

  ErrorCode NestedRefine::generate_mesh_hierarchy(int *level_degrees, int num_level, EntityHandle *hm_set)
  {
    ErrorCode error;

    error = generate_hm(level_degrees, num_level, hm_set);
    if (error != MB_SUCCESS) return error;

 //   delete [] hmest;
    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::get_connectivity(EntityHandle ent, int level, std::vector<EntityHandle> &conn)
  {
    ErrorCode error;
    EntityType type = mbImpl->type_from_handle(ent);
    EntityHandle start_ent ;
    if (level>= 0)
      {
        if (type == MBEDGE)
          {
            conn.reserve(2);
            start_ent = level_mesh[level].start_edge;
            conn.push_back( level_mesh[level].edge_conn[2*(ent-start_ent)]);
            conn.push_back( level_mesh[level].edge_conn[2*(ent-start_ent)+1]);
          }
        else if (type == MBTRI || type == MBQUAD)
          {
            int num_corners = ahf->local_maps_2d(ent);
            start_ent = level_mesh[level].start_face;
            for (int i=0; i<num_corners; i++)
              conn.push_back(level_mesh[level].face_conn[num_corners*(ent-start_ent)+i]);
          }
        else if (type == MBTET || type == MBPRISM || type == MBHEX)
          {
            int index = ahf->get_index_from_type(*_incells.begin());
            int num_corners = ahf->lConnMap3D[index].num_verts_in_cell;
            start_ent = level_mesh[level].start_cell;
            for (int i=0; i<num_corners; i++)
              conn.push_back(level_mesh[level].cell_conn[num_corners*(ent-start_ent)+i]);
          }
        else
          {
            std::cout<<"Unsupport element type"<<std::endl;
            return MB_FAILURE;
          }
      }
      else
      {
        error = mbImpl->get_connectivity(&ent, 1, conn);
        if (error != MB_SUCCESS) return error;
      }

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::get_coordinates(std::vector<EntityHandle> &verts, int num_verts, int cur_level, double *coords)
  {
    EntityHandle vstart = level_mesh[cur_level].start_vertex;
    for (int i=0; i< num_verts; i++)
      {
        EntityHandle vid = verts[i];
        coords[3*i] = level_mesh[cur_level].coordinates[0][vid-vstart];
        coords[3*i+1] =  level_mesh[cur_level].coordinates[1][vid-vstart];
        coords[3*i+2] =  level_mesh[cur_level].coordinates[2][vid-vstart];

        //std::cout<<" x = "<<coords[3*i]<<", y = "<<coords[3*i+1]<<", z = "<<coords[3*i+2]<<std::endl;
      }
    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::get_adjacencies(const EntityHandle source_entity,
                            const unsigned int target_dimension,
                            std::vector<EntityHandle> &target_entities)

  {
    ErrorCode error;
    error = ahf->get_adjacencies(source_entity, target_dimension, target_entities);
    if (error != MB_SUCCESS) return error;

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
    error = count_subentities(set, cur_level-1, &nedges, &nfaces);
    if (error != MB_SUCCESS) return error;

    int d = get_index_from_degree(level_degree);
    int nverts = refTemplates[MBEDGE-1][d].nv_edge*nedges;
    hmest[0] = nverts_prev + nverts;
    hmest[1] = nedges_prev*refTemplates[MBEDGE-1][d].total_new_ents;
    hmest[2] = 0;
    hmest[3] = 0;

    int findex, cindex;
    if (_infaces.size())
      {
        findex = mbImpl->type_from_handle(*(_infaces.begin()))-1;
        hmest[0] += refTemplates[findex][d].nv_face*nfaces_prev;
        hmest[2] = nfaces_prev*refTemplates[findex][d].total_new_ents;
      }
    if (_incells.size())
      {
        cindex = mbImpl->type_from_handle(*(_incells.begin()))-1;
        hmest[0] += refTemplates[cindex][d].nv_face*nfaces;
        hmest[0] += refTemplates[cindex][d].nv_cell*ncells_prev;
        hmest[3] = ncells_prev*refTemplates[cindex][d].total_new_ents;
      }

    //debug
      std::cout<<"vtotal ="<<hmest[0]<<", etotal = "<<hmest[1]<<", ftotal = "<<hmest[2]<<", ctotal = "<<hmest[3]<<std::endl;

  /*  int nverts_prev, nedges_prev, nfaces_prev, ncells_prev;
    for (int l=1; l<num_level; l++)
      {
        nverts_prev = hmest[4*(l-1)];
        nedges_prev = hmest[4*(l-1)+1];
        nfaces_prev = hmest[4*(l-1)+2];
        ncells_prev = hmest[4*(l-1)+3];

        d = get_index_from_degree(level_degrees[l]);
        nverts = refTemplates[MBEDGE-1][d].nv_edge*nedges_prev;
        hmest[4*l] = nverts_prev + nverts;
        hmest[4*l+1] = nedges_prev*refTemplates[MBEDGE-1][d].total_new_ents;
        hmest[4*l+2] = 0;
        hmest[4*l+3] = 0;

        if (_infaces.size())
          {
            hmest[4*l] += refTemplates[findex][d].nv_face*nfaces_prev;
            hmest[4*l+2] = nfaces_prev*refTemplates[findex][d].total_new_ents;
          }
        if (_incells.size())
          {
            hmest[4*l] +=  refTemplates[cindex][d].nv_cell*ncells_prev;
            hmest[4*l+3] = ncells_prev*refTemplates[cindex][d].total_new_ents;
          }
      }*/

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::create_hm_storage_single_level(EntityHandle *set, int cur_level, int estL[4])
  {
    //Obtain chunks of memory for the current level. Add them to a particular meshset.
    EntityHandle set_handle;
    ErrorCode error = mbImpl->create_meshset(MESHSET_SET, set_handle);
    if (error != MB_SUCCESS) return error;
    *set = set_handle;

    std::cout<<"Created mesh set for level["<<cur_level<<"] = "<<*set<<std::endl;

    ReadUtilIface *read_iface;
    error = mbImpl->query_interface(read_iface);
    if (error != MB_SUCCESS) return error;

    //Vertices
    error = read_iface->get_node_coords(3, estL[0], 0,  level_mesh[cur_level].start_vertex , level_mesh[cur_level].coordinates);
    if (error != MB_SUCCESS) return error;
    level_mesh[cur_level].num_verts = estL[0];

    EntityHandle *varray = new EntityHandle[estL[0]];
    for (int i=0; i< estL[0]; i++)
      varray[i] = level_mesh[cur_level].start_vertex + i;

    error = mbImpl->add_entities(*set, varray, estL[0]);
    if (error != MB_SUCCESS) return error;

    delete [] varray;

    // Edges
    if (estL[1])
      {
        error = read_iface->get_element_connect(estL[1], 2, MBEDGE, 0, level_mesh[cur_level].start_edge, level_mesh[cur_level].edge_conn);
        if (error != MB_SUCCESS) return error; 

        level_mesh[cur_level].num_edges = estL[1];

        EntityHandle *earray = new EntityHandle[estL[1]];
        for (int i=0; i< estL[1]; i++)
          earray[i] = level_mesh[cur_level].start_edge + i;

        error = mbImpl->add_entities(*set, earray, estL[1]);
        if (error != MB_SUCCESS) return error;

        delete [] earray;
      }
    else
      level_mesh[cur_level].num_edges = 0;

    //Faces
    if (estL[2])
      {
        EntityType type = mbImpl->type_from_handle(*(_infaces.begin()));
        int nvpf = ahf->local_maps_2d(*(_infaces.begin()));
        error = read_iface->get_element_connect(estL[2], nvpf, type, 0, level_mesh[cur_level].start_face, level_mesh[cur_level].face_conn);
        if (error != MB_SUCCESS) return error;
        level_mesh[cur_level].num_faces = estL[2];

        EntityHandle *farray = new EntityHandle[estL[2]];
        for (int i=0; i< estL[2]; i++)
          farray[i] = level_mesh[cur_level].start_face + i;

        error = mbImpl->add_entities(*set, farray, estL[2]);
        if (error != MB_SUCCESS) return error;

        delete [] farray;
      }
    else
      level_mesh[cur_level].num_faces = 0;

    //Cells
    if (estL[3])
      {
        EntityType type = mbImpl->type_from_handle(*(_incells.begin()));
        int index = ahf->get_index_from_type(*_incells.begin());
        int nvpc = ahf->lConnMap3D[index].num_verts_in_cell;
        error = read_iface->get_element_connect(estL[3], nvpc, type, 0, level_mesh[cur_level].start_cell, level_mesh[cur_level].cell_conn);
        if (error != MB_SUCCESS) return error;
        level_mesh[cur_level].num_cells = estL[3];

        EntityHandle *carray = new EntityHandle[estL[3]];
        for (int i=0; i< estL[3]; i++)
          carray[i] = level_mesh[cur_level].start_cell + i;

        error = mbImpl->add_entities(*set, carray, estL[3]);
        if (error != MB_SUCCESS) return error;

        delete [] carray;
      }
    else
      level_mesh[cur_level].num_cells = 0;

    return MB_SUCCESS;
  }



  /**********************************
   *   Hierarchical Mesh Generation  *
   * *********************************/

  ErrorCode NestedRefine::generate_hm(int *level_degrees, int num_level, EntityHandle *hm_set)
  {
    ErrorCode error;

    for (int l = 0; l<num_level; l++)
      {
        // Estimate storage
        int hmest[4] = {0,0,0,0};
        EntityHandle set;
        if (l)
          set = hm_set[l-1];
        else
          set = 0;
        error = estimate_hm_storage(set, level_degrees[l], l, hmest);
        if (error != MB_SUCCESS) return error;
        std::cout<<"Estimated storage successfully"<<std::endl;

        //Create arrays for storing the current level
        error = create_hm_storage_single_level(&hm_set[l], l, hmest);
        if (error != MB_SUCCESS) return error;

        std::cout<<"Created storage successfully"<<std::endl;

        //Copy the old vertices along with their coordinates
        error = copy_vertices_from_prev_level(l);
        if (error != MB_SUCCESS) return error;

        std::cout<<"Copied vertices from previous level successfully"<<std::endl;

        //Create the new entities and new vertices
        error = construct_hm_entities(l, level_degrees[l]);
        if (error != MB_SUCCESS) return error;

        std::cout<<"Created current level mesh successfully"<<std::endl;
      }
    return MB_SUCCESS;
  }


  ErrorCode NestedRefine::construct_hm_entities(int cur_level, int deg)
  {
    ErrorCode error;

    //Generate mesh for current level by refining previous level.
    if (ahf->thismeshtype == CURVE)
      {
        error = construct_hm_1D(cur_level, deg);
        if (MB_SUCCESS != error) return error;
      }
     else if(ahf->thismeshtype == SURFACE || ahf->thismeshtype == SURFACE_MIXED)
      {
        error = construct_hm_2D(cur_level, deg);
        if (MB_SUCCESS != error) return error;
      }
    else
      {
        error = construct_hm_3D(cur_level, deg);
        if (MB_SUCCESS != error) return error;
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

    //std::vector<EntityHandle> conn(2);
     std::vector<EntityHandle> conn;
    int count_nents = 0;
    int count_verts = nverts_prev;

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int eid = 0; eid < nents_prev; eid++)
      {
        conn.clear();
        //vbuffer.clear();

        // EntityHandle of the working edge
        EntityHandle edge;
        if (cur_level)
          edge = level_mesh[cur_level-1].start_edge + eid; // Makes the assumption initial mesh in contiguous in memory
        else
          edge = *_inedges.begin()+eid; // Makes the assumption initial mesh in contiguous in memory

        error = get_connectivity(edge, cur_level-1, conn);
        if (error != MB_SUCCESS) return error;

        //Add the vertex handles to vbuffer for the current level for the working edge

        // Since the old vertices are copied first, their local indices do not change as new levels are added.
        // Clearly the local indices of the new vertices introduced in the current level is still they the same
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

        error = update_local_ahf(deg, MBEDGE,  vbuffer, ent_buffer, etotal);
        if (error != MB_SUCCESS) return error;

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

    error = update_global_ahf(MBEDGE, cur_level, deg);
    if (error != MB_SUCCESS) return error;

 //   error = print_tags_1D(cur_level);
//    if (error != MB_SUCCESS) return error;

    delete [] vbuffer;

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
    int nepf = ahf->local_maps_2d(*_infaces.begin());
    EntityType type = mbImpl->type_from_handle(*(_infaces.begin()));
    int findex = type-1;

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
          face = *_infaces.begin() + fid;

  //      std::cout<<"Working face = "<<face<<std::endl;

        error = get_connectivity(face, cur_level-1, conn);
        if (error != MB_SUCCESS) return error;

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

    //    for (int i=0; i<etotal; i++)
    //      std::cout<<"ent_buffer["<<i<<"] = "<<ent_buffer[i]<<std::endl;


        // Step 3: Update the local AHF maps
        error = update_local_ahf(deg, type, vbuffer, ent_buffer, etotal);
        if (error != MB_SUCCESS) return error;

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
            error = ahf->get_up_adjacencies_2d(face, i, false, sibfids, &sibleids, &siborient);
            if (error != MB_SUCCESS) return error;

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

   //     for (int i=0; i< nents_prev*nepf*nve; i++)
   //       std::cout<<"trackF["<<i<<"] = "<<trackvertsF[i]<<std::endl;

        //Step 5: Compute the coordinates of the new vertices, avoids computing more than once via the flag_verts array.
        double *corner_coords = new double[nepf*3];
        error = get_coordinates(cur_conn, nepf, cur_level, corner_coords);
        if (error != MB_SUCCESS) return error;

        error = compute_coordinates(cur_level, deg, type, vbuffer, vtotal, corner_coords, flag_verts, nverts_prev);
        if (error != MB_SUCCESS) return error;

        delete [] corner_coords;

     //   error = print_tags_2D(cur_level, type);
    //    if (error != MB_SUCCESS) return error;
      }

    // Step 6: Update the global maps
    error = update_global_ahf(type, cur_level, deg);
    if (error != MB_SUCCESS) return error;

   // error = print_tags_2D(cur_level, type);
  //  if (error != MB_SUCCESS) return error;

    delete [] vbuffer;
    delete [] ent_buffer;

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::construct_hm_3D(int cur_level, int deg)
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
    EntityType type = mbImpl->type_from_handle(*(_incells.begin()));
    int cindex = type -1;
    int d = get_index_from_degree(deg);
    int ne = refTemplates[cindex][d].nv_edge;
    int nvf = refTemplates[cindex][d].nv_face;
    int nvtotal = refTemplates[cindex][d].total_new_verts;

    int index = ahf->get_index_from_type(*(_incells.begin()));
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
    std::vector<EntityHandle> conn;

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int cid = 0; cid < nents_prev; cid++)
      {
        conn.clear();
        for (int i=0; i<vtotal; i++)
          vbuffer[i] = 0;

        //EntityHandle of the working cell
        EntityHandle cell;
        if (cur_level)
          cell = level_mesh[cur_level-1].start_cell + cid;
        else
          cell = *_incells.begin() + cid;

        error = get_connectivity(cell, cur_level-1, conn);
        if (error != MB_SUCCESS) return error;

        //Step 1: Add vertices from the current level for the working face that will be used for subdivision.
        // Add the corners to vbuffer
        for (int i=0; i<(int)conn.size(); i++)
          {
            if (cur_level)
              vbuffer[i] =  level_mesh[cur_level].start_vertex + (conn[i]-level_mesh[cur_level-1].start_vertex);
            else
              vbuffer[i] = level_mesh[cur_level].start_vertex +  (conn[i]-*_inverts.begin()) ;
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

        std::cout<<"vbuffer = [";
        for (int i=0; i<vtotal; i++)
          std::cout<<" "<<vbuffer[i]<<", ";
        std::cout<<" ] "<<std::endl;

        //Step 2: Use the template to obtain the subentities. The coordinates and local ahf maps are also constructed.
        error = subdivide_cells(type, conn, cur_level, deg, vbuffer, vtotal, &count_ents, flag_verts, nverts_prev);
        if (error != MB_SUCCESS) return error;

        //Step 3: Update tracking information
        error = update_tracking_verts(cell, cur_level, deg, trackvertsC_edg, trackvertsC_face, vbuffer);
        if (error != MB_SUCCESS) return error;
      }

    //Step 4: Update the global maps
    error = update_global_ahf(type, cur_level, deg);
    if (error != MB_SUCCESS) return error;

    error = print_tags_3D(cur_level, type);
    if (error != MB_SUCCESS) return error;

    delete [] vbuffer;

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::subdivide_cells(EntityType type, std::vector<EntityHandle> &conn, int cur_level, int deg, EntityHandle *vbuffer, int vtotal,  int *count_ents, std::vector<int> &flag_verts, int nverts_prev)
  {
    //Subdivide the cell into children cells. Also update the sibling half-facet maps for the new cells.
    //Obtain incident edges/faces on this cell
    ErrorCode error;
    int index = ahf->get_index_from_type(*_incells.begin());
    int nvpc = ahf->lConnMap3D[index].num_verts_in_cell;
    int d = get_index_from_degree(deg);

    //Coordinates of the new vertices
    double *corner_coords = new double[nvpc*3];
    error = get_coordinates(conn, nvpc, cur_level, corner_coords);
    if (error != MB_SUCCESS) return error;

    error = compute_coordinates(cur_level, deg, type, vbuffer, vtotal, corner_coords, flag_verts, nverts_prev);
    if (error != MB_SUCCESS) return error;

    //Choose pattern. Special treatment is required for tetrahedral.
    int pat_id;
    if (type==MBTET)
      {
        int diag = find_shortest_diagonal_octahedron(cur_level, deg, vbuffer);
        pat_id = diag + 2;
      }
    else
      pat_id = type-1;

    //Connectivity of the children
    int etotal = refTemplates[pat_id][d].total_new_ents;
    EntityHandle  *ent_buffer = new EntityHandle[etotal];

    for (int i = 0; i < etotal; i++)
      {
        for (int k = 0; k < nvpc; k++)
          {
            int idx = refTemplates[type-1][d].ents_conn[i][k];
            level_mesh[cur_level].cell_conn[nvpc*count_ents[0]+k] = vbuffer[idx];
          }
        ent_buffer[i] =  level_mesh[cur_level].start_cell + count_ents[0];
        count_ents[0] += 1;
      }

    //Update local ahf maps
    error = update_local_ahf(deg, type, vbuffer, ent_buffer, etotal);
    if (error != MB_SUCCESS) return error;

    delete [] ent_buffer;
    delete [] corner_coords;

    return MB_SUCCESS;
  }

/*
  ErrorCode NestedRefine::subdivide_tets(std::vector<EntityHandle> conn, int cur_level, int deg, std::vector<EntityHandle> vbuffer, int *count_ents)
  {
    ErrorCode error;
    int index = get_index_from_type(*(_cells.begin()));
    int nvpc = lConnMap3D[index].num_verts_in_cell;

    double *corner_coords;
    error = get_coordinates(conn, nvpc, cur_level, corner_coords);
    if (error != MB_SUCCESS) return error;

    //Find or get the coordinates of the vertices in vbuffer
    int d = get_index_from_degree(deg);
    int vtotal = refTemplates[4][d].total_new_verts;

    error = compute_coordinates(cur_level, deg, MBTET, vbuffer, vtotal, corner_coords);
    if (error != MB_SUCCESS) return error;

    //Subdivide the tet into mixed entities = tets+octs
    int nents = refTemplates[4][d].num_ents;
    EntityHandle mixed_rep[nents][8];
    std::vector<int> nents_flag(nents, 0);
    std::vector<int> idx_buffer(nents,0);

    for (int i=0; i<nents; i++)
      {
        for (int j=0; j<8; j++)
          {
            mixed_rep[i][j] = vbuffer[refTemplates[4][d].ents_conn[i][j]];
            if ((j=4) && mixed_rep[i][j])
              nents_flag[i] = 1; //This is an oct
          };
      };

//    int count=0;
    //Add tets and tessellate octs based on shortest diagonal.
    for (int i=0; i<nents; i++)
      {
        idx_buffer[i] = count_ents[0];

        if (!nents_flag[i])//This is a tet
          {
            level_mesh[cur_level].cell_conn[4*count_ents[0]] = mixed_rep[i][0];
            level_mesh[cur_level].cell_conn[4*count_ents[0]+1] = mixed_rep[i][1];
            level_mesh[cur_level].cell_conn[4*count_ents[0]+2] = mixed_rep[i][2];
            level_mesh[cur_level].cell_conn[4*count_ents[0]+3] = mixed_rep[i][3];
            count_ents[0] += 1;
          }
        else
          {
            double ocoords[18];
            for (int j=0; j<6; j++)
              {
                int id = mixed_rep[i][j] -level_mesh[cur_level].start_vertex;
                ocoords[3*i] =  level_mesh[cur_level].coordinates[0][id] ;
                ocoords[3*i+1] =  level_mesh[cur_level].coordinates[1][id] ;
                ocoords[3*i+2] =  level_mesh[cur_level].coordinates[2][id] ;
              }

            int diag = find_shortest_diagonal_octahedron(ocoords);
            nents_flag[i] = diag;
            for (int j=0; j<4; j++)
              {
               for (int k=0; k<4; k++)
                 {
                   int id = oct_tessellation[diag].tet_conn[j][k];
                   level_mesh[cur_level].cell_conn[4*count_ents[0+k]] = mixed_rep[i][id];
                   count_ents[0] += 1;
                 }
              }
          }
      }


    //Construct the local ahf maps
    error = update_local_ahf(cur_level, deg, nents_flag, idx_buffer);
    if (error != MB_SUCCESS) return error;

    return MB_SUCCESS;
  } */

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
          //  std::cout<<"vflag["<<vbuffer[i] - vstart - nverts_prev<<"] = "<<vflag[vbuffer[i] - vstart - nverts_prev]<<std::endl;

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
         //   std::cout<<"vflag["<<vbuffer[i] - vstart - nverts_prev<<"] = "<<vflag[vbuffer[i] - vstart - nverts_prev]<<std::endl;
          }
      }
    return MB_SUCCESS;
  }

  /**********************************
   *          Update AHF maps           *
   * ********************************/

ErrorCode NestedRefine::update_local_ahf(int deg, EntityType type, EntityHandle *vbuffer, EntityHandle *ent_buffer, int etotal)
{
  ErrorCode error;
  int nhf, nv, total_new_verts;
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
      nhf = ahf->local_maps_2d(*_infaces.begin());
      nv = nhf;
      total_new_verts = refTemplates[type-1][d].total_new_verts;
    }
  else if (type == MBTET || type == MBHEX)
    {
      int index = ahf->get_index_from_type(*_incells.begin());
      nhf = ahf->lConnMap3D[index].num_faces_in_cell;
      nv =  ahf->lConnMap3D[index].num_verts_in_cell;
      total_new_verts = refTemplates[type-1][d].total_new_verts;
    }

  //Update the vertex to half-facet map
  for (int i=0; i<total_new_verts; i++)
    {
      EntityHandle ent;
      int lid;

      EntityHandle vid = vbuffer[i+nv];
      error = ahf->get_incident_tag(type, vid, &ent, &lid);
      if (error != MB_SUCCESS) return error;

      if (ent)
        continue;

      int id = refTemplates[type-1][d].v2hf[i+nv][0]-1;
      ent = ent_buffer[id];
      lid = refTemplates[type-1][d].v2hf[i+nv][1];

      error = ahf->set_incident_tag(type, vid, &ent, &lid);
      if (error != MB_SUCCESS) return error;
    }

  //Update the sibling half-facet map
  for (int i=0; i< etotal; i++)
    {
      EntityHandle  *sib_entids = new EntityHandle[nhf];
      int *sib_lids = new int[nhf];

      error = ahf->get_sibling_tag(type, ent_buffer[i], &sib_entids[0], &sib_lids[0]);
      if (error != MB_SUCCESS) return error;

      for (int l=0; l< nhf; l++)
        {
          if (sib_entids[l])
            continue;

          // Fill out the sibling values
          int id = refTemplates[type-1][d].ents_opphfs[i][2*l];
          if (id)
            {
              sib_entids[l] = ent_buffer[id-1];
              sib_lids[l] = refTemplates[type-1][d].ents_opphfs[i][2*l+1];
            }
          else
            {
              sib_entids[l] = 0;
              sib_lids[l] = 0;
            }
        }

      error = ahf->set_sibling_tag(type, ent_buffer[i], &sib_entids[0], &sib_lids[0]);
      if  (error != MB_SUCCESS) return error;

      for (int l=0; l< nhf; l++)
        {
          if (sib_entids[l]){
              EntityHandle *set_entids = new EntityHandle[nhf];
              int *set_lids = new int[nhf];

              error = ahf->get_sibling_tag(type, sib_entids[l], &set_entids[0], &set_lids[0]);
              if (error != MB_SUCCESS) return error;

              set_entids[sib_lids[l]] = ent_buffer[i];
              set_lids[sib_lids[l]] = l;

              error = ahf->set_sibling_tag(type, sib_entids[l], set_entids, set_lids);
              if (error != MB_SUCCESS) return error;

              delete [] set_entids;
              delete [] set_lids;

            }
        }

      delete [] sib_entids;
      delete [] sib_lids;

    }
  return MB_SUCCESS;
}

/*ErrorCode NestedRefine::update_local_ahf(int cur_level, int deg, std::vector<int> nents_flag, std::vector<int> idx_buffer)
{
  ErrorCode error;
  int d = get_index_from_degree(deg);

  // Update the local ahf maps for tet mesh.
  int nents = refTemplates[4][d].num_ents;
  EntityHandle start_cell = level_mesh[cur_level].start_cell;
  int nfpt = 4;

  for (int i=0; i<nents; i++) //Loop over entities of the mixed representation
    {
      if (!nents_flag[i]) // This is a tet
        {
          std::vector<EntityHandle> sib_cids;
          std::vector<int>  sib_lfids;

          EntityHandle cid = start_cell + idx_buffer[i];

          error = get_sibling_tag(MBTET, &cid, 1, sib_cids, sib_lfids);
          if (MB_SUCCESS != error) return error;


          for (int lf = 0; lf < nfpt; lf++) //Loop over local faces
            {
              if (sib_cids[lf])
                continue;

              int tid = refTemplates[4][d].ents_opphfs[i][2*lf];
              int lid = refTemplates[4][d].ents_opphfs[i][2*lf+1];

              if (!nents_flag[tid]) // Opposite cell is also a tet
                {
                  sib_cids[lf] = start_cell + idx_buffer[tid];
                  sib_lfids[lf] = lid;
                }
              else // Opposite cell is an oct
                {
                  int diag = nents_flag[tid]; //This is the chosen diagonal
                  int tet_id = oct_tessellation[diag].olfid_to_tlfid[lid][0];
                  int tet_lfid = oct_tessellation[diag].olfid_to_tlfid[lid][1];
                  sib_cids[lf] = start_cell + idx_buffer[tid] + tet_id;
                  sib_lfids[lf] = tet_lfid;
                }
            }

          // Set the sibling half-face map for the current tet
          error = set_sibling_tag(MBTET, &cid, 1, sib_cids, sib_lfids);
          if (MB_SUCCESS != error) return error;

          std::vector<EntityHandle> set_cids;
          std::vector<int>  set_lfids;
          //Set the sibling info for the sibling tet
          for (int lf = 0; lf< nfpt; lf++){
              if (sib_cids[lf]) {
                  error = get_sibling_tag(MBTET, &sib_cids[lf], 1, set_cids, set_lfids);
                  if (MB_SUCCESS != error) return error;

                  set_cids[sib_lfids[lf]] = cid;
                  set_lfids[sib_lfids[lf]] = lf;

                  error = set_sibling_tag(MBTET, &sib_cids[lf], 1, set_cids, set_lfids);
                  if (MB_SUCCESS != error) return error;
                }
            }
        }
      else //This is originally an oct
        {
          int diag = nents_flag[i];
          for (int t = 0; t< 4; t++) // Loop over the 4 sub-tets
            {
              std::vector<EntityHandle> sib_cids;
             std::vector<int>  sib_lfids;

              EntityHandle cid = start_cell + idx_buffer[i] + t;
              error = get_sibling_tag(MBTET, &cid, 1, sib_cids, sib_lfids);
              if (MB_SUCCESS != error) return error;

              for (int l= 0; l < 4; l++) // Loop over the four faces of the working sub tet
                {
                  if (sib_cids[l])
                    continue;

                  int olid = oct_tessellation[diag].tlfid_to_olfid[t][l];
                  if (olid) // This face maps to a face of the parent octahedron
                  {
                    int tet_id = refTemplates[4][d].ents_opphfs[i][2*olid];
                    sib_cids[l] = start_cell + idx_buffer[tet_id];
                    sib_lfids[l] = refTemplates[4][d]. ents_opphfs[i][2*olid+1];
                  }
                  else
                    {
                      int tet_id = oct_tessellation[diag].tet_opphfs[t][2*l];
                      sib_cids[l] = start_cell + idx_buffer[i] + tet_id;
                      sib_lfids[l] = oct_tessellation[diag].tet_opphfs[t][2*l+1];
                    }
                }

              // Set the sibling half-face map for the current tet
              error = set_sibling_tag(MBTET, &cid, 1, sib_cids, sib_lfids);
              if (MB_SUCCESS != error) return error;

              std::vector<EntityHandle> set_cids;
              std::vector<int>  set_lfids;

              //Set the sibling info for the sibling tet
              for (int lf = 0; lf< nfpt; lf++){
                  if (sib_cids[lf]) {
                      error = get_sibling_tag(MBTET, &sib_cids[lf], 1, set_cids, set_lfids);
                      if (MB_SUCCESS != error) return error;

                      set_cids[sib_lfids[lf]] = cid;
                      set_lfids[sib_lfids[lf]] = lf;

                      error = set_sibling_tag(MBTET, &sib_cids[lf], 1, set_cids, set_lfids);
                      if (MB_SUCCESS != error) return error;
                    }
                }
            }
        }
    }
return MB_SUCCESS;
}
*/

ErrorCode NestedRefine::update_global_ahf(EntityType type, int cur_level, int deg)
{

  ErrorCode error;

  //Get the number of half-facets and number of children of each type
  if (type == MBEDGE)
    {
      error = update_global_ahf_1D(cur_level, deg);
      if (error != MB_SUCCESS) return error;
    }
  else if (type == MBTRI || type == MBQUAD)
    {
      error = update_global_ahf_2D(cur_level, deg);
      if (error != MB_SUCCESS) return error;
    }
  else
    {
      error = update_global_ahf_3D(cur_level, deg);
      if (error != MB_SUCCESS) return error;
    }

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

   //Update the vertex to half-facet maps for duplicate vertices
   for (int i=0; i<nverts_prev; i++)
     {
       //Vertex id in the previous mesh
       EntityHandle vid;
       if (cur_level)
         vid = level_mesh[cur_level-1].start_vertex + i;
       else
         vid = _inverts[i];

       //Get the incident half-vert in the previous mesh
       EntityHandle inci_ent=0;
       int inci_lid=0;

       error = ahf->get_incident_tag(MBEDGE, vid, &inci_ent, &inci_lid);
       if (error != MB_SUCCESS) return error;

       // Obtain the corresponding incident child in the current mesh
       int lvid = get_local_vid(vid, inci_ent, cur_level-1);
       int chid = refTemplates[0][d].v2hf[lvid][0]-1;

       int pid;
       if (cur_level)
         pid = inci_ent - level_mesh[cur_level-1].start_edge;
       else
         pid = inci_ent - *_inedges.begin();

       int ind = nchilds*pid;

       EntityHandle child_ent = level_mesh[cur_level].start_edge + ind+chid ;
       int child_lid = refTemplates[0][d].v2hf[lvid][1];

       EntityHandle cur_vid = level_mesh[cur_level].start_vertex + i;
       error = ahf->set_incident_tag(MBEDGE, cur_vid, &child_ent, &child_lid);
       if (error != MB_SUCCESS) return error;
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

       error = ahf->get_sibling_tag(MBEDGE, ent, &sib_entids[0], &sib_lids[0]);
       if (error != MB_SUCCESS) return error;

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

           error = ahf->get_sibling_tag(MBEDGE, child_ent, &sib_childs[0], &sib_chlids[0]);
           if (error != MB_SUCCESS) return error;

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

           error = ahf->set_sibling_tag(MBEDGE, child_ent, &sib_childs[0], &sib_chlids[0]);
           if (error != MB_SUCCESS) return error;

           delete [] sib_childs;
           delete [] sib_chlids;
         }
       delete [] sib_entids;
       delete [] sib_lids;
     }

   return MB_SUCCESS;
 }

 ErrorCode NestedRefine::update_global_ahf_2D(int cur_level, int deg)
 {
   ErrorCode error;

   EntityType type = mbImpl->type_from_handle(*_infaces.begin());
   int nhf, nchilds, nverts_prev, nents_prev;

   nhf = ahf->local_maps_2d(*_infaces.begin());
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

   //Update the vertex to half-edge maps for old/duplicate vertices
   for (int i=0; i<nverts_prev; i++)
     {
       //Vertex id in the previous mesh
       EntityHandle vid;
       if (cur_level)
         vid = level_mesh[cur_level-1].start_vertex + i;
       else
         vid = _inverts[i];

       //Get the incident half-vert in the previous mesh
       EntityHandle inci_ent=0;
       int inci_lid=0;

       error = ahf->get_incident_tag(type, vid, &inci_ent, &inci_lid);
       if (error != MB_SUCCESS) return error;

       // Obtain the corresponding incident child in the current mesh
       int lvid = get_local_vid(vid, inci_ent, cur_level-1);
       int chid = refTemplates[type-1][d].v2hf[lvid][0]-1;

       int pid;
       if (cur_level)
         pid = inci_ent - level_mesh[cur_level-1].start_face;
       else
         pid = inci_ent - *_infaces.begin();

       int ind = nchilds*pid;

       EntityHandle child_ent = level_mesh[cur_level].start_face + ind+chid ;
       int child_lid = refTemplates[type-1][d].v2hf[lvid][1];

       EntityHandle cur_vid = level_mesh[cur_level].start_vertex + i;

       error = ahf->set_incident_tag(type, cur_vid, &child_ent, &child_lid);
       if (error != MB_SUCCESS) return error;
     }


   EntityHandle fedge[2];
   int * next = new int[nhf];
   int  * prev = new int[nhf];
   error = ahf->local_maps_2d(nhf, next, prev);
   if (error != MB_SUCCESS) return error;

   //Update the sibling half-facet maps across entities
   for (int i=0; i< nents_prev; i++)
     {
       EntityHandle ent;
       if (cur_level)
         ent = level_mesh[cur_level-1].start_face + i;
       else
         ent = _infaces[i];

      // std::cout<<"Working ent = "<<ent<<std::endl;

       std::vector<EntityHandle> fid_conn;
       error = get_connectivity(ent, cur_level-1, fid_conn);
       if (MB_SUCCESS != error) return error;

       EntityHandle *sib_entids = new EntityHandle[nhf];
       int *sib_lids = new int[nhf];

       error = ahf->get_sibling_tag(type, ent, &sib_entids[0], &sib_lids[0]);
       if (error != MB_SUCCESS) return error;

       int id, idx;

       for (int l=0; l < nhf; l++)
         {
           if (!sib_entids[l])
             continue;

           fedge[0] = fid_conn[l];
           fedge[1] = fid_conn[next[l]];

           EntityHandle sfid = sib_entids[l];
           int slid = sib_lids[l];

           std::vector<EntityHandle> conn;
           error = get_connectivity(sfid, cur_level-1, conn);
           if (MB_SUCCESS != error) return error;

           bool orient;
           if ((fedge[0] == conn[slid])&&(fedge[1] == conn[next[slid]]))
             orient = true;
           else if ((fedge[1] == conn[slid])&&(fedge[0] == conn[next[slid]]))
             orient = false;

           //Find the childrens incident on the half-facet
           int nch = refTemplates[type-1][d].ents_on_pent[l][0];
           idx = nchilds*i;

           //Loop over all the incident childrens
           for (int k=0; k<nch; k++)
             {
               id = refTemplates[type-1][d].ents_on_pent[l][k+1]-1;
               EntityHandle child_ent = level_mesh[cur_level].start_face + idx+id ;
               int ch_lid = l;

               //Find the sibling of the child
               EntityHandle *sib_childs = new EntityHandle[nhf];
               int *sib_chlids = new int[nhf];

               error = ahf->get_sibling_tag(type, child_ent, &sib_childs[0], &sib_chlids[0]);
               if (error != MB_SUCCESS) return error;

               //If the sibling already exists, dont do anything
               if (sib_childs[ch_lid])
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
               sib_childs[ch_lid] = psib_child;
               sib_chlids[ch_lid] = psib_chlid;

               error = ahf->set_sibling_tag(type, child_ent, &sib_childs[0], &sib_chlids[0]);
               if (error != MB_SUCCESS) return error;

               delete [] sib_childs;
               delete [] sib_chlids;
             }
         }
       delete [] sib_entids;
       delete [] sib_lids;
     }

   delete [] next;
   delete [] prev;

   return MB_SUCCESS;

 }

ErrorCode NestedRefine::update_global_ahf_3D(int cur_level, int deg)
{
  ErrorCode error;
  int nhf, nchilds, nverts_prev, nents_prev;

  EntityType type = mbImpl->type_from_handle(*_incells.begin());
  int index = ahf->get_index_from_type(*_incells.begin());
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

  //Step 1: Update the V2HF maps for old/duplicate vertices
  for (int i=0; i<nverts_prev; i++)
    {
      //Vertex id in the previous mesh
      EntityHandle vid;
      if (cur_level)
        vid = level_mesh[cur_level-1].start_vertex + i;
      else
        vid = _inverts[i];

      //Get the incident half-vert in the previous mesh
      EntityHandle inci_ent=0;
      int inci_lid=0;

      error = ahf->get_incident_tag(type, vid, &inci_ent, &inci_lid);
      if (error != MB_SUCCESS) return error;

      // Obtain the corresponding incident child in the current mesh
      int lvid = get_local_vid(vid, inci_ent, cur_level-1);
      int chid = refTemplates[type-1][d].v2hf[lvid][0]-1;

      int pid;
      if (cur_level)
        pid = inci_ent - level_mesh[cur_level-1].start_cell;
      else
        pid = inci_ent - *_incells.begin();

      int ind = nchilds*pid;

      EntityHandle child_ent = level_mesh[cur_level].start_cell + ind+chid ;
      int child_lid = refTemplates[type-1][d].v2hf[lvid][1];

      EntityHandle cur_vid = level_mesh[cur_level].start_vertex + i;

      error = ahf->set_incident_tag(type, cur_vid, &child_ent, &child_lid);
      if (error != MB_SUCCESS) return error;
    }

  //Step 2: Update SIBHFS maps
  for (int i=0; i< nents_prev; i++)
    {
      EntityHandle ent;
      if (cur_level)
        ent = level_mesh[cur_level-1].start_cell + i;
      else
        ent = _incells[i];

      // std::cout<<"Working ent = "<<ent<<std::endl;

      EntityHandle *sib_entids = new EntityHandle[nhf];
      int *sib_lids = new int[nhf];

      error = ahf->get_sibling_tag(type, ent, &sib_entids[0], &sib_lids[0]);
      if (error != MB_SUCCESS) return error;

      int id, idx;

      for (int l=0; l < nhf; l++)
        {
          if (!sib_entids[l])
            continue;

          //Get the number of children incident on this half-face
          int nch = refTemplates[type-1][d].ents_on_pent[l][0];

          //Get the order of children indices incident on this half-face
          int *id_sib = new int[nch];
          error = reorder_indices(cur_level, deg, ent, l, sib_entids[l], sib_lids[l], 1, id_sib);
          if (error != MB_SUCCESS) return error;

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
              EntityHandle *child_sibents= new EntityHandle[nhf];
              int *child_siblids = new int[nhf];

              error = ahf->get_sibling_tag(type, child_ent, &child_sibents[0], &child_siblids[0]);
              if (error != MB_SUCCESS) return error;

              //If the sibling already exists, dont do anything
              if (child_sibents[child_lid])
                continue;

              //Get the correponding child of the sibling of the current parent
              // We have already computed the order the children on incident corresponding to the working half-face
              id = refTemplates[type-1][d].ents_on_pent[plid][id_sib[k]]-1;

              EntityHandle psib_child = level_mesh[cur_level].start_cell + sidx+id ;
              int psib_chlid = plid;

              //Set the siblings of children incident on current half-face
              child_sibents[child_lid] = psib_child;
              child_siblids[child_lid] = psib_chlid;

              error = ahf->set_sibling_tag(type, child_ent, child_sibents, child_siblids);
              if (error != MB_SUCCESS) return error;

              //Set the sibling of the sibling of the children to the children
              EntityHandle *sibch_ents= new EntityHandle[nhf];
              int *sibch_lids = new int[nhf];

              error = ahf->get_sibling_tag(type, psib_child, &sibch_ents[0], &sibch_lids[0]);
              if (error != MB_SUCCESS) return error;

              sibch_ents[psib_chlid] = child_ent;
              sibch_lids[psib_chlid]  = child_lid;

              error = ahf->set_sibling_tag(type, psib_child, sibch_ents, sibch_lids);
              if (error != MB_SUCCESS) return error;

              delete [] child_sibents;
              delete [] child_siblids;
              delete [] sibch_ents;
              delete [] sibch_lids;
            }

          delete [] id_sib;
        }

      delete [] sib_entids;
      delete [] sib_lids;

    }

  return MB_SUCCESS;
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
      error = mbImpl->get_coords(_inverts, vcoords);
      if (error != MB_SUCCESS) return error;

      for (int i = 0; i < nverts_in; i++)
        {
          // std::cout<<"x= "<<vcoords[3*i]<<", y = "<<vcoords[3*i+1]<<", z = "<<vcoords[3*i+2]<<std::endl;
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

  int index = ahf->get_index_from_type(*(_incells.begin()));
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

      error = ahf->get_up_adjacencies_edg_3d(cid, i, inc_cids,  &inc_leids, &inc_orient);
      if (error != MB_SUCCESS) return error;

      if (inc_cids.size() ==1)
        continue;

      //Add the vertices to the edges of the incident cells
      for (int k=1; k<(int)inc_cids.size(); k++)
        {
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
          error = ahf->get_up_adjacencies_face_3d(cid, i, sib_cids, &sib_lfids);
          if (error != MB_SUCCESS) return error;

          if (sib_cids.size()==1)
            continue;

          int *id_sib = new int[nvf];
          error = reorder_indices(cur_level, deg, sib_cids[1], sib_lfids[1], cid, i, 0, id_sib);
          if (error != MB_SUCCESS) return error;

        /*  for (int l=0; l<nvf; l++)
            std::cout<<"face_buf[ "<<l<<"] = "<<face_vbuf[l]<<", ";
          std::cout<<std::endl;

          std::cout<<"id_sib = [";
          for (int l=0; l<nvf; l++)
            std::cout<<id_sib[l]<<" ";
          std::cout<<" ] "<<std::endl;*/

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
  int idx = ahf->get_index_from_type(*_incells.begin());
  int nvF = ahf->lConnMap3D[idx].hf2v_num[lfid];
  int nco = permutation[nvF-3].num_comb;

  if (!index && nvF==4 && deg==2)
    {
      id_sib[0] = 1;
    }
  else
    {
      //Get connectivity of the cell and its sibling cell
      std::vector<EntityHandle> conn, sib_conn;
      error = get_connectivity(cell, cur_level, conn);
      if (error != MB_SUCCESS) return error;

      error = get_connectivity(sib_cell, cur_level, sib_conn);
      if (error != MB_SUCCESS) return error;

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

   /*   std::cout<<"lface = [ ";
      for (int i=0; i<nvF; i++)
        std::cout<<lface[i]<<" ";
      std::cout<<" ] "<<std::endl;
      std::cout<<"lface_sib = [ ";
      for (int i=0; i<nvF; i++)
        std::cout<<lface_sib[i]<<" ";
      std::cout<<" ] "<<std::endl;*/

      //Find the combination
      int c;
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

      //Get the ordered indices
      if (((!index)&&(nvF==4)&&(deg==3))||(deg==2))
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

ErrorCode NestedRefine::count_subentities(EntityHandle set, int cur_level, int *nedges, int *nfaces)
{
  ErrorCode error;

  if (cur_level>=0){
      Range edges, faces, cells;

      error = mbImpl->get_entities_by_dimension( set, 1, edges);
      if (MB_SUCCESS != error) return error;

      error = mbImpl->get_entities_by_dimension( set, 2, faces);
      if (MB_SUCCESS != error) return error;

      error = mbImpl->get_entities_by_dimension( set, 3, cells);
      if (MB_SUCCESS != error) return error;

      error = ahf->count_subentities(edges, faces, cells, nedges, nfaces);
      if (MB_SUCCESS != error) return error;
    }
  else
    {
      error = ahf->count_subentities(_inedges, _infaces, _incells, nedges, nfaces);
      if (MB_SUCCESS != error) return error;
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::get_octahedron_corner_coords(int cur_level, int deg, EntityHandle *vbuffer, double *ocoords)
{
  int lid[6];

  if (deg==2)
    {
      lid[0] = 5; lid[1] = 8; lid[2] = 9;
      lid[3] = 6; lid[4] = 4; lid[5] = 7;
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
  assert(error == MB_SUCCESS);

  int diag_map[6] = {1,3,2,4,5,0};
  double length = std::numeric_limits<double>::max();

  int diag = 0;
  double dlen = 0;
  double x = 0, y = 0, z = 0;

  for (int d=0; d<3; d++)
    {
      int id1 = diag_map[2*d];
      int id2 = diag_map[2*d+1];
      x = coords[3*id1] - coords[3*id2];
      y = coords[3*id1+1] - coords[3*id2+1];
      z = coords[3*id1+2] - coords[3*id2+2];
      dlen =  sqrt(x*x+y*y+z*z);
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

  error = get_connectivity(ent, level, conn);
  assert(error == MB_SUCCESS);
//  if (error != MB_SUCCESS) return error;

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
  int d;
  if (degree==2)
    d=0;
  else if (degree==3)
    d = 1;
  else if (degree == 5)
    d=2;
  else
    std::cout<<"Unsupported degree"<<std::endl;
  return d;
}

ErrorCode NestedRefine::print_tags_1D(int level)
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
      error = ahf->get_incident_tag(MBEDGE, vid, &eid, &lvid);
      if (MB_SUCCESS != error) return error;

      std::cout<<"For vertex = "<<vid<<"::Incident halfvertex "<<eid<<"  "<<lvid<<std::endl;
    }

  //SIBHVS
  std::cout<<"start_edge = "<<start_edge<<std::endl;
  std::cout<<"<SIBHVS_EID,SIBHVS_LVID>"<<std::endl;
  for (int i=0; i<ne; i++)
    {
      EntityHandle ent = start_edge+i;

      EntityHandle eid[2];  int lvid[2];
      error = ahf->get_sibling_tag(MBEDGE, ent, &eid[0], &lvid[0]);
      if (MB_SUCCESS != error) return error;

      std::cout<<"<"<<eid[0]<<","<<lvid[0]<<">"<<"      "<<"<"<<eid[1]<<","<<lvid[1]<<">"<<std::endl;
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::print_tags_2D(int level, EntityType type)
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
      error = ahf->get_incident_tag(type, vid, &fid, &leid);
      if (MB_SUCCESS != error) return error;

      std::cout<<"For vertex = "<<vid<<"::Incident halfedge "<<fid<<"  "<<leid<<std::endl;
    }

  //SIBHES
  std::cout<<"start_face = "<<start_face<<std::endl;
  std::cout<<"<SIBHES_FID,SIBHES_LEID>"<<std::endl;
  int nepf = ahf->local_maps_2d(*_infaces.begin());

  for (int i=0; i<nf; i++)
    {
      EntityHandle ent = start_face+i;

      EntityHandle *fid = new EntityHandle[nepf];
      int *leid = new int [nepf];
      error = ahf->get_sibling_tag(type, ent, &fid[0], &leid[0]);
      if (MB_SUCCESS != error) return error;

      for (int j=0; j<nepf; j++){
          std::cout<<"<"<<fid[j]<<","<<leid[j]<<">"<<"      ";
      }
      std::cout<<std::endl;

      delete [] fid;
      delete [] leid;
    }

  return MB_SUCCESS;
}

ErrorCode NestedRefine::print_tags_3D(int level, EntityType type)
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
      error = ahf->get_incident_tag(type, vid, &cid, &lfid);
      if (MB_SUCCESS != error) return error;

      std::cout<<"For vertex = "<<vid<<"::Incident halfface "<<cid<<"  "<<lfid<<std::endl;
    }

  //SIBHFS
  std::cout<<"start_cell = "<<start_cell<<std::endl;
  std::cout<<"<SIBHFS_CID,SIBHFS_LFID>"<<std::endl;
  int index = ahf->get_index_from_type(start_cell);
  int nfpc = ahf->lConnMap3D[index].num_faces_in_cell;

  for (int i=0; i<nc; i++)
    {
      EntityHandle ent = start_cell+i;

      EntityHandle *cid = new EntityHandle[nfpc];
      int *lfid = new int [nfpc];
      error = ahf->get_sibling_tag(type, ent, &cid[0], &lfid[0]);
      if (MB_SUCCESS != error) return error;

      for (int j=0; j<nfpc; j++){
          std::cout<<"<"<<cid[j]<<","<<lfid[j]<<">"<<"      ";
      }
      std::cout<<std::endl;

      delete [] cid;
      delete [] lfid;
    }

  return MB_SUCCESS;
}

}//namesapce moab

