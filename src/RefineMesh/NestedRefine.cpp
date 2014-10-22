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

  NestedRefine::NestedRefine(Core *impl):HalfFacetRep(impl)
  {
    ErrorCode error;
    error = initialize();
    assert(error==MB_SUCCESS);
  }
  /************************************************************
   *     Interface Functions                                  *
   ************************************************************/

  ErrorCode NestedRefine::generate_mesh_hierarchy(int *level_degrees, int num_level, EntityHandle *hm_set)
  {
    ErrorCode error;
    int *hmest = new int[4*num_level];
    error = estimate_hm_storage(level_degrees, num_level, hmest);
    if (error != MB_SUCCESS) return error;

    error = generate_hm(level_degrees, num_level, hmest, hm_set);
    if (error != MB_SUCCESS) return error;

  //  error = add_entities( hm_set, num_level, hmest);
 //   if (error != MB_SUCCESS) return error;

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::get_connectivity(EntityHandle ent, int level, std::vector<EntityHandle> &conn)
  {
    ErrorCode error;
    EntityType type = mb->type_from_handle(ent);
    EntityHandle start_ent ;
    if (level>= 0)
      {
        if (type == MBEDGE)
          {
            start_ent = level_mesh[level].start_edge;
            conn[0] = level_mesh[level].edge_conn[2*(ent-start_ent)];
            conn[1] = level_mesh[level].edge_conn[2*(ent-start_ent)+1];
          }
        else if (type == MBTRI || type == MBQUAD)
          {
            int num_corners = local_maps_2d(ent);
            start_ent = level_mesh[level].start_face;
            for (int i=0; i<num_corners; i++)
              conn[i] = level_mesh[level].face_conn[num_corners*(ent-start_ent)+i];
          }
        else if (type == MBTET || type == MBPRISM || type == MBHEX)
          {
            int index = get_index_from_type(*_cells.begin());
            int num_corners = lConnMap3D[index].num_verts_in_cell;
            start_ent = level_mesh[level].start_cell;
            for (int i=0; i<num_corners; i++)
              conn[i] = level_mesh[level].cell_conn[num_corners*(ent-start_ent)+i];
          }
        else
          {
            std::cout<<"Unsupport element type"<<std::endl;
            return MB_FAILURE;
          }
      }
      else
      {
        error = mb->get_connectivity(&ent, 1, conn);
        if (error != MB_SUCCESS) return error;
      }

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::get_coordinates(std::vector<EntityHandle> verts, int num_verts, int cur_level, double *coords)
  {
    EntityHandle vstart = level_mesh[cur_level].start_vertex;
    for (int i=0; i< num_verts; i++)
      {
        EntityHandle vid = verts[i];
        coords[3*i] = level_mesh[cur_level].coordinates[vid-vstart][0];
        coords[3*i+1] =  level_mesh[cur_level].coordinates[vid-vstart][1];
        coords[3*i+2] =  level_mesh[cur_level].coordinates[vid-vstart][2];
      }
    return MB_SUCCESS;
  }

  /***********************************************
   *  Basic functionalities: generate HM         *
   ***********************************************/

  ErrorCode NestedRefine::estimate_hm_storage(int *level_degrees, int num_level, int *hmest)
  {
    ErrorCode error;

    //Obtain the size of input mesh.
    int nverts_in = _verts.size();
    int nedges_in = _edges.size();
    int nfaces_in = _faces.size();
    int ncells_in = _cells.size();

    //Estimate mesh size of level 1 mesh.
    int nedges=0, nfaces=0;
    error = count_subentities(_edges, _faces, _cells, &nedges, &nfaces);
    if (error != MB_SUCCESS) return error;

    int d = get_index_from_degree(level_degrees[0]);
    int nverts = refTemplates[MBEDGE-1][d].nv_edge*nedges;
    hmest[0] = nverts_in + nverts;
    hmest[1] = nedges_in*refTemplates[MBEDGE-1][d].total_new_ents;
    hmest[2] = 0;
    hmest[3] = 0;

    int findex, cindex;
    if (_faces.size())
      {
        findex = mb->type_from_handle(*(_faces.begin()))-1;
        nverts += refTemplates[findex][d].nv_face*nfaces;
        hmest[2] = nfaces_in*refTemplates[findex][d].total_new_ents;
      }
    if (_cells.size())
      {
        cindex = mb->type_from_handle(*(_cells.begin()))-1;
        nverts +=  refTemplates[cindex][d].nv_cell*ncells_in;
        hmest[3] = ncells_in*refTemplates[cindex][d].total_new_ents;
      }

    int nverts_prev, nedges_prev, nfaces_prev, ncells_prev;
    for (int l=1; l<num_level; l++)
      {
        nverts_prev = hmest[4*(l-1)];
        nedges_prev = hmest[4*(l-1)+1];
        nfaces_prev = hmest[4*(l-1)+2];
        ncells_prev = hmest[4*(l-1)+3];

        d = get_index_from_degree(level_degrees[l]);
        nverts = refTemplates[MBEDGE-1][d].nv_edge*nedges;
        hmest[4*l] = nverts_prev + nverts;
        hmest[4*l+1] = nedges_prev*refTemplates[MBEDGE-1][d].total_new_ents;
        hmest[4*l+2] = 0;
        hmest[4*l+3] = 0;

        if (_faces.size())
          {
            nverts += refTemplates[findex][d].nv_face*nfaces;
            hmest[4*l+2] = nfaces_prev*refTemplates[findex][d].total_new_ents;
          }
        if (_cells.size())
          {
            nverts +=  refTemplates[cindex][d].nv_cell*ncells_in;
            hmest[4*l+3] = ncells_prev*refTemplates[cindex][d].total_new_ents;
          }
      }

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::create_hm_storage_single_level(EntityHandle *set, int cur_level, int *estL)
  {
    //Obtain chunks of memory for the current level. Add them to a particular meshset.
    EntityHandle set_handle;
    ErrorCode error = mb->create_meshset(MESHSET_SET, set_handle);
    if (error != MB_SUCCESS) return error;
    *set = set_handle;

    ReadUtilIface *read_iface;
    error = mb->query_interface(read_iface);
    if (error != MB_SUCCESS) return error;

    //Vertices
    error = read_iface->get_node_coords(3, estL[4*cur_level], 0,  level_mesh[cur_level].start_vertex , level_mesh[cur_level].coordinates);
    if (error != MB_SUCCESS) return error;
    level_mesh[cur_level].num_verts = estL[4*cur_level];

    EntityHandle *varray = new EntityHandle[estL[4*cur_level]];
    for (int i=0; i< estL[4*cur_level]; i++)
      varray[i] = level_mesh[cur_level].start_vertex + i;

    error = mb->add_entities(*set, varray, estL[4*cur_level]);
    if (error != MB_SUCCESS) return error;

    // Edges
    if (estL[4*cur_level+1])
      {
        error = read_iface->get_element_connect(estL[4*cur_level+1], 2, MBEDGE, 0, level_mesh[cur_level].start_edge, level_mesh[cur_level].edge_conn);
        if (error != MB_SUCCESS) return error;
        level_mesh[cur_level].num_edges = estL[4*cur_level+1];

        EntityHandle *earray = new EntityHandle[4*cur_level+1];
        for (int i=0; i< estL[4*cur_level+1]; i++)
          earray[i] = level_mesh[cur_level].start_edge + i;

        error = mb->add_entities(*set, earray, estL[4*cur_level+1]);
        if (error != MB_SUCCESS) return error;
      }

    //debug
    std::cout<<"estL[0] ="<<estL[4*cur_level]<<", estL[1] = "<<estL[4*cur_level+1]<<", estL[2] = "<<estL[4*cur_level+2]<<", estL[3] = "<<estL[4*cur_level+3]<<std::endl;
    Range verts, ents;

    error = mb->get_entities_by_type(*set, MBVERTEX, verts);
    error = mb->get_entities_by_type(*set, MBEDGE, ents);

    if (verts.empty() || ents.empty())
      std::cout<<"something not right"<<std::endl;

    for (Range::iterator v = verts.begin(); v!= verts.end(); v++)
      std::cout<<"v["<<(*v-*verts.begin())<<"] = "<<*v<<std::endl;

    for (Range::iterator e = ents.begin(); e!= ents.end(); e++)
      std::cout<<"ent["<<(*e-*ents.begin())<<"] = "<<*e<<std::endl;



    //Faces
    if (estL[4*cur_level+2])
      {
        EntityType type = mb->type_from_handle(*(_faces.begin()));
        int nvpf = local_maps_2d(*(_faces.begin()));
        error = read_iface->get_element_connect(estL[4*cur_level+2], nvpf, type, 0, level_mesh[cur_level].start_face, level_mesh[cur_level].face_conn);
        if (error != MB_SUCCESS) return error;
        level_mesh[cur_level].num_faces = estL[4*cur_level+2];

        EntityHandle *farray = new EntityHandle[4*cur_level+2];
        for (int i=0; i< estL[4*cur_level+2]; i++)
          farray[i] = level_mesh[cur_level].start_face + i;

        error = mb->add_entities(*set, farray, estL[4*cur_level+2]);
        if (error != MB_SUCCESS) return error;
      }

    //Cells
    if (estL[4*cur_level+3])
      {
        EntityType type = mb->type_from_handle(*(_cells.begin()));
        int index = get_index_from_type(*_cells.begin());
        int nvpc = lConnMap3D[index].num_verts_in_cell;
        error = read_iface->get_element_connect(estL[4*cur_level+3], nvpc, type, 0, level_mesh[cur_level].start_cell, level_mesh[cur_level].cell_conn);
        if (error != MB_SUCCESS) return error;
        level_mesh[cur_level].num_cells = estL[4*cur_level+3];

        EntityHandle *carray = new EntityHandle[4*cur_level+3];
        for (int i=0; i< estL[4*cur_level+3]; i++)
          carray[i] = level_mesh[cur_level].start_cell + i;

        error = mb->add_entities(*set, carray, estL[4*cur_level+3]);
        if (error != MB_SUCCESS) return error;
      }

    return MB_SUCCESS;
  }


  ErrorCode NestedRefine::add_entities(EntityHandle *hm_set, int num_level, int *estL)
  {
    //Add the entities created by the hierarchy to its set

    ErrorCode error;
    for (int l=0; l< num_level; l++)
      {

        // Vertices
        error = mb->add_entities(hm_set[l], &level_mesh[l].start_vertex, estL[4*l]);
        if (error != MB_SUCCESS) return error;

        // Edges
        if (estL[4*l+1])
          {
            error = mb->add_entities(hm_set[l], &level_mesh[l].start_edge, estL[4*l+1]);
            if (error != MB_SUCCESS) return error;
          }

        //Faces
        if (estL[4*l+2])
          {
            error = mb->add_entities(hm_set[l], &level_mesh[l].start_face, estL[4*l+2]);
            if (error != MB_SUCCESS) return error;
          }

        //Cells
        if (estL[4*l+3])
          {
            error = mb->add_entities(hm_set[l], &level_mesh[l].start_cell, estL[4*l+3]);
            if (error != MB_SUCCESS) return error;
          }

        //debug
        std::cout<<"estL[0] ="<<estL[4*l]<<", estL[1] = "<<estL[4*l+1]<<", estL[2] = "<<estL[4*l+2]<<", estL[3] = "<<estL[4*l+3]<<std::endl;
        Range verts, ents;

        error = mb->get_entities_by_type(hm_set[l], MBVERTEX, verts);
        error = mb->get_entities_by_type(hm_set[l], MBEDGE, ents);

        if (verts.empty() || ents.empty())
          std::cout<<"something not right"<<std::endl;

        for (Range::iterator v = verts.begin(); v!= verts.end(); ++v)
          std::cout<<"v["<<(*v-*verts.begin())<<"] = "<<*v<<std::endl;

        for (Range::iterator e = ents.begin(); e!= ents.end(); ++e)
          std::cout<<"ent["<<(*e-*ents.begin())<<"] = "<<*e<<std::endl;
      }

    return MB_SUCCESS;
  }


  ErrorCode NestedRefine::generate_hm(int *level_degrees, int num_level, int *hmest, EntityHandle *hm_set)
  {
    ErrorCode error;

    for (int l = 0; l<num_level; l++)
      {
        //Create arrays for storing the current level
        error = create_hm_storage_single_level(&hm_set[l], l, hmest);
        if (error != MB_SUCCESS) return error;

        //Copy the old vertices along with their coordinates
        error = copy_vertices_from_prev_level(l);
        if (error != MB_SUCCESS) return error;

        //Create the new entities and new vertices
        error = construct_hm_entities(l, level_degrees[l]);
        if (error != MB_SUCCESS) return error;
      }
    return MB_SUCCESS;
  }


  ErrorCode NestedRefine::construct_hm_entities(int cur_level, int deg)
  {
    ErrorCode error;

    //Generate mesh for current level by refining previous level.
    if (thismeshtype == CURVE)
      {
        error = construct_hm_1D(cur_level, deg);
        if (MB_SUCCESS != error) return error;
      }
     else if(thismeshtype == SURFACE || thismeshtype == SURFACE_MIXED)
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
        nverts_prev = _verts.size();
        nents_prev = _edges.size();
      }

    int d = get_index_from_degree(deg);
    int vtotal = 2+ refTemplates[0][d].total_new_verts;
    EntityHandle *vbuffer = new EntityHandle[vtotal];

    std::vector<EntityHandle> conn(2);
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
          edge = *_edges.begin()+eid; // Makes the assumption initial mesh in contiguous in memory

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
              vbuffer[i] = level_mesh[cur_level].start_vertex + (conn[i] - *_verts.begin());

            std::cout<<"vbuffer["<<i<<"] = "<<vbuffer[i]<<std::endl;
          }

        //Adding rest of the entityhandles to working buffer for vertices.
        int num_new_verts = refTemplates[0][d].total_new_verts;

        for (int i=0; i<num_new_verts; i++)
          {
            vbuffer[i+2] = level_mesh[cur_level].start_vertex+count_verts;
            count_verts += 1;
            std::cout<<"vbuffer["<<i+2<<"] = "<<vbuffer[i+2]<<std::endl;
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
            std::cout<<"ent_buffer["<<i<<"] = "<<ent_buffer[i]<<std::endl;
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
                id1 = conn[1]-level_mesh[cur_level-1].start_vertex;
              }
            else
              {
                id1 = _verts.index(conn[0]);
                id2 = _verts.index(conn[1]);
              }

            level_mesh[cur_level].coordinates[0][idx] = (1-xi)*level_mesh[cur_level].coordinates[0][id1] + xi*level_mesh[cur_level].coordinates[0][id2];
            level_mesh[cur_level].coordinates[1][idx] = (1-xi)*level_mesh[cur_level].coordinates[1][id1] + xi*level_mesh[cur_level].coordinates[1][id2];
            level_mesh[cur_level].coordinates[2][idx] = (1-xi)*level_mesh[cur_level].coordinates[2][id1] + xi*level_mesh[cur_level].coordinates[2][id2];
          }

        delete [] ent_buffer;
      }

    error = update_global_ahf(MBEDGE, cur_level, deg);
    if (error != MB_SUCCESS) return error;

    error = print_tags_1D(cur_level);
    if (error != MB_SUCCESS) return error;

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
        nverts_prev = _verts.size();
        nents_prev = _faces.size();
      }

    //Create some book-keeping arrays over the old mesh to avoid introducing duplicate vertices
    int nepf = local_maps_2d(*_faces.begin());
    EntityType type = mb->type_from_handle(*(_faces.begin()));
    int findex = type-1;

    int d = get_index_from_degree(deg);
    int tnv = refTemplates[findex][d].total_new_verts;
    int vtotal = nepf + tnv;
    EntityHandle *vbuffer = new EntityHandle[vtotal];

    int nve = refTemplates[findex][d].nv_edge;
    int nvf = refTemplates[findex][d].nv_face;
    std::vector<EntityHandle> trackvertsF(nepf*nve*nents_prev, 0);
    std::vector<int> flag_ents(nents_prev, 0);

    int count_nverts = nverts_prev;
    int count_nents = 0;
    std::vector<EntityHandle> conn;

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int fid = 0; fid < nents_prev; fid++)
      {        
        flag_ents[fid] = 1;
        conn.clear();

        //EntityHandle of the working face
        EntityHandle face;
        if (cur_level)
          face = level_mesh[cur_level-1].start_face + fid;
        else
          face = *_faces.begin() + fid;

        error = get_connectivity(face, cur_level-1, conn);
        if (error != MB_SUCCESS) return error;

        // Add the corners to vbuffer
        for (int i=0; i<(int)conn.size(); i++)
          vbuffer[i] = level_mesh[cur_level].start_vertex+(conn[i]-level_mesh[cur_level-1].start_vertex);

        //Gather vertices already added to tracking array due to refinement of the sibling faces
        for (int i = 0; i < nepf; i++)
          {
            for (int j = 0; j < nve; j++)
              {
                int id = refTemplates[findex][d].vert_on_edges[i][j];
                vbuffer[id] = trackvertsF[fid*nve*nepf+nve*i+j];
              }
          }

        //Add the remaining vertex handles to vbuffer for the current level for the working face
        for (int i=0; i<nvf; i++)
          {
            if (~vbuffer[i+conn.size()]){
                vbuffer[i+conn.size()] = level_mesh[cur_level].start_vertex+count_nverts;
                count_nverts += 1;
              }
          }

        //Use the template to obtain the subentities
        int idx;
        int etotal = refTemplates[findex][d].total_new_ents;
        EntityHandle *ent_buffer = new EntityHandle[etotal];

        for (int i = 0; i < etotal; i++)
          {
            for (int k = 0; k < nepf; k++)
              {
                idx = refTemplates[findex][d].ents_conn[i][k];
                level_mesh[cur_level].face_conn[nepf*count_nents+k] = vbuffer[idx];
                ent_buffer[i] = count_nents;
                count_nents += 1;
              }
          }

        // Update the local AHF maps
        error = update_local_ahf(deg, type, vbuffer, ent_buffer, etotal);
        if (error != MB_SUCCESS) return error;

        //Add the new vertices to the tracking array
        std::vector<EntityHandle> sib_fids;
        std::vector<int> sib_leids;
        std::vector<int> sib_orient;
        int id;
        for (int i = 0; i < nepf; i++)
          {
            for (int j = 0; j < nve; j++)
              {
                // Add the vertices to trackvertsF for fid
                id = refTemplates[findex][d].vert_on_edges[i][j];
                trackvertsF[fid*nve*nepf+nve*i+j] = vbuffer[id];

                error = get_up_adjacencies_2d(face, i, false, sib_fids, true,  &sib_leids, true, &sib_orient);
                if (error != MB_SUCCESS) return error;

                //Add the vertices to trackvertsF for siblings of fid
                for (int s = 0; s < (int)sib_fids.size(); s++)
                  {
                    if (sib_orient[s]) // Same half-edge direction as the current half-edge
                      {
                        id = refTemplates[findex][d].vert_on_edges[sib_leids[s]][j];
                        trackvertsF[sib_fids[s]*nve*nepf+nve*sib_leids[s]+j] = vbuffer[id];
                      }
                    else
                      {
                        id = refTemplates[findex][d].vert_on_edges[sib_leids[s]][nve-j];
                        trackvertsF[sib_fids[s]*nve*nepf+nve*sib_leids[s]+j] = vbuffer[id];
                      }
                  }
              }
          }
        // Update the global maps
        error = update_global_ahf(type, cur_level, deg);
        if (error != MB_SUCCESS) return error;

        //Compute the coordinates of the new vertices
        double *corner_coords = new double[nepf*3];
        error = get_coordinates(conn, nepf, cur_level, corner_coords);
        if (error != MB_SUCCESS) return error;

        error = compute_coordinates(cur_level, deg, type, vbuffer, vtotal, corner_coords);
        if (error != MB_SUCCESS) return error;

        delete [] ent_buffer;
        delete [] corner_coords;
      }

    delete [] vbuffer;

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
        nverts_prev = _verts.size();
        nents_prev = _cells.size();
      }

    //Create some book-keeping arrays over the parent mesh to avoid introducing duplicate vertices
    EntityType type = mb->type_from_handle(*(_cells.begin()));
    int cindex = type -1;
    int d = get_index_from_degree(deg);
    int ne = refTemplates[cindex][d].nv_edge;
    int nvf = refTemplates[cindex][d].nv_face;
    int nvtotal = refTemplates[cindex][d].total_new_verts;

    int index = get_index_from_type(*(_cells.begin()));
    int nvpc = lConnMap3D[index].num_verts_in_cell;
    int nepc = lConnMap3D[index].num_edges_in_cell;
    int nfpc = lConnMap3D[index].num_faces_in_cell;

    int vtotal = nvpc + nvtotal;
   // std::vector<EntityHandle> vbuffer(vtotal,0);
    EntityHandle *vbuffer = new EntityHandle[vtotal];

    std::vector<EntityHandle> trackvertsC_edg(nepc*ne*nents_prev, 0);
    std::vector<EntityHandle> trackvertsC_face(nfpc*nvf*nents_prev, 0);
    std::vector<int> flag_ents(nents_prev, 0);

    int count_nverts = nverts_prev;
    int count_ents = 0;
    std::vector<EntityHandle> conn;

    //Step 1: Create the subentities via refinement of the previous mesh
    for (int cid = 0; cid < nents_prev; cid++)
      {
        flag_ents[cid] = 1;
        conn.clear();

        //EntityHandle of the working cell
        EntityHandle cell;
        if (cur_level)
          cell = level_mesh[cur_level-1].start_cell + cid;
        else
          cell = _cells[cid];

        error = get_connectivity(cell, cur_level, conn);
        if (error != MB_SUCCESS) return error;

        // Add the corners to vbuffer
        for (int i=0; i<(int)conn.size(); i++)
          {
            if (cur_level)
              vbuffer[i] = (conn[i]-level_mesh[cur_level-1].start_vertex) + level_mesh[cur_level].start_vertex;
            else
              vbuffer[i] = (conn[i]-*_verts.begin()) + level_mesh[0].start_vertex;
          }

        //Add already created vertices to vbuffer
        int idx;
        for (int i=0; i<nepc; i++){
            for (int j=0; j<ne; j++)
              {
                idx = refTemplates[cur_level][d].vert_on_edges[i][j];
                vbuffer[idx] = trackvertsC_edg[cid*nepc*ne+ne*i+j];
              };
        };

        for (int i=0; i<nfpc; i++){
            for (int j=0; j<nvf; j++)
              {
                idx = refTemplates[cur_level][d].vert_on_faces[i][j];
                vbuffer[idx] = trackvertsC_face[cid*nfpc*nvf+nvf*i+j];
              };
          };

        //Add the remaining vertex handles to vbuffer for the current level for the working cell
        for (int i=0; i<nvtotal; i++){
            if (~vbuffer[i+conn.size()]){
                vbuffer[i+conn.size()] = level_mesh[cur_level].start_vertex+count_nverts;
                count_nverts += 1;
              };
          };

        //Use the template to obtain the subentities
        error = subdivide_cells(type, conn, cur_level, deg, vbuffer, &count_ents);
        if (error != MB_SUCCESS) return error;

        //Update tracking information
        error = update_tracking_verts(cell, cur_level, deg, trackvertsC_edg, trackvertsC_face, vbuffer);
        if (error != MB_SUCCESS) return error;
      }

    // Update the global maps
    error = update_global_ahf(type, cur_level, deg);
    if (error != MB_SUCCESS) return error;

    return MB_SUCCESS;
  }

  ErrorCode NestedRefine::subdivide_cells(EntityType type, std::vector<EntityHandle> conn, int cur_level, int deg, EntityHandle *vbuffer,  int *count_ents)
  {
    //Subdivide the cell into children cells. Also update the sibling half-facet maps for the new cells.
    //Obtain incident edges/faces on this cell
    ErrorCode error;
    int idx;
    if (type == MBTET)
      {
        // Subdivide the tet
       // error = subdivide_tets(conn, cur_level, deg, vbuffer, count_ents);
        //if (error != MB_SUCCESS) return error;
      }
    else {
        // Connectivity
        int index = get_index_from_type(*(_cells.begin()));
        int nvpc = lConnMap3D[index].num_verts_in_cell;
        int d = get_index_from_degree(deg);
        int nents_total = refTemplates[type-1][d].total_new_ents;
        EntityHandle  *ent_buffer = new EntityHandle[nents_total];
        for (int i = 0; i < nents_total; i++)
          {

            for (int k = 0; k < nvpc; k++)
              {
                idx = refTemplates[type-1][d].ents_conn[i][k];
                level_mesh[cur_level].cell_conn[nvpc*count_ents[0]+k] = vbuffer[idx];
                ent_buffer[i] = count_ents[0];
                count_ents[0] += 1;
              }

          }

        //Update local ahf maps
        int vtotal = refTemplates[type-1][d].total_new_verts;
        error = update_local_ahf(deg, type, vbuffer, ent_buffer, nents_total);
        if (error != MB_SUCCESS) return error;

        //Coordinates
        double *corner_coords = new double[nvpc*3];
        error = get_coordinates(conn, nvpc, cur_level, corner_coords);
        if (error != MB_SUCCESS) return error;

     //   EntityType type = mb->type_from_handle(cell);
        error = compute_coordinates(cur_level, deg, type, vbuffer, vtotal, corner_coords);
        if (error != MB_SUCCESS) return error;

        delete [] ent_buffer;
        delete [] corner_coords;
      }

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

  ErrorCode NestedRefine::compute_coordinates(int cur_level, int deg, EntityType type, EntityHandle *vbuffer, int vtotal, double *corner_coords)
  {
    EntityHandle vstart = level_mesh[cur_level].start_vertex;
    double x, y, z;
    int d = get_index_from_degree(deg);
    if (type == MBTRI)
      {
        double xi, eta, N[3];
        int findex = mb->type_from_handle(*(_faces.begin()))-1;

        for (int i=0; i<vtotal; i++)
          {
            xi = refTemplates[findex][d].vert_nat_coord[i][0];
            eta = refTemplates[findex][d].vert_nat_coord[i][1];
            N[0] = 1-xi-eta; N[1] = xi; N[2] = eta;

            for (int j=0; j<3; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;
          }
      }
      else if (type == MBQUAD)
      {
        double xi, eta, N[4];
        int findex = mb->type_from_handle(*(_faces.begin()))-1;

        for (int i=0; i<vtotal; i++)
          {
            xi = refTemplates[findex][d].vert_nat_coord[i][0];
            eta = refTemplates[findex][d].vert_nat_coord[i][1];
            N[0] = (1-xi)*(1-eta); N[1] = xi*(1-eta); N[2] = xi*eta, N[3] = (1-xi)*eta;

            for (int j=0; j<4; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;
          }

      }
      else if (type == MBTET)
      {
        double xi, eta, zeta, N[4];
        int cindex = mb->type_from_handle(*(_cells.begin()))-1;

        for (int i=0; i<vtotal; i++)
          {
            xi = refTemplates[cindex][d].vert_nat_coord[i][0];
            eta = refTemplates[cindex][d].vert_nat_coord[i][1];
            zeta = refTemplates[cindex][d].vert_nat_coord[i][2];

            N[0] = 1-xi-eta-zeta; N[1] = xi; N[2] = eta, N[3] = zeta;

            for (int j=0; j<4; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;
          }

      }
      else if (type == MBPRISM)
      {
        double xi, eta, zeta, N[6];
        int cindex = mb->type_from_handle(*(_cells.begin()))-1;

        for (int i=0; i<vtotal; i++)
          {
            xi = refTemplates[cindex][d].vert_nat_coord[i][0];
            eta = refTemplates[cindex][d].vert_nat_coord[i][1];
            zeta = refTemplates[cindex][d].vert_nat_coord[i][2];

            N[0] = (1-xi-eta)*(1-zeta), N[1] = xi*(1-zeta), N[2] = eta*(1-zeta), N[3] = (1-xi-eta)*zeta, N[4] = xi*zeta, N[5] = eta*zeta;

            for (int j=0; j<6; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;
          }

      }
      else if (type == MBHEX)
      {
        double xi, eta, zeta, N[8];
        double t1, t2, t3, t12, t13, t23, s12, s13, s23;
        int cindex = mb->type_from_handle(*(_cells.begin()))-1;

        for (int i=0; i<vtotal; i++)
          {
            xi = refTemplates[cindex][d].vert_nat_coord[i][0];
            eta = refTemplates[cindex][d].vert_nat_coord[i][1];
            zeta = refTemplates[cindex][d].vert_nat_coord[i][2];

            t1 = 1-xi; t2 = 1- eta; t3 = 1-zeta;
            t12 = t1*t2; t13 = t1*t3; t23 = t2*t3;
            s12 = xi*eta; s13 = xi*zeta; s23 = eta*zeta;
            N[0] = t12*t3; N[1] = xi*t23; N[2] = s12*t3; N[3] = t13*eta; N[4] = t12*zeta; N[5] = s13*t2; N[6] = s12*zeta; N[7] = t1*s23;


            for (int j=0; j<8; j++)
              {
                x += N[j]*corner_coords[3*j];
                y += N[j]*corner_coords[3*j+1];
                z += N[j]*corner_coords[3*j+2];
              }

            level_mesh[cur_level].coordinates[0][vbuffer[i]-vstart] = x;
            level_mesh[cur_level].coordinates[1][vbuffer[i]-vstart] = y;
            level_mesh[cur_level].coordinates[2][vbuffer[i]-vstart] = z;
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
      nhf = local_maps_2d(*_faces.begin());
      nv = nhf;
      total_new_verts = refTemplates[type-1][d].total_new_verts;
    }
  else if (type == MBPRISM || type == MBHEX)
    {
      int index = get_index_from_type(*_cells.begin());
      nhf = lConnMap3D[index].num_faces_in_cell;
      nv =  lConnMap3D[index].num_verts_in_cell;
      total_new_verts = refTemplates[type-1][d].total_new_verts;
    }

  //Update the vertex to half-facet map
  for (int i=0; i<total_new_verts; i++)
    {
      EntityHandle ent=0;
      int lid=0;

      std::cout<<"vbuffer["<<i+2<<"] = "<<vbuffer[i+nv]<<std::endl;

      error = get_incident_tag(type, vbuffer[i+nv], &ent, &lid);
      if (error != MB_SUCCESS) return error;

      if (ent)
        continue;

      int id = refTemplates[type-1][d].v2hf[i+nv][0]-1;
      ent = ent_buffer[id];
      lid = refTemplates[type-1][d].v2hf[i+nv][1];

      std::cout<<"ent = "<<ent<<", lid = "<<lid<<std::endl;
      error = set_incident_tag(type, vbuffer[i+2], &ent, &lid);
      if (error != MB_SUCCESS) return error;

      //Check if value set is correct
      EntityHandle get_ent;
      int get_lid;
      error = get_incident_tag(type, vbuffer[i+nv], &get_ent, &get_lid);
      if (error != MB_SUCCESS) return error;
      std::cout<<"get_ent = "<<get_ent<<", get_lid = "<<get_lid<<std::endl;
    }

  //Update the sibling half-facet map
  for (int i=0; i< etotal; i++)
    {
      EntityHandle  *sib_entids = new EntityHandle[nhf];
      int *sib_lids = new int[nhf];

      error = get_sibling_tag(type, ent_buffer[i], sib_entids, sib_lids);
      if (error != MB_SUCCESS) return error;

      for (int l=0; l< nhf; l++)
        {
          std::cout<<"sib_entids["<<l<<"]= "<<sib_entids[l]<<std::endl;
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
          std::cout<<"sib_entids["<<l<<"]= "<<sib_entids[l]<<std::endl;
        }

      error = set_sibling_tag(type, ent_buffer[i], sib_entids, sib_lids);
      if  (error != MB_SUCCESS) return error;

      EntityHandle *set_entids = new EntityHandle[nhf];
      int *set_lids = new int[nhf];

      for (int l=0; l< nhf; l++)
        {
          if (sib_entids[l]){
              error = get_sibling_tag(type, sib_entids[l], set_entids, set_lids);
              if (error != MB_SUCCESS) return error;

              set_entids[sib_lids[l]] = ent_buffer[i];
              set_lids[sib_lids[l]] = l;

              error = set_sibling_tag(type, sib_entids[l], set_entids, set_lids);
              if (error != MB_SUCCESS) return error;
            }
        }

      delete [] sib_entids;
      delete [] sib_lids;
      delete [] set_entids;
      delete [] set_lids;
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
       nverts_prev = _verts.size();
       nents_prev = _edges.size();
     }

   //Update the vertex to half-facet maps for duplicate vertices
   for (int i=0; i<nverts_prev; i++)
     {
       //Vertex id in the previous mesh
       EntityHandle vid;
       if (cur_level)
         vid = level_mesh[cur_level-1].start_vertex + i;
       else
         vid = _verts[i];

       //Get the incident half-vert in the previous mesh
       EntityHandle inci_ent;
       int inci_lid;

       error = get_incident_tag(MBEDGE, vid, &inci_ent, &inci_lid);
       if (error != MB_SUCCESS) return error;

     //  if (!inci_ent)
     //    continue;

       // Obtain the corresponding incident child in the current mesh
       int lvid = get_local_vid(vid, inci_ent, cur_level-1);
       int chid = refTemplates[0][d].v2hf[lvid][0]-1;

       int pid;
       if (cur_level)
         pid = inci_ent - level_mesh[cur_level].start_edge;
       else
         pid = inci_ent - *_edges.begin();

       int ind = nchilds*pid;

       EntityHandle child_ent = level_mesh[cur_level].start_edge + ind+chid ;
       int child_lid = refTemplates[0][d].v2hf[lvid][1];

       EntityHandle cur_vid = level_mesh[cur_level].start_vertex + i;
       error = set_incident_tag(MBEDGE, cur_vid, &child_ent, &child_lid);
       if (error != MB_SUCCESS) return error;
     }

   //Update the sibling half-facet maps across entities
   for (int i=0; i< nents_prev; i++)
     {
       EntityHandle ent;
       if (cur_level)
         ent = level_mesh[cur_level-1].start_edge + i;
       else
         ent = _edges[i];

       EntityHandle *sib_entids = new EntityHandle[nhf];
       int *sib_lids = new int[nhf];

       error = get_sibling_tag(MBEDGE, ent, sib_entids, sib_lids);
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

           error = get_sibling_tag(MBEDGE, child_ent, sib_childs, sib_chlids);
           if (error != MB_SUCCESS) return error;

           //If the sibling already exists, dont do anything
           if (sib_childs[ch_lid])
             continue;

           //Get the correponding child of the sibling of the current parent
           int psib;
           if (cur_level)
             psib = sib_entids[l] - level_mesh[cur_level].start_edge;
           else
             psib = sib_entids[l] - *_edges.begin();

           int plid = sib_lids[l];

           id = refTemplates[0][d].ents_on_pent[plid][1]-1;
           idx = nchilds*psib;

           EntityHandle psib_child = level_mesh[cur_level].start_edge + idx+id ;
           int psib_chlid = plid;

           //Set the siblings
           sib_childs[ch_lid] = psib_child;
           sib_chlids[ch_lid] = psib_chlid;

           error = set_sibling_tag(MBEDGE, child_ent, sib_childs, sib_chlids);
           if (error != MB_SUCCESS) return error;
         }
     }

   return MB_SUCCESS;
 }

 ErrorCode NestedRefine::update_global_ahf_2D(int cur_level, int deg)
 {
   ErrorCode error;

   EntityType type = mb->type_from_handle(*_faces.begin());
   int nhf, nchilds, nents_prev;

   nhf = local_maps_2d(*_faces.begin());
   int d = get_index_from_degree(deg);
   nchilds = refTemplates[type-1][d].total_new_ents;

   if (cur_level)
     nents_prev = level_mesh[cur_level-1].num_faces;
   else
     nents_prev = _faces.size();



   //Loop over all entities and the half-facets of the previous mesh
   for (int i=0; i< nents_prev; i++)
     {
       EntityHandle ent;
       if (cur_level)
         ent = level_mesh[cur_level-1].start_face + i;
       else
         ent = _faces[i];

       std::vector<EntityHandle> conn;
       error = get_connectivity(ent, cur_level, conn);
       if (error != MB_SUCCESS) return error;

       EntityHandle *sib_entids = new EntityHandle[nhf];
       int *sib_lids = new int[nhf];

       error = get_sibling_tag(type, ent, sib_entids, sib_lids);
       if (error != MB_SUCCESS) return error;

       int nc, id, idx;

       for (int l=0; l < nhf; l++)
         {
           if (!sib_entids[l])
             continue;


           //Find the children incident on the half-facet
           nc = refTemplates[type][d].ents_on_pent[l][0];

           for (int k=0; k<nc; k++)
             {
               id = refTemplates[type][d].ents_on_pent[l][k+1];
               idx = nchilds*i;
               EntityHandle child_ent = level_mesh[cur_level].start_face + idx+id ;
               int ch_lid = l;

               //Find the sibling of the child
               EntityHandle *sib_childs= new EntityHandle[nhf];
               int *sib_chlids = new int[nhf];

               error = get_sibling_tag(type, child_ent, sib_childs, sib_chlids);
               if (error != MB_SUCCESS) return error;

               //If the sibling already exists, dont do anything
               if (sib_childs[ch_lid])
                 continue;

               //Get the correponding child of the sibling of the current parent
               int psib;
               if (cur_level)
                 psib = sib_entids[l] - level_mesh[cur_level-1].start_face;
               else
                 psib = sib_entids[l] - *_faces.begin();

               int plid = sib_lids[l];

               id = refTemplates[type][d].ents_on_pent[plid][2*k];
               idx = nchilds*psib;

               EntityHandle psib_child = level_mesh[cur_level].start_face + idx+id ;
               int psib_chlid = refTemplates[type][d].ents_on_pent[plid][k+1];

               //Set the siblings
               sib_childs[ch_lid] = psib_child;
               sib_chlids[ch_lid] = psib_chlid;

               error = set_sibling_tag(MBEDGE, child_ent, sib_childs, sib_chlids);
               if (error != MB_SUCCESS) return error;

               delete [] sib_childs;
               delete [] sib_chlids;
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
  int nhf, nchilds, nents_prev;
  EntityType type = mb->type_from_handle(*_cells.begin());
  int index = get_index_from_type(*_cells.begin());
  nhf = lConnMap3D[index].num_faces_in_cell;

  int d = get_index_from_degree(deg);
  nchilds = refTemplates[type-1][d].total_new_ents;
  if (cur_level)
    nents_prev = level_mesh[cur_level-1].num_edges;
  else
    nents_prev = _edges.size();

  //Loop over all entities and the half-facets of the previous mesh
  for (int i=0; i< nents_prev; i++)
    {
      EntityHandle ent;
      if (cur_level)
        ent = level_mesh[cur_level-1].start_edge + i;
      else
        ent = _edges[i];

      EntityHandle *sib_entids = new EntityHandle[nhf];
      int *sib_lids = new int[nhf];

      error = get_sibling_tag(MBEDGE, ent, sib_entids, sib_lids);
      if (error != MB_SUCCESS) return error;

      int id, idx;

      for (int l=0; l < nhf; l++)
        {
          //Find the child incident on the half-facet
          id = refTemplates[0][d].ents_on_pent[l][0];
          idx = nchilds*i;
          EntityHandle child_ent = level_mesh[cur_level].start_edge + idx+id ;
          int ch_lid = refTemplates[0][d].ents_on_pent[l][1];

          //Find the sibling of the child
          EntityHandle *sib_childs= new EntityHandle[nhf];
          int *sib_chlids = new int[nhf];

          error = get_sibling_tag(MBEDGE, child_ent, sib_childs, sib_chlids);
          if (error != MB_SUCCESS) return error;

          //If the sibling already exists, dont do anything
          if (sib_childs[ch_lid])
            continue;

          //Get the correponding child of the sibling of the current parent
          int psib;
          if (cur_level)
            psib = sib_entids[l] - level_mesh[cur_level].start_edge;
          else
            psib = sib_entids[l] - *_edges.begin();

          int plid = sib_lids[l];

          id = refTemplates[0][d].ents_on_pent[plid][0];
          idx = nchilds*psib;

          EntityHandle psib_child = level_mesh[cur_level].start_edge + idx+id ;
          int psib_chlid = refTemplates[0][d].ents_on_pent[plid][1];

          //Set the siblings
          sib_childs[ch_lid] = psib_child;
          sib_chlids[ch_lid] = psib_chlid;

          error = set_sibling_tag(MBEDGE, child_ent, sib_childs, sib_chlids);
          if (error != MB_SUCCESS) return error;
        }
    }

  return MB_SUCCESS;
}

/*****************************
 *    Tag Management            *
******************************/

ErrorCode NestedRefine::get_sibling_tag(EntityType type, EntityHandle ent,  EntityHandle *sib_entids, int *sib_lids)
{
  ErrorCode error;

  if (type == MBEDGE)
    {
      error = mb->tag_get_data(sibhvs_eid, &ent, 1, sib_entids);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_get_data(sibhvs_lvid, &ent, 1, sib_lids);
      if (MB_SUCCESS != error) return error;
    }
  else if (type == MBTRI || type == MBQUAD)
    {
      error = mb->tag_get_data(sibhes_fid, &ent, 1, sib_entids);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_get_data(sibhes_leid, &ent, 1, sib_lids);
      if (MB_SUCCESS != error) return error;
    }
  else
    {
      error = mb->tag_get_data(sibhfs_cid, &ent, 1, sib_entids);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_get_data(sibhfs_lfid, &ent, 1, sib_lids);
      if (MB_SUCCESS != error) return error;
    }
  return MB_SUCCESS;
}

ErrorCode NestedRefine::set_sibling_tag(EntityType type, EntityHandle ent, EntityHandle *set_entids, int *set_lids)
{
  ErrorCode error;
  if (type == MBEDGE)
    {
      error = mb->tag_set_data(sibhvs_eid, &ent, 1, set_entids);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_set_data(sibhvs_lvid, &ent, 1, set_lids);
      if (MB_SUCCESS != error) return error;
    }
  else if (type == MBTRI || type == MBQUAD)
    {
      error = mb->tag_set_data(sibhes_fid, &ent, 1, set_entids);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_set_data(sibhes_leid, &ent, 1, set_lids);
      if (MB_SUCCESS != error) return error;
    }
  else
    {
      error = mb->tag_set_data(sibhfs_cid, &ent, 1, set_entids);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_set_data(sibhfs_lfid, &ent, 1, set_lids);
      if (MB_SUCCESS != error) return error;
    }
  return MB_SUCCESS;
}

ErrorCode NestedRefine::get_incident_tag(EntityType type, EntityHandle vid, EntityHandle *inci_entid, int  *inci_lid)
{
  ErrorCode error;

  if (type == MBEDGE)
    {
      error = mb->tag_get_data(v2hv_eid, &vid, 1, inci_entid);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_get_data(v2hv_lvid, &vid, 1, inci_lid);
      if (MB_SUCCESS != error) return error;
    }
  else if (type == MBTRI || type == MBQUAD)
    {
      error = mb->tag_get_data(v2he_fid, &vid, 1, inci_entid);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_get_data(v2he_leid, &vid, 1, inci_lid);
      if (MB_SUCCESS != error) return error;
    }
  else
    {
      error = mb->tag_get_data(v2hf_cid, &vid, 1, inci_entid);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_get_data(v2hf_lfid, &vid, 1, inci_lid);
      if (MB_SUCCESS != error) return error;
    }
  return MB_SUCCESS;
}

ErrorCode NestedRefine::set_incident_tag(EntityType type, EntityHandle vid, EntityHandle *set_entid, int *set_lid)
{
  ErrorCode error;
  if (type == MBEDGE)
    {
      error = mb->tag_set_data(v2hv_eid, &vid, 1, set_entid);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_set_data(v2hv_lvid, &vid, 1, set_lid);
      if (MB_SUCCESS != error) return error;
    }
  else if (type == MBTRI || type == MBQUAD)
    {
      error = mb->tag_set_data(v2he_fid, &vid, 1,&set_entid);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_set_data(v2he_leid, &vid, 1, &set_lid);
      if (MB_SUCCESS != error) return error;
    }
  else
    {
      error = mb->tag_set_data(v2hf_cid, &vid, 1,&set_entid);
      if (MB_SUCCESS != error) return error;

      error = mb->tag_set_data(v2hf_lfid, &vid, 1, &set_lid);
      if (MB_SUCCESS != error) return error;
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
      int nverts_in = _verts.size();
      double* x_coords = new double[nverts_in] ;
      double* y_coords = new double[nverts_in] ;
      double* z_coords = new double[nverts_in] ;

      error = mb->get_coords(_verts, x_coords, y_coords, z_coords);
      if (error != MB_SUCCESS) return error;

      for (int i = 0; i < nverts_in; i++)
        {
          std::cout<<"x= "<<x_coords[i]<<", y = "<<y_coords[i]<<", z = "<<z_coords[i]<<std::endl;
          level_mesh[cur_level].coordinates[0][i] = x_coords[i];
          level_mesh[cur_level].coordinates[1][i] = y_coords[i];
          level_mesh[cur_level].coordinates[2][i] = z_coords[i];
        }

      delete [] x_coords;
      delete [] y_coords;
      delete [] z_coords;

    }
  return MB_SUCCESS;
  //To add: Map from old vertices to new duplicates: NOT NEEDED
}

ErrorCode NestedRefine::update_tracking_verts(EntityHandle cid, int cur_level, int deg, std::vector<EntityHandle> trackvertsC_edg, std::vector<EntityHandle> trackvertsC_face, EntityHandle *vbuffer)
{
  ErrorCode error;
  //There are two steps here.
  //First we need to add the vertices in vbuffer to the tracking arrays of edges and faces for the current cell.
  //Secondly, we need to do the same for the sibling half-faces.
  EntityType cell_type = mb->type_from_handle(*(_cells.begin()));
  int cindex = cell_type -1;
  int d = get_index_from_degree(deg);

  int nve = refTemplates[cindex][d].nv_edge;
  int nvf = refTemplates[cindex][d].nv_face;

  int index = get_index_from_type(*(_cells.begin()));
  int nepc = lConnMap3D[index].num_edges_in_cell;
  int nfpc = lConnMap3D[index].num_faces_in_cell;

  std::vector<EntityHandle> inc_cids;
  std::vector<int> inc_leids, inc_orient;
  int id;

  for (int i=0; i<nepc; i++)
    {
      error = get_up_adjacencies_edg_3d(cid, i, inc_cids, true, &inc_leids, true, &inc_orient);
      if (error != MB_SUCCESS) return error;

      for (int j=0; j< nve; j++)
        {
          //Add the vertices to edges of the current cell
          id = refTemplates[cindex][d].vert_on_edges[i][j];
          trackvertsC_edg[cid*nve*nepc+nve*i+j] = vbuffer[id];

          //Add the vertices to the edges of the incident cells
          for (int k=0; k<(int)inc_cids.size(); k++)
            {
              if (inc_orient[k]) // Same edge direction as the current edge
                {
                  id = refTemplates[cindex][d].vert_on_edges[inc_leids[k]][j];
                  trackvertsC_edg[inc_cids[k]*nve*nepc+nve*inc_leids[k]+j] = vbuffer[id];
                }
              else
                {
                  id = refTemplates[cindex][d].vert_on_edges[inc_leids[k]][nve-j];
                  trackvertsC_edg[inc_cids[k]*nve*nepc+nve*inc_leids[k]+j] = vbuffer[id];
                }
            }
        }
    }

  if (nvf) {
      std::vector<EntityHandle> sib_cids;
      std::vector<int>  sib_lfids;
      for (int i=0; i< nfpc; i++)
        {
          error = get_up_adjacencies_face_3d(cid, i, sib_cids, true, &sib_lfids);
          if (error != MB_SUCCESS) return error;

          if (sib_cids[0]) {
              int *id_sib = new int[nvf];
              error = match_and_reorder_vertices(cell_type, cur_level, deg, cid, i, sib_cids[0], sib_lfids[0], id_sib);

              for (int j=0; j< nvf; j++)
                {
                  //Add vertices to the tracking array of vertices on faces for the current cell
                  id = refTemplates[cindex][d].vert_on_faces[i][j];
                  trackvertsC_face[cid*nvf*nfpc+nvf*i+j] = vbuffer[id];

                  if (sib_cids[0]){
                      //Add vertices to the tracking array of vertices on faces for the sibling cell of the current cell
                      trackvertsC_face[sib_cids[0]*nvf*nfpc+nvf*sib_lfids[0]+j] = vbuffer[id_sib[j]];
                    }
                }
            }
        }
    }
  return MB_SUCCESS;
}

ErrorCode NestedRefine::match_and_reorder_vertices(EntityType type, int cur_level, int deg, EntityHandle cell, int lfid, EntityHandle sib_cell, int sib_lfid, int *id_sib)
{
  assert(deg ==2 || deg == 3);

  ErrorCode error;
  int index = get_index_from_type(*_cells.begin());
  int nvF = lConnMap3D[index].hf2v_num[lfid];
  int nco = permute_matrix[nvF-3].num_comb;
  int d = get_index_from_degree(deg);

  if (deg == 3)
    {
      std::vector<EntityHandle> cell_conn, sib_cell_conn;
      error = get_connectivity(cell, cur_level, cell_conn);
      if (error != MB_SUCCESS) return error;

      error = get_connectivity(sib_cell, cur_level, sib_cell_conn);
      if (error != MB_SUCCESS) return error;

      int *lface = new int[nvF];
      int *lface_sib = new int[nvF];
      for (int i=0; i<nvF; i++)
        {
          int id = lConnMap3D[index].hf2v[lfid][i];
          lface[i] = cell_conn[id];

          int idx = lConnMap3D[index].hf2v[sib_lfid][i];
          lface_sib[i] = sib_cell_conn[idx];
        }

      for (int i=0; i<nco; i++)
        {
          int count = 0;
          for (int j=0; j<nvF; j++)
            {
              int k = permute_matrix[nvF-3].mat[i][j];
              if (lface[j] == lface_sib[k])
                count += 1;
            }
          if (count == nvF)
            {
              int *id_buffer = new int[nvF];
              for (int j=0; j< nvF; j++)
                {
                  id_buffer[j] = refTemplates[type-1][d].vert_on_faces[sib_lfid][j];
                }
              for (int j=0; i<nvF; j++)
                {
                  int k =  permute_matrix[nvF-3].mat[i][j];
                  id_sib[j] = id_buffer[k];
                }
              break;
            }
        }
    }
  else
    {
      if (nvF == 3)
        return MB_SUCCESS;
      else
        id_sib[0] = refTemplates[type-1][d].vert_on_faces[sib_lfid][0];
    }
  return MB_SUCCESS;
}


int NestedRefine::find_shortest_diagonal_octahedron(double *coords)
{
    int diag;
    int diag_map[6] = {0,5,1,3,2,4};
    double len_sqr;
    double length = std::numeric_limits<double>::max();

    double x,y,z;
    for (int d=0; d<3; d++)
      {
        int id1 = diag_map[2*d];
        int id2 = diag_map[2*d+1];
        x = coords[3*id1] - coords[3*id2];
        y = coords[3*id1+1] - coords[3*id2+1];
        z = coords[3*id1+2] - coords[3*id2+2];
        len_sqr =  x*x+y*y+z*z;
        if ( sqrt(len_sqr) < length)
          {
            length = sqrt(len_sqr);
            diag = d;
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
  if (error != MB_SUCCESS) return error;

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
      EntityHandle eid; int lvid;
      EntityHandle vid = level_mesh[level].start_vertex+i;
      error = mb->tag_get_data(v2hv_eid, &vid, 1, &eid);
      if (MB_SUCCESS != error) return error;
      error = mb->tag_get_data(v2hv_lvid, &vid, 1, &lvid);
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
      error = mb->tag_get_data(sibhvs_eid, &ent, 1, eid);
      if (MB_SUCCESS != error) return error;
      error = mb->tag_get_data(sibhvs_lvid, &ent, 1, lvid);
      if (MB_SUCCESS != error) return error;

      std::cout<<"<"<<eid[0]<<","<<lvid[0]<<">"<<"      "<<"<"<<eid[1]<<","<<lvid[1]<<">"<<std::endl;
    }

  return MB_SUCCESS;
}

}//namesapce moab

