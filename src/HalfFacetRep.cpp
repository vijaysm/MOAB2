/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include "moab/HalfFacetRep.hpp"
#include <iostream>
#include <assert.h>
#include <vector>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"

namespace moab {

  const int MAXSIZE = 500;

  HalfFacetRep::HalfFacetRep(Core *impl)
  {
    assert(NULL != impl);
    mb = impl;
    mInitAHFmaps = false;
  }

  HalfFacetRep::~HalfFacetRep()
  {
    ErrorCode result;
    result = deinitialize();
    //if (MB_SUCCESS != result) return result;
  }

  MESHTYPE HalfFacetRep::get_mesh_type(int nverts, int nedges, int nfaces, int ncells)
  {
    MESHTYPE mesh_type;

    if (nverts && nedges && (!nfaces) && (!ncells))
      mesh_type = CURVE;
    else if (nverts && !nedges && nfaces && !ncells)
      mesh_type = SURFACE;
    else if (nverts && nedges && nfaces && !ncells)
      mesh_type = SURFACE_MIXED;
    else if (nverts && !nedges && !nfaces && ncells)
      mesh_type = VOLUME;
    else if (nverts && nedges && !nfaces && ncells)
      mesh_type = VOLUME_MIXED_1;
    else if (nverts && !nedges && nfaces && ncells)
      mesh_type = VOLUME_MIXED_2;
    else if (nverts && nedges && nfaces && ncells)
      mesh_type = VOLUME_MIXED;

    return mesh_type;
  }

  /*******************************************************
   * initialize                                          *
   ******************************************************/

  ErrorCode HalfFacetRep::initialize()
  {
    mInitAHFmaps = true;

    ErrorCode error;

    error = mb->get_entities_by_dimension( 0, 0, _verts);
    if (MB_SUCCESS != error) return error;

    error = mb->get_entities_by_dimension( 0, 1, _edges);
    if (MB_SUCCESS != error) return error;

    error = mb->get_entities_by_dimension( 0, 2, _faces);
    if (MB_SUCCESS != error) return error;
    
    error = mb->get_entities_by_dimension( 0, 3, _cells);
    if (MB_SUCCESS != error) return error;
    
    int nverts = _verts.size();
    int nedges = _edges.size();
    int nfaces = _faces.size();
    int ncells = _cells.size();

    MESHTYPE mesh_type = get_mesh_type(nverts, nedges, nfaces, ncells);
    std::cout<<"MeshType = "<<mesh_type<<std::endl;
 
    if (mesh_type == CURVE){
        error = init_curve();
        if (MB_SUCCESS != error) return error;
      }
    
    else if (mesh_type == SURFACE){
      error = init_surface();
      if (MB_SUCCESS != error) return error;
      }

    else if (mesh_type == SURFACE_MIXED){
        error = init_curve();
        if (MB_SUCCESS != error) return error;
        error = init_surface();
        if (MB_SUCCESS != error) return error;
      }
    else if (mesh_type == VOLUME){
        error = init_volume();
        if (MB_SUCCESS != error) return error;
      }

    else if (mesh_type == VOLUME_MIXED_1){
        error = init_curve();
        if (MB_SUCCESS != error) return error;
        error = init_volume();
        if (MB_SUCCESS != error) return error;
    }

    else if (mesh_type == VOLUME_MIXED_2){
        error = init_surface();
        if (MB_SUCCESS != error) return error;
        error = init_volume();
        if (MB_SUCCESS != error) return error;
      }

    else if (mesh_type == VOLUME_MIXED){
        error = init_curve();
        if (MB_SUCCESS != error) return error;
        error = init_surface();
        if (MB_SUCCESS != error) return error;
        error = init_volume();
        if (MB_SUCCESS != error) return error;
      }

    return MB_SUCCESS;
  }

  ErrorCode HalfFacetRep::init_curve()
  {
    ErrorCode error;

    EntityHandle sdefval[2] = {0,0};  int sval[2] = {0,0};
    EntityHandle idefval = 0; int ival = 0;

    error = mb->tag_get_handle("__SIBHVS_EID", 2, MB_TYPE_HANDLE, sibhvs_eid, MB_TAG_DENSE | MB_TAG_CREAT, sdefval);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_handle("__SIBHVS_LVID", 2, MB_TYPE_INTEGER, sibhvs_lvid, MB_TAG_DENSE | MB_TAG_CREAT, sval);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_handle("__V2HV_EID", 1, MB_TYPE_HANDLE, v2hv_eid, MB_TAG_DENSE | MB_TAG_CREAT, &idefval);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_handle("__V2HV_LVID", 1, MB_TYPE_INTEGER, v2hv_lvid, MB_TAG_DENSE | MB_TAG_CREAT, &ival);
    if (MB_SUCCESS != error) return error;

    error = determine_sibling_halfverts(_edges);
    if (MB_SUCCESS != error) return error;
    error = determine_incident_halfverts(_edges);
    if (MB_SUCCESS != error) return error;

    return MB_SUCCESS;
  }

  ErrorCode HalfFacetRep::init_surface()
  {
    ErrorCode error;

    int nepf = local_maps_2d(*_faces.begin());
    EntityHandle * sdefval = new EntityHandle[nepf];
    int * sval = new int[nepf];
    EntityHandle  idefval = 0;
    int ival = 0;
    for (int i=0; i<nepf; i++){
        sdefval[i]=0;
        sval[i]=0;
      }

    // Create tag handles
    error = mb->tag_get_handle("__SIBHES_FID", nepf, MB_TYPE_HANDLE, sibhes_fid, MB_TAG_DENSE | MB_TAG_CREAT, sdefval);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_handle("__SIBHES_LEID", nepf, MB_TYPE_INTEGER, sibhes_leid, MB_TAG_DENSE | MB_TAG_CREAT, sval);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_handle("__V2HE_FID", 1, MB_TYPE_HANDLE, v2he_fid, MB_TAG_DENSE | MB_TAG_CREAT, &idefval);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_handle("__V2HE_LEID", 1, MB_TYPE_INTEGER, v2he_leid, MB_TAG_DENSE | MB_TAG_CREAT, &ival);
    if (MB_SUCCESS != error) return error;

    // Construct ahf maps
    error = determine_sibling_halfedges(_faces);
    if (MB_SUCCESS != error) return error;
    error = determine_incident_halfedges(_faces);
    if (MB_SUCCESS != error) return error;

    delete [] sdefval;
    delete [] sval;

    return MB_SUCCESS;
  }

  ErrorCode HalfFacetRep::init_volume()
  {
    ErrorCode error;

    int index = get_index_from_type(*_cells.begin());
    int nfpc = lConnMap3D[index].num_faces_in_cell;
    EntityHandle * sdefval = new EntityHandle[nfpc];
    int * sval = new int[nfpc];
    EntityHandle idefval = 0;
    int ival = 0;
    for (int i = 0; i < nfpc; i++){
        sdefval[i] = 0;
        sval[i] = 0;
      }

    error = mb->tag_get_handle("__SIBHFS_CID", nfpc, MB_TYPE_HANDLE, sibhfs_cid, MB_TAG_DENSE | MB_TAG_CREAT, sdefval);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_handle("__SIBHFS_LFID", nfpc, MB_TYPE_INTEGER, sibhfs_lfid, MB_TAG_DENSE | MB_TAG_CREAT, sval);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_handle("__V2HF_CID", 1, MB_TYPE_HANDLE, v2hf_cid, MB_TAG_DENSE | MB_TAG_CREAT, &idefval);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_handle("__V2HF_LFID", 1, MB_TYPE_INTEGER, v2hf_lfid, MB_TAG_DENSE | MB_TAG_CREAT, &ival);
    if (MB_SUCCESS != error) return error;

    error = determine_sibling_halffaces(_cells);
    if (MB_SUCCESS != error) return error;
    error = determine_incident_halffaces(_cells);
    if (MB_SUCCESS != error) return error;

    delete [] sdefval;
    delete [] sval;

    return MB_SUCCESS;
  }

  /*******************************************************
   * deinitialize                                        *
   ******************************************************/
   ErrorCode HalfFacetRep::deinitialize()
   {
     ErrorCode error;
     MESHTYPE mesh_type = get_mesh_type(_verts.size(), _edges.size(), _faces.size(), _cells.size());

     if (mesh_type == CURVE){
         error = deinit_curve();
         if (MB_SUCCESS != error) return error;
       }
     else if (mesh_type == SURFACE){
         error = deinit_surface();
         if (MB_SUCCESS != error) return error;
       }
     else if (mesh_type == SURFACE_MIXED){
         error = deinit_curve();
         if (MB_SUCCESS != error) return error;
         error = deinit_surface();
         if (MB_SUCCESS != error) return error;
       }
     else if (mesh_type == VOLUME){
         error = deinit_volume();
         if (MB_SUCCESS != error) return error;
       }
     else if (mesh_type == VOLUME_MIXED_1){
         error = deinit_curve();
         if (MB_SUCCESS != error) return error;
         error = deinit_volume();
         if (MB_SUCCESS != error) return error;
       }
     else if (mesh_type == VOLUME_MIXED_2){
         error = deinit_surface();
         if (MB_SUCCESS != error) return error;
         error = deinit_volume();
         if (MB_SUCCESS != error) return error;         
       }
     else if (mesh_type == VOLUME_MIXED){
         error = deinit_curve();
         if (MB_SUCCESS != error) return error;
         error = deinit_surface();
         if (MB_SUCCESS != error) return error;
         error = deinit_volume();
         if (MB_SUCCESS != error) return error;
       }

     return MB_SUCCESS;
   }

   ErrorCode HalfFacetRep::deinit_curve(){
     ErrorCode error;
     error = mb->tag_delete(sibhvs_eid);
     if (MB_SUCCESS != error) return error;
     error = mb->tag_delete(sibhvs_lvid);
     if (MB_SUCCESS != error) return error;
     error = mb->tag_delete(v2hv_eid);
     if (MB_SUCCESS != error) return error;
     error = mb->tag_delete(v2hv_lvid);
     if (MB_SUCCESS != error) return error;

     return MB_SUCCESS;
   }

   ErrorCode HalfFacetRep::deinit_surface(){
     ErrorCode error;
     error = mb->tag_delete(sibhes_fid);
     if (MB_SUCCESS != error) return error;
     error = mb->tag_delete(sibhes_leid);
     if (MB_SUCCESS != error) return error;
     error = mb->tag_delete(v2he_fid);
     if (MB_SUCCESS != error) return error;
     error = mb->tag_delete(v2he_leid);
     if (MB_SUCCESS != error) return error;

     return MB_SUCCESS;
   }

   ErrorCode HalfFacetRep::deinit_volume(){
     ErrorCode error;
     error = mb->tag_delete(sibhfs_cid);
     if (MB_SUCCESS != error) return error;
     error = mb->tag_delete(sibhfs_lfid);
     if (MB_SUCCESS != error) return error;
     error = mb->tag_delete(v2hf_cid);
     if (MB_SUCCESS != error) return error;
     error = mb->tag_delete(v2hf_lfid);
     if (MB_SUCCESS != error) return error;

     return MB_SUCCESS;
   }

   //////////////////////////////////////////////////
   ErrorCode HalfFacetRep::print_tags()
   {
     ErrorCode error;
     int nepf = local_maps_2d(*_faces.begin());
     int index = get_index_from_type(*_cells.begin());
     int nfpc = lConnMap3D[index].num_faces_in_cell;

     //////////////////////////
     // Print out the tags
     EntityHandle start_edge = *_edges.begin();
     EntityHandle start_face = *_faces.begin();
     EntityHandle start_cell = *_cells.begin();
     std::cout<<"start_edge = "<<start_edge<<std::endl;
     std::cout<<"<SIBHVS_EID,SIBHVS_LVID>"<<std::endl;

     for (Range::iterator i = _edges.begin(); i != _edges.end(); ++i){
       EntityHandle eid[2];  int lvid[2];
       error = mb->tag_get_data(sibhvs_eid, &*i, 1, eid);
       if (MB_SUCCESS != error) return error;
       error = mb->tag_get_data(sibhvs_lvid, &*i, 1, lvid);
       if (MB_SUCCESS != error) return error;

       if ((eid[0] == 0)&&(eid[1] != 0))
         std::cout<<"<"<<eid[0]<<","<<lvid[0]<<">"<<"      "<<"<"<<(eid[1]-start_edge)<<","<<lvid[1]<<">"<<std::endl;
       else if ((eid[1] == 0)&&(eid[0] != 0))
         std::cout<<"<"<<(eid[0]-start_edge)<<","<<lvid[0]<<">"<<"      "<<"<"<<eid[1]<<","<<lvid[1]<<">"<<std::endl;
     }

     std::cout<<"<V2HV_EID, V2HV_LVID>"<<std::endl;

     for (Range::iterator i = _verts.begin(); i != _verts.end(); ++i){
       EntityHandle eid; int lvid;
       error = mb->tag_get_data(v2hv_eid, &*i, 1, &eid);
       if (MB_SUCCESS != error) return error;
       error = mb->tag_get_data(v2hv_lvid, &*i, 1, &lvid);
       if (MB_SUCCESS != error) return error;

       if (eid != 0)
         std::cout<<"For vertex = "<<*i<<"::Incident halfvertex "<<(eid-start_edge+1)<<"  "<<lvid<<std::endl;
       else
         std::cout<<"For vertex = "<<*i<<"::Incident halfvertex "<<eid<<"  "<<lvid<<std::endl;
     }

     std::cout<<"start_face = "<<start_face<<std::endl;
     std::cout<<"<SIBHES_FID,SIBHES_LEID>"<<std::endl;

     for (Range::iterator i = _faces.begin(); i != _faces.end(); ++i){
       std::vector<EntityHandle> fid(nepf);
       std::vector<int> leid(nepf);
       error = mb->tag_get_data(sibhes_fid, &*i, 1, &fid[0]);
       if (MB_SUCCESS != error) return error;
       error = mb->tag_get_data(sibhes_leid, &*i, 1, &leid[0]);
       if (MB_SUCCESS != error) return error;

       for (int j=0; j<nepf; j++){
         if (fid[j] == 0)
           std::cout<<"<"<<fid[j]<<","<<leid[j]<<">"<<" ";
         else
           std::cout<<"<"<<(fid[j]-start_face+1)<<","<<leid[j]<<">"<<" ";
       }
       std::cout<<std::endl;
     }

     std::cout<<"<V2HE_FID, V2HE_LEID>"<<std::endl;

     for (Range::iterator i = _verts.begin(); i != _verts.end(); ++i){
       EntityHandle fid; int lid;
       error = mb->tag_get_data(v2he_fid, &*i, 1, &fid);
       if (MB_SUCCESS != error) return error;
       error = mb->tag_get_data(v2he_leid, &*i, 1, &lid);
       if (MB_SUCCESS != error) return error;

       if (fid == 0)
         std::cout<<"For vertex = "<<*i<<"::Incident halfedge "<<fid<<"  "<<lid<<std::endl;
       else
         std::cout<<"For vertex = "<<*i<<"::Incident halfedge "<<(fid-start_face+1)<<"  "<<lid<<std::endl;
       }

     std::cout<<"start_cell = "<<start_cell<<std::endl;
     std::cout<<"<SIBHES_CID,SIBHES_LFID>"<<std::endl;

     for (Range::iterator i = _cells.begin(); i != _cells.end(); ++i){
       std::vector<EntityHandle> cid(nfpc);
       std::vector<int> lfid(nfpc);
       error = mb->tag_get_data(sibhfs_cid, &*i, 1, &cid[0]);
       if (MB_SUCCESS != error) return error;
       error = mb->tag_get_data(sibhfs_lfid, &*i, 1, &lfid[0]);
       if (MB_SUCCESS != error) return error;

       for(int j = 0; j < nfpc; j++)
         if (cid[j] == 0)
           std::cout<<"<"<<cid[j]<<","<<lfid[j]<<">"<<" ";
         else
           std::cout<<"<"<<(cid[j]-start_cell+1)<<","<<lfid[j]<<">"<<" ";
       std::cout<<std::endl;
     }
     std::cout<<" Finished printing AHF:sibhfs "<<std::endl;

     std::cout<<"<V2HF_CID, V2HF_LFID>"<<std::endl;

     for (Range::iterator i = _verts.begin(); i != _verts.end(); ++i){
       EntityHandle fid; int lid;
       error = mb->tag_get_data(v2hf_cid, &*i, 1, &fid);
       if (MB_SUCCESS != error) return error;
       error = mb->tag_get_data(v2hf_lfid, &*i, 1, &lid);
       if (MB_SUCCESS != error) return error;

       if (fid==0)
         std::cout<<"For vertex = "<<*i<<"::Incident halfface "<<fid<<"  "<<lid<<std::endl;
       else
         std::cout<<"For vertex = "<<*i<<"::Incident halfface "<<(fid-start_cell+1)<<"  "<<lid<<std::endl;
     }
     return MB_SUCCESS;
   }

   /**********************************************************
   *      User interface for adjacency functions            *
   ********************************************************/

  ErrorCode HalfFacetRep::get_adjacencies(const EntityHandle source_entity,
                                          const unsigned int target_dimension,
                                          std::vector<EntityHandle> &target_entities)
  {

      ErrorCode error;

      unsigned int source_dimension = mb->dimension_from_handle(source_entity);

      if (mInitAHFmaps == false)
      {
       error = initialize();
       if (MB_SUCCESS != error) return error;
      }

      if ((source_dimension == 0) && (target_dimension == 1))
      {
          error = get_up_adjacencies_1d(source_entity, target_entities);
          if (MB_SUCCESS != error) return error;
      }

      else if ((source_dimension == 1) && (target_dimension == 2))
      {
          error = get_up_adjacencies_2d(source_entity, target_entities);
          if (MB_SUCCESS != error) return error;
      }
      else if ((source_dimension == 1) && (target_dimension == 3))
      {
          error = get_up_adjacencies_edg_3d(source_entity, target_entities);
          if (MB_SUCCESS != error) return error;
      }
      else if ((source_dimension == 2) && (target_dimension ==3))
      {
          error = get_up_adjacencies_face_3d(source_entity, target_entities);
          if (MB_SUCCESS != error) return error;
      }
      else if (source_dimension == target_dimension)
      {
          if (target_dimension == 1)
          {
              error = get_neighbor_adjacencies_1d(source_entity, target_entities);
              if (MB_SUCCESS != error) return error;
          }

          else if (target_dimension == 2)
          {
              error = get_neighbor_adjacencies_2d(source_entity, target_entities);
              if (MB_SUCCESS != error) return error;
          }
          else if (target_dimension == 3)
          {
              error = get_neighbor_adjacencies_3d(source_entity, target_entities);
              if (MB_SUCCESS != error) return error;
          }
      }
      return MB_SUCCESS;
  }


   ErrorCode HalfFacetRep::get_up_adjacencies(EntityHandle ent,
                                              int out_dim,
                                              std::vector<EntityHandle> &adjents,
                                              bool local_id,
                                              std::vector<int> * lids)
   {
    ErrorCode error;
    int in_dim = mb->dimension_from_handle(ent);

    if ((in_dim == 0) && (out_dim == 1))
      {
        error = get_up_adjacencies_1d(ent, adjents, local_id, lids);
        if (MB_SUCCESS != error) return error;
      }

    else if ((in_dim == 1) && (out_dim == 2))
      {
        error = get_up_adjacencies_2d(ent, adjents, local_id, lids);
        if (MB_SUCCESS != error) return error;
      }
    else if ((in_dim == 1) && (out_dim == 3))
      {
        error = get_up_adjacencies_edg_3d(ent, adjents, local_id, lids);
        if (MB_SUCCESS != error) return error;
      }
    else if ((in_dim == 2) && (out_dim ==3))
      {
        error = get_up_adjacencies_face_3d(ent, adjents, local_id, lids);
        if (MB_SUCCESS != error) return error;
      }
    return MB_SUCCESS;
   }

   ErrorCode HalfFacetRep::get_neighbor_adjacencies(EntityHandle ent,
                                                    std::vector<EntityHandle> &adjents)
   {
     ErrorCode error;
     int in_dim = mb->dimension_from_handle(ent);

     if (in_dim == 1)
       {
         error = get_neighbor_adjacencies_1d(ent, adjents);
         if (MB_SUCCESS != error) return error;
       }

     else if (in_dim == 2)
       {
         error = get_neighbor_adjacencies_2d(ent, adjents);
         if (MB_SUCCESS != error) return error;
       }
     else if (in_dim == 3)
       {
         error = get_neighbor_adjacencies_3d(ent, adjents);
         if (MB_SUCCESS != error) return error;
       }
     return MB_SUCCESS;
   }

  /******************************************************** 
  * 1D: sibhvs, v2hv, incident and neighborhood queries   *
  *********************************************************/
  ErrorCode HalfFacetRep::determine_sibling_halfverts( Range &edges)
  {
    ErrorCode error;

    EntityHandle start_edge = *edges.begin();
    int count, *sibhvs_lvid_ptr; EntityHandle *sibhvs_eid_ptr;
    error = mb->tag_iterate(sibhvs_eid, edges.begin(), edges.end(), count, (void*&)sibhvs_eid_ptr);
    if (MB_SUCCESS != error) return error;
    assert(count == (int) edges.size());
    error = mb->tag_iterate(sibhvs_lvid, edges.begin(), edges.end(), count, (void*&)sibhvs_lvid_ptr);
    if (MB_SUCCESS != error) return error;
    assert(count == (int) edges.size());

    //Step 1: Create an index list storing the starting position for each vertex
    int nv = _verts.size();
    int *is_index = new int[nv+1];
    for (int i =0; i<nv+1; i++)
      is_index[i] = 0;

    std::vector<EntityHandle> conn(2);
    for (Range::iterator eid = edges.begin(); eid != edges.end(); ++eid)
      {
        conn.clear();
        error = mb->get_connectivity(&*eid, 1, conn);
        if (MB_SUCCESS != error) return error;

        is_index[conn[0]-*_verts.begin()+1] += 1;
        is_index[conn[1]-*_verts.begin()+1] += 1;
      }
    is_index[0] = 0;

    for (int i=0; i<nv; i++)
      is_index[i+1] = is_index[i] + is_index[i+1];

    //Step 2: Define two arrays v2hv_eid, v2hv_lvid storing every half-facet on a vertex
    EntityHandle * v2hv_map_eid = new EntityHandle[2*edges.size()];
    int * v2hv_map_lvid = new int[2*edges.size()];

    for (Range::iterator eid = edges.begin(); eid != edges.end(); ++eid)
      {
        conn.clear();
        error = mb->get_connectivity(&*eid, 1, conn);
        if (MB_SUCCESS != error) return error;

        for (int j = 0; j< 2; j++)
          {
            int v = conn[j] - *_verts.begin();
            v2hv_map_eid[is_index[v]] = *eid;
            v2hv_map_lvid[is_index[v]] = j;
            is_index[v] += 1;
          }
      }

    for (int i=nv-2; i>=0; i--)
      is_index[i+1] = is_index[i];
    is_index[0] = 0;

    //Step 3: Fill up sibling half-verts map
    for (Range::iterator eid = edges.begin(); eid != edges.end(); ++eid)
      {
        conn.clear();
        error = mb->get_connectivity(&*eid, 1, conn);
        if (MB_SUCCESS != error) return error;

        for (int k =0; k<2; k++)
          {
            if (sibhvs_eid_ptr[2*(*eid-start_edge)+k] != 0)
              continue;

            int v = conn[k] - *_verts.begin();
            int last = is_index[v+1] - 1;
            if (last > is_index[v])
              {
                EntityHandle prev_eid = v2hv_map_eid[last];
                int prev_lvid = v2hv_map_lvid[last];
                for (int i=is_index[v]; i<=last; i++)
                  {
                    EntityHandle cur_eid = v2hv_map_eid[i];
                    int cur_lvid = v2hv_map_lvid[i];
                    sibhvs_eid_ptr[2*(prev_eid-start_edge) + prev_lvid] = cur_eid;
                    sibhvs_lvid_ptr[2*(prev_eid-start_edge) + prev_lvid] = cur_lvid;

                    prev_eid = cur_eid;
                    prev_lvid = cur_lvid;
                  }
              }
          }
      }

    std::cout<<"Finished creating sibling half-verts map"<<std::endl;

    delete [] is_index;
    delete [] v2hv_map_eid;
    delete [] v2hv_map_lvid;
    return MB_SUCCESS;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode HalfFacetRep::determine_incident_halfverts( Range &edges){
    ErrorCode error;

    for (Range::iterator e_it = edges.begin(); e_it != edges.end(); ++e_it){      
      std::vector<EntityHandle> conn(2);
      error = mb->get_connectivity(&*e_it, 1, conn);
      if (MB_SUCCESS != error) return error;

      for(int i=0; i<2; ++i){
	EntityHandle v = conn[i], eid = 0;	 
	error = mb->tag_get_data(v2hv_eid, &v, 1, &eid); 
	if (MB_SUCCESS != error) return error;

	if (eid==0){
	  EntityHandle edg = *e_it;	  
	  int lvid = i;

	  error = mb->tag_set_data(v2hv_eid, &v, 1, &edg); 
	  if (MB_SUCCESS != error) return error;
	  error = mb->tag_set_data(v2hv_lvid, &v, 1, &lvid);
	  if (MB_SUCCESS != error) return error;
	}
      }      
    }  
    std::cout<<"Finished creating incident half-verts map"<<std::endl;

    return MB_SUCCESS;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode  HalfFacetRep::get_up_adjacencies_1d( EntityHandle vid,
                                                  std::vector< EntityHandle > &adjents,
                                                  bool local_id,
                                                  std::vector<int> * lvids){
    ErrorCode error;
    EntityHandle start_eid, eid, sibeid[2];
    int start_lid, lid, siblid[2];    
   
    error = mb->tag_get_data(v2hv_eid, &vid, 1, &start_eid);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_data(v2hv_lvid, &vid, 1, &start_lid);    
    if (MB_SUCCESS != error) return error;

    eid = start_eid; lid = start_lid;

    if (eid != 0){
      adjents.push_back(eid);
      if (local_id)
        lvids->push_back(lid);

      while (eid !=0) {	 
	  error = mb->tag_get_data(sibhvs_eid, &eid, 1, sibeid);
	  if (MB_SUCCESS != error) return error;
	  error = mb->tag_get_data(sibhvs_lvid, &eid, 1, siblid);
	  if (MB_SUCCESS != error) return error;

	  eid = sibeid[lid];
	  lid = siblid[lid];
	  if ((!eid)||(eid == start_eid))
	    break;
	  adjents.push_back(eid);
	  if (local_id)
	    lvids->push_back(lid);
      }
    }

    return MB_SUCCESS;
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode  HalfFacetRep::get_neighbor_adjacencies_1d( EntityHandle eid,
                                                        std::vector<EntityHandle> &adjents){

    ErrorCode error;

    EntityHandle sib_eids[2], sibhv_eid;
    int sib_lids[2], sibhv_lid;    

    error = mb->tag_get_data(sibhvs_eid, &eid, 1, sib_eids);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_data(sibhvs_lvid, &eid, 1, sib_lids);
    if (MB_SUCCESS != error) return error;
    
    for (int lid = 0;lid < 2; ++lid){
      sibhv_eid = sib_eids[lid];
      sibhv_lid = sib_lids[lid];

      if (sibhv_eid != 0){
        adjents.push_back(sibhv_eid);
	
	EntityHandle next_eid[2], hv_eid; 
	int next_lid[2], hv_lid;
	
	error = mb->tag_get_data(sibhvs_eid, &sibhv_eid, 1, next_eid);
	if (MB_SUCCESS != error) return error;
	error = mb->tag_get_data(sibhvs_lvid, &sibhv_eid, 1, next_lid);
	if (MB_SUCCESS != error) return error;
	
	hv_eid = next_eid[sibhv_lid];
	hv_lid = next_lid[sibhv_lid];
	
	while (hv_eid != 0){	    
	  if (hv_eid != eid)
	    adjents.push_back(hv_eid);
	  
	  error = mb->tag_get_data(sibhvs_eid, &hv_eid, 1, next_eid);
	  if (MB_SUCCESS != error) return error;
	  error = mb->tag_get_data(sibhvs_lvid, &hv_eid, 1, next_lid);
	  if (MB_SUCCESS != error) return error;

	  if (next_eid[hv_lid] == sibhv_eid)
	    break;
	  hv_eid = next_eid[hv_lid];
	  hv_lid = next_lid[hv_lid];      
	}
      }
    } 

    return MB_SUCCESS;   
  }
  
  /*******************************************************
  * 2D: sibhes, v2he, incident and neighborhood queries  *
  ********************************************************/
  int HalfFacetRep::local_maps_2d( EntityHandle face)
  {
    // nepf: Number of edges per face
    EntityType type = mb->type_from_handle(face);

    int nepf;
    if (type == MBTRI)  nepf = 3;
    else if (type ==MBQUAD)   nepf = 4;
    
    return nepf;
  }
  /////////////////////////////////////////////////////////////////////////////////
  ErrorCode HalfFacetRep::local_maps_2d( int nepf, int *next, int *prev)
  {
    // nepf: Number of edges per face
    // next: Local ids of next edges
    // prev: Local ids of prev edges

    if (nepf == 3)
    {  
      next[0] = 1; next[1] = 2; next[2] = 0;
      prev[0] = 2; prev[1] = 0; prev[2] = 1;
    }
    else if (nepf == 4)
    {
      next[0] = 1; next[1] = 2; next[2] = 3; next[3] = 0;
      prev[0] = 3; prev[1] = 0; prev[2] = 1; prev[3] = 2;
    }

    return MB_SUCCESS;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ErrorCode HalfFacetRep::determine_sibling_halfedges( Range &faces)
  {
    ErrorCode error;
    EntityHandle start_face = *faces.begin();

    // Create pointers to the tags for direct access to memory
    int count, *sibhes_leid_ptr; EntityHandle *sibhes_fid_ptr;
    error = mb->tag_iterate(sibhes_fid, faces.begin(), faces.end(), count, (void*&)sibhes_fid_ptr);
    if (MB_SUCCESS != error) return error;
    assert(count == (int) faces.size());
    error = mb->tag_iterate(sibhes_leid, faces.begin(), faces.end(), count, (void*&)sibhes_leid_ptr);
    if (MB_SUCCESS != error) return error;
    assert(count == (int) faces.size());

    int nepf = local_maps_2d(start_face);
    int * next = new int[nepf];
    int * prev = new int[nepf];
    error = local_maps_2d(nepf, next, prev);
    if (MB_SUCCESS != error) return error;

    //Step 1: Create an index list storing the starting position for each vertex
    int nv = _verts.size();
    int *is_index = new int[nv+1];
    for (int i =0; i<nv+1; i++)
      is_index[i] = 0;

    std::vector<EntityHandle> conn(nepf);
    for (Range::iterator fid = faces.begin(); fid != faces.end(); ++fid)
       {
        conn.clear();
         error = mb->get_connectivity(&*fid, 1, conn);
         if (MB_SUCCESS != error) return error;

         for (int i = 0; i<nepf; i++)
           {
             is_index[conn[i]-*_verts.begin()+1] += 1;
           }
       }
     is_index[0] = 0;

     for (int i=0; i<nv; i++)
       is_index[i+1] = is_index[i] + is_index[i+1];

     //Step 2: Define two arrays v2hv_eid, v2hv_lvid storing every half-facet on a vertex
     EntityHandle * v2nv = new EntityHandle[nepf*faces.size()];
     EntityHandle * v2he_map_fid = new EntityHandle[nepf*faces.size()];
     int * v2he_map_leid = new int[nepf*faces.size()];

     for (Range::iterator fid = faces.begin(); fid != faces.end(); ++fid)
       {
         conn.clear();
         error = mb->get_connectivity(&*fid, 1, conn);
         if (MB_SUCCESS != error) return error;

         for (int j = 0; j< nepf; j++)
           {
             int v = conn[j] - *_verts.begin();
             v2nv[is_index[v]] = conn[next[j]];
             v2he_map_fid[is_index[v]] = *fid;
             v2he_map_leid[is_index[v]] = j;
             is_index[v] += 1;
           }
       }

     for (int i=nv-2; i>=0; i--)
       is_index[i+1] = is_index[i];
     is_index[0] = 0;

     //Step 3: Fill up sibling half-verts map
     for (Range::iterator fid = faces.begin(); fid != faces.end(); ++fid)
       {
         conn.clear();
         error = mb->get_connectivity(&*fid, 1, conn);
         if (MB_SUCCESS != error) return error;

         for (int k =0; k<nepf; k++)
           {
             if (sibhes_fid_ptr[nepf*(*fid-start_face)+k] != 0)
               continue;

             int v = conn[k] - *_verts.begin();
             int vn = conn[next[k]] - *_verts.begin();

             EntityHandle first_fid = *fid;
             int first_leid = k;

             EntityHandle prev_fid = *fid;
             int prev_leid = k;

             for (int index = is_index[vn]; index <= is_index[vn+1]-1; index++)
               {
                 if (v2nv[index] == conn[k])
                   {
                     EntityHandle cur_fid = v2he_map_fid[index];
                     int cur_leid = v2he_map_leid[index];
                     sibhes_fid_ptr[nepf*(prev_fid-start_face)+prev_leid] = cur_fid;
                     sibhes_leid_ptr[nepf*(prev_fid-start_face)+prev_leid] = cur_leid;
                     prev_fid = cur_fid;
                     prev_leid = cur_leid;
                   }
               }

             for (int index = is_index[v]; index <= is_index[v+1]-1; index++)
               {
                 if ((v2nv[index] == conn[next[k]])&&(v2he_map_fid[index] != *fid))
                   {
                     EntityHandle cur_fid = v2he_map_fid[index];
                     int cur_leid = v2he_map_leid[index];
                     sibhes_fid_ptr[nepf*(prev_fid-start_face)+prev_leid] = cur_fid;
                     sibhes_leid_ptr[nepf*(prev_fid-start_face)+prev_leid] = cur_leid;
                     prev_fid = cur_fid;
                     prev_leid = cur_leid;
                   }
               }

             if (prev_fid != first_fid){
               sibhes_fid_ptr[nepf*(prev_fid-start_face) + prev_leid] = first_fid;
               sibhes_leid_ptr[nepf*(prev_fid-start_face) + prev_leid] = first_leid;
             }

           }
       }

     delete [] next;
     delete [] prev;
     delete [] is_index;
     delete [] v2nv;
     delete [] v2he_map_fid;
     delete [] v2he_map_leid;

     std::cout<<"Finished creating sibling half-edges map"<<std::endl;
     return MB_SUCCESS;

  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode HalfFacetRep::determine_incident_halfedges( Range &faces)
  {
    ErrorCode error;    

    int nepf = local_maps_2d(*faces.begin());

    std::vector<EntityHandle> conn(nepf);
    std::vector<EntityHandle> sibfid(nepf);

    for (Range::iterator it = faces.begin(); it != faces.end(); ++it){      
      conn.clear();
      error = mb->get_connectivity(&*it, 1, conn); 
      if (MB_SUCCESS != error) return error;

      sibfid.clear();
      error = mb->tag_get_data(sibhes_fid, &*it, 1, &sibfid[0]);
      if (MB_SUCCESS != error) return error;

      for(int i=0; i<nepf; ++i){
	EntityHandle v = conn[i],  vfid = 0;
	error = mb->tag_get_data(v2he_fid, &v, 1, &vfid);
	if (MB_SUCCESS != error) return error;

	if ((vfid==0)||((vfid!=0) && (sibfid[i]==0))){
	  EntityHandle fid = *it;	  
	  int lid = i;
	  error = mb->tag_set_data(v2he_fid, &v, 1, &fid);
	  if (MB_SUCCESS != error) return error;
	  error = mb->tag_set_data(v2he_leid, &v, 1, &lid);
	  if (MB_SUCCESS != error) return error;
	}
      }      
    }     
    std::cout<<"Finished creating incident half-edges map"<<std::endl;
    
    return MB_SUCCESS;
  }
  ///////////////////////////////////////////////////////////////////
 /* ErrorCode HalfFacetRep::get_up_adjacencies_2d(EntityHandle vid, std::vector<EntityHandle> &adjents)
  {
    ErrorCode error;

    EntityHandle fid; int lid;
    error = mb->tag_get_data(v2he_fid, &vid, 1, &fid);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_data(v2he_leid, &vid, 1, &lid);
    if (MB_SUCCESS != error) return error;

    if (fid != 0)
      {
        adjents.push_back(fid);

        EntityHandle queue_fid[MAXSIZE], trackfaces[MAXSIZE];
        int queue_lid[MAXSIZE];
        for (int i = 0; i< MAXSIZE; i++)
          {
            queue_fid[i] = 0;
            queue_lid[i] = 0;
            trackfaces[i] = 0;
          }
        int qsize = 0, tcount = -1;
        int num_qvals = 0;
        error = gather_halfedges(vid, fid, lid, queue_fid, queue_lid, &qsize, trackfaces, &tcount);
        if (MB_SUCCESS != error) return error;

        while (num_qvals < *qsize)
          {
            EntityHandle curfid = queue_fid[num_qvals];
            int curlid = queue_lid[num_qvals];
            num_qvals += 1;

            EntityHandle he2_fid; int he2_lid;
            error = another_halfedge(vid, curfid, curlid, &he2_fid, &he2_lid);
            if (MB_SUCCESS != error) return error;

            bool found_ent = find_match_in_array(he2_fid, trackfaces, tcount[0]);

            if (found_ent)
              continue;
            tcount[0] += 1;
            trackfaces[tcount[0]] = he2_fid;

            error = get_up_adjacencies_2d(he2_fid, he2_lid, queue_fid, queue_lid, qsize, trackfaces, tcount);
            if (MB_SUCCESS != error) return error;

            adjents.push_back(he2_fid);

          }
      }
  }
*/
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode  HalfFacetRep::get_up_adjacencies_2d( EntityHandle eid,
                                                  std::vector<EntityHandle> &adjents,
                                                  bool local_id,
                                                  std::vector<int>  * leids)
{

  // Given an explicit edge eid, find the incident faces.

    ErrorCode error;
    EntityHandle he_fid=0; int he_lid=0;

    // Step 1: Given an explicit edge, find a corresponding half-edge from the surface mesh.
    bool found = find_matching_halfedge(eid, &he_fid, &he_lid);

    // Step 2: If there is a corresponding half-edge, collect all sibling half-edges and store the incident faces.
    if (found)
      { 
        error = get_up_adjacencies_2d(he_fid, he_lid, true, adjents, local_id, leids);
        if (MB_SUCCESS != error) return error;
      }

    return MB_SUCCESS;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode HalfFacetRep::get_up_adjacencies_2d( EntityHandle fid,
                                                 int leid,
                                                 bool add_inent,
                                                 std::vector<EntityHandle> &fids,
                                                 bool local_id,
                                                 std::vector<int> * leids)
  {
    // Given an implicit half-edge <fid, leid>, find the incident half-edges.
    ErrorCode error;
    int nepf = local_maps_2d(fid);

    if (!fid)
      return MB_FAILURE;

    if (add_inent)
      {
        fids.push_back(fid);
        if (local_id)
          leids->push_back(leid);
      }

    std::vector<EntityHandle> sib_fids(nepf);
    std::vector<int> sib_lids(nepf);
    error = mb->tag_get_data(sibhes_fid, &fid, 1, &sib_fids[0]);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_data(sibhes_leid, &fid, 1, &sib_lids[0]);
    if (MB_SUCCESS != error) return error;

    EntityHandle curfid = sib_fids[leid];
    int curlid = sib_lids[leid];
    
    while ((curfid != fid)&&(curfid != 0)){//Should no go into the loop when no sibling exists
        fids.push_back(curfid);
        if (local_id)
          leids->push_back(curlid);

      sib_fids.clear();
      sib_lids.clear();
      error = mb->tag_get_data(sibhes_fid, &curfid, 1, &sib_fids[0]);
      if (MB_SUCCESS != error) return error;
      error = mb->tag_get_data(sibhes_leid, &curfid, 1, &sib_lids[0]);
      if (MB_SUCCESS != error) return error;
      
      curfid = sib_fids[curlid];
      curlid = sib_lids[curlid];
    }

    return MB_SUCCESS;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ErrorCode HalfFacetRep::get_up_adjacencies_2d(EntityHandle fid,
                                                 int lid,
                                                 EntityHandle *queue_fid,
                                                 int *queue_lid,
                                                 int *qsize,
                                                 EntityHandle *trackfaces,
                                                 int *tcount)
  {

    ErrorCode error; 
    int nepf = local_maps_2d(fid);
    
    std::vector<EntityHandle> sib_fids(nepf);
    std::vector<int> sib_lids(nepf);
    error = mb->tag_get_data(sibhes_fid, &fid, 1, &sib_fids[0]);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_data(sibhes_leid, &fid, 1, &sib_lids[0]);
    if (MB_SUCCESS != error) return error;
    EntityHandle curfid = sib_fids[lid];
    int curlid = sib_lids[lid];

    if (curfid == 0){
        int index = 0;
        bool found_ent = find_match_in_array(fid, queue_fid, qsize[0]-1, true, &index);
        if ((!found_ent)||((found_ent) && (queue_lid[index] != lid)))

          {
            queue_fid[qsize[0]] = fid;
            queue_lid[qsize[0]] = lid;
            qsize[0] += 1;
          }
      }

    while ((curfid != fid)&&(curfid != 0)) {

        bool found = find_match_in_array(curfid, trackfaces, tcount[0]);

        if (!found){
            queue_fid[qsize[0]] = curfid;
            queue_lid[qsize[0]] = curlid;
            qsize[0] += 1;

          }

        sib_fids.clear();
        sib_lids.clear();
        error = mb->tag_get_data(sibhes_fid, &curfid, 1, &sib_fids[0]);
        if (MB_SUCCESS != error) return error;
        error = mb->tag_get_data(sibhes_leid, &curfid, 1, &sib_lids[0]);
        if (MB_SUCCESS != error) return error;

        curfid = sib_fids[curlid];
        curlid = sib_lids[curlid];
      }
    
    return MB_SUCCESS;
   }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool HalfFacetRep::find_matching_halfedge( EntityHandle eid,
                                             EntityHandle *hefid,
                                             int *helid){
    ErrorCode error;
    std::vector<EntityHandle> conn(2);
    error = mb->get_connectivity(&eid, 1, conn);
    if (MB_SUCCESS != error) return error;

    EntityHandle fid; int lid;
    error = mb->tag_get_data(v2he_fid, &conn[0], 1, &fid);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_data(v2he_leid, &conn[0], 1, &lid);
    if (MB_SUCCESS != error) return error;

    bool found = false;

    if (fid!=0){
        EntityHandle queue_fid[MAXSIZE], trackfaces[MAXSIZE];
        int queue_lid[MAXSIZE];
        for (int i = 0; i< MAXSIZE; i++)
          {
            queue_fid[i] = 0;
            queue_lid[i] = 0;
            trackfaces[i] = 0;
          }
        int qsize = 0, tcount = -1;

        EntityHandle vid = conn[0];

        error = gather_halfedges(vid, fid, lid, queue_fid, queue_lid, &qsize, trackfaces, &tcount);
        if (MB_SUCCESS != error) return error;

        found =  collect_and_compare(conn, queue_fid, queue_lid, &qsize, trackfaces, &tcount, hefid, helid);
      }

    return found;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ErrorCode HalfFacetRep::gather_halfedges( EntityHandle vid,
                                             EntityHandle he_fid,
                                             int he_lid,
                                             EntityHandle *queue_fid,
                                             int *queue_lid,
                                             int *qsize,
                                             EntityHandle *trackfaces,
                                             int *tcount)
  {  
    ErrorCode error;
    EntityHandle he2_fid = 0; int he2_lid = 0;

    error = another_halfedge(vid, he_fid, he_lid, &he2_fid, &he2_lid);
    if (MB_SUCCESS != error) return error;

    queue_fid[*qsize] = he_fid; queue_lid[*qsize] = he_lid;
    *qsize += 1;
    queue_fid[*qsize] = he2_fid; queue_lid[*qsize] = he2_lid;
    *qsize += 1;
    tcount[0] += 1;
    trackfaces[tcount[0]] = he_fid;
    
    error = get_up_adjacencies_2d(he_fid, he_lid, queue_fid, queue_lid, qsize, trackfaces, tcount);
    if (MB_SUCCESS != error) return error;
    error = get_up_adjacencies_2d(he2_fid, he2_lid, queue_fid, queue_lid, qsize, trackfaces, tcount);
    if (MB_SUCCESS != error) return error;

    return MB_SUCCESS;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode HalfFacetRep::another_halfedge( EntityHandle vid,
                                            EntityHandle he_fid,
                                            int he_lid,
                                            EntityHandle *he2_fid,
                                            int *he2_lid)
  {    
    ErrorCode error;
    int nepf;
    nepf = local_maps_2d(he_fid);
    int * next = new int[nepf];
    int * prev = new int[nepf];
    error = local_maps_2d(nepf, next, prev);
    if (MB_SUCCESS != error) return error;
    
    std::vector<EntityHandle> conn(nepf);
    error = mb->get_connectivity(&he_fid, 1, conn);
    if (MB_SUCCESS != error) return error;

    *he2_fid = he_fid;
    if (conn[he_lid] == vid)
      *he2_lid = prev[he_lid];
    else
      *he2_lid = next[he_lid];

    delete [] next;
    delete [] prev;
    return MB_SUCCESS;
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool HalfFacetRep::collect_and_compare(std::vector<EntityHandle> &edg_vert,
                                         EntityHandle *queue_fid,
                                         int *queue_lid,
                                         int *qsize,
                                         EntityHandle *trackfaces,
                                         int *tcount,
                                         EntityHandle *he_fid,
                                         int *he_lid)
  {
    ErrorCode error;
    int nepf = local_maps_2d(*_faces.begin());
    int *next = new int[nepf];
    int *prev = new int[nepf];
    error = local_maps_2d(nepf, next, prev);
    if (MB_SUCCESS != error) return error;

    bool found = false;
    int num_qvals = 0, counter = 0;

    while (num_qvals < *qsize && counter <MAXSIZE )
      {
        EntityHandle curfid = queue_fid[num_qvals];
        int curlid = queue_lid[num_qvals];
        num_qvals += 1;       

	std::vector<EntityHandle> conn(nepf);
	error = mb->get_connectivity(&curfid, 1, conn);
	if (MB_SUCCESS != error) return error;

	int id = next[curlid];
	if (((conn[curlid]==edg_vert[0])&&(conn[id]==edg_vert[1]))||((conn[curlid]==edg_vert[1])&&(conn[id]==edg_vert[0]))){
	    *he_fid = curfid;
	    *he_lid = curlid;
	    found = true;
	    break;
	}

	bool found_ent = find_match_in_array(curfid, trackfaces, tcount[0]);

	if (found_ent)
	  continue;
	tcount[0] += 1;
	trackfaces[tcount[0]] = curfid;

	EntityHandle he2_fid; int he2_lid;
	error = another_halfedge(edg_vert[0], curfid, curlid, &he2_fid, &he2_lid);

	if (MB_SUCCESS != error) return error;
	error = get_up_adjacencies_2d(he2_fid, he2_lid, queue_fid, queue_lid, qsize, trackfaces, tcount);
	if (MB_SUCCESS != error) return error;      

	counter += 1;
      }
      
    delete [] next; 
    delete [] prev;
    return found;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode  HalfFacetRep::get_neighbor_adjacencies_2d( EntityHandle fid,
                                                        std::vector<EntityHandle> &adjents)
  {
    ErrorCode error; 
    if (fid != 0){
      
      int nepf = local_maps_2d(fid);
      for (int lid = 0; lid < nepf; ++lid){
        error = get_up_adjacencies_2d(fid, lid, false, adjents);
	if (MB_SUCCESS != error) return error;
      }
    }
    
    return MB_SUCCESS;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  int HalfFacetRep::find_total_edges_2d(Range &faces)
  {
    ErrorCode error;
    int firstF = *faces.begin();
    int nepf = local_maps_2d(firstF);
    int nfaces = faces.size();

    int total_edges = nepf*nfaces;
    
    bool *trackF = new bool(total_edges);
    std::vector<EntityHandle> adj_fids;
    std::vector<int> adj_lids;

    for (Range::iterator f = faces.begin(); f != faces.end(); f++){
      for (int l = 0; l < nepf; l++){

	adj_fids.clear();
	adj_lids.clear();

	int id = nepf*(*f-firstF)+l;
	if (!trackF[id])
	  {
	    error = get_up_adjacencies_2d(*f,l, false, adj_fids, true, &adj_lids);
	    if (MB_SUCCESS != error) return error;

	    total_edges -= adj_fids.size();

	    for (int i = 0; i < (int)adj_fids.size(); i++)
	      trackF[nepf*(adj_fids[i]-firstF)+adj_lids[i]] = true;
	  };
      };
   };
    
    delete [] trackF;        
    return total_edges;
  }

  /*******************************************************
  * 3D: sibhfs, v2hf, incident and neighborhood queries  *
  ********************************************************/

  int HalfFacetRep::get_index_from_type(EntityHandle cid)
  {
    int index;
    EntityType type = mb->type_from_handle(cid);
    if (type == MBTET)
      index = 0;
    else if (type == MBPYRAMID)
      index = 1;
    else if (type == MBPRISM)
      index = 2;
    else if (type == MBHEX)
      index = 3;

    return index;   
  }

  /////////////////////////////////////////////
   HalfFacetRep::LocalMaps3D HalfFacetRep::lConnMap3D[4] = 
    {
      // Tet
      {4, 6, 4, {3,3,3,3}, {{0,1,3},{1,2,3},{2,0,3},{0,2,1}},   {3,3,3,3},   {{0,2,3},{0,1,3},{1,2,3},{0,1,2}},   {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}},   {{3,0},{3,1},{3,2},{0,2},{0,1},{1,2}},   {{0,4,3},{1,5,4},{2,3,5},{2,1,0}},     {{-1,0,2,3},{0,-1,1,4},{2,1,-1,5},{3,4,5,-1}}},

      // Pyramid: Note: In MOAB pyramid follows the CGNS convention. Look up src/MBCNArrays.hpp
      {5, 8, 5, {4,3,3,3,3}, {{0,3,2,1},{0,1,4},{1,2,4},{2,3,4},{3,0,4}},  {3,3,3,3,4},   {{0,1,4},{0,1,2},{0,2,3},{0,3,4},{1,2,3,4}},   {{0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4}},   {{0,1},{0,2},{0,3},{0,4},{1,4},{1,2},{2,3},{3,4}},    {{3,2,1,0},{0,5,4},{1,6,5},{2,7,6},{3,4,7}},    {{-1,0,-1,3,4},{0,-1,1,-1,5},{-1,1,-1,2,6},{3,-1,2,-1,7},{4,5,6,7,-1}}},

      // Prism
      {6, 9, 5, {4,4,4,3,3}, {{0,1,4,3},{1,2,5,4},{0,3,5,2},{0,2,1},{3,4,5}},  {3,3,3,3,3,3}, {{0,2,3},{0,1,3},{1,2,3},{0,2,4},{0,1,4},{1,4,2}},    {{0,1},{1,2},{2,0},{0,3},{1,4},{2,5},{3,4},{4,5},{5,3}},    {{0,3},{1,3},{2,3},{0,2},{0,1},{1,2},{0,4},{1,4},{2,4}},     {{0,4,6,3},{1,5,7,4},{2,3,8,5},{2,1,0},{6,7,8}},    {{-1,0,2,3,-1,-1},{0,-1,1,-1,4,-1},{2,1,-1,-1,-1,5},{3,-1,-1,-1,6,8},{-1,4,-1,6,-1,7},{-1,-1,5,8,7,-1}}},

      // Hex
      {8, 12, 6, {4,4,4,4,4,4}, {{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{0,3,2,1},{4,5,6,7}},   {3,3,3,3,3,3,3,3},   {{0,3,4},{0,1,4},{1,2,4},{2,3,4},{0,3,5},{0,1,5},{1,2,5},{2,3,5}},    {{0,1},{1,2},{2,3},{3,0},{0,4},{1,5},{2,6},{3,7},{4,5},{5,6},{6,7},{7,4}},     {{0,4},{1,4},{2,4},{3,4},{0,3},{0,1},{1,2},{2,3},{0,5},{1,5},{2,5},{3,5}},     {{0,5,8,4},{1,6,9,5},{2,7,10,6},{3,4,11,7},{3,2,1,0},{8,9,10,11}},     {{-1,0,-1,3,4,-1,-1,-1},{0,-1,1,-1,-1,5,-1,-1},{-1,1,-1,2,-1,-1,6,-1},{3,-1,2,-1,-1,-1,-1,7},{4,-1,-1,-1,-1,8,-1,11},{-1,5,-1,-1,8,-1,9,-1},{-1,-1,6,-1,-1,9,-1,10},{-1,-1,-1,7,11,-1,10,-1}}}
      
    };

  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ErrorCode HalfFacetRep::determine_sibling_halffaces( Range &cells)
  {
    ErrorCode error;
    EntityHandle start_cell = *cells.begin();

    int count, *sibhfs_lfid_ptr;
    EntityHandle *sibhfs_cid_ptr;
    error = mb->tag_iterate(sibhfs_cid, cells.begin(), cells.end(), count, (void*&)sibhfs_cid_ptr);
    if (MB_SUCCESS != error) return error;
    assert(count == (int) cells.size());
    error = mb->tag_iterate(sibhfs_lfid, cells.begin(), cells.end(), count, (void*&)sibhfs_lfid_ptr);
    if (MB_SUCCESS != error) return error;
    assert(count == (int) cells.size());

    int index = get_index_from_type(start_cell);
    int nvpc = lConnMap3D[index].num_verts_in_cell;
    int nfpc = lConnMap3D[index].num_faces_in_cell;

    //Step 1: Create an index list storing the starting position for each vertex
    int nv = _verts.size();
    int *is_index = new int[nv+1];
    for (int i =0; i<nv+1; i++)
      is_index[i] = 0;

    std::vector<EntityHandle> conn(nvpc);
    for (Range::iterator cid = cells.begin(); cid != cells.end(); ++cid)
       {
        conn.clear();
        error = mb->get_connectivity(&*cid, 1, conn);
        if (MB_SUCCESS != error) return error;

        for (int i = 0; i<nfpc; ++i)
          {
            int nvF = lConnMap3D[index].hf2v_num[i];
            EntityHandle v = 0;
            for (int k = 0; k<nvF; k++)
              {
                int id = lConnMap3D[index].hf2v[i][k];
                if (v <= conn[id])
                  v = conn[id];
              }
            is_index[v-*_verts.begin()+1] += 1;
          }
      }
     is_index[0] = 0;

     for (int i=0; i<nv; i++)
       is_index[i+1] = is_index[i] + is_index[i+1];

     //Step 2: Define four arrays v2hv_eid, v2hv_lvid storing every half-facet on a vertex
     EntityHandle * v2oe_v1 = new EntityHandle[is_index[nv]];
     EntityHandle * v2oe_v2 = new EntityHandle[is_index[nv]];
     EntityHandle * v2hf_map_cid = new EntityHandle[is_index[nv]];
     int * v2hf_map_lfid = new int[is_index[nv]];

     for (Range::iterator cid = cells.begin(); cid != cells.end(); ++cid)
       {
         conn.clear();
         error = mb->get_connectivity(&*cid, 1, conn);
         if (MB_SUCCESS != error) return error;

         for (int i = 0; i< nfpc; i++)
           {
             int nvF = lConnMap3D[index].hf2v_num[i];
             EntityHandle *vs = new EntityHandle[nvF];
             EntityHandle vmax = 0;
             int lv = -1;
             for (int k = 0; k<nvF; k++)
               {
                 int id = lConnMap3D[index].hf2v[i][k];
                 vs[k] = conn[id];
                 if (vmax <= conn[id])
                   {
                     vmax = conn[id];
                     lv = k;
                   }
               }

             int * next = new int[nvF];
             int * prev = new int[nvF];
             error = local_maps_2d(nvF, next, prev);
             if (MB_SUCCESS != error) return error;

             int v = vmax - *_verts.begin();
             v2oe_v1[is_index[v]] = vs[next[lv]];
             v2oe_v2[is_index[v]] = vs[prev[lv]];
             v2hf_map_cid[is_index[v]] = *cid;
             v2hf_map_lfid[is_index[v]] = i;
             is_index[v] += 1;

             delete [] next;
             delete [] prev;
             delete [] vs;
           }
       }

     for (int i=nv-2; i>=0; i--)
       is_index[i+1] = is_index[i];
     is_index[0] = 0;

     //Step 3: Fill up sibling half-verts map
     for (Range::iterator cid = cells.begin(); cid != cells.end(); ++cid)
       {
         conn.clear();
         error = mb->get_connectivity(&*cid, 1, conn);
         if (MB_SUCCESS != error) return error;

         for (int i =0; i<nfpc; i++)
           {
             if (sibhfs_cid_ptr[nfpc*(*cid-start_cell)+i] != 0)
               continue;


             int nvF = lConnMap3D[index].hf2v_num[i];
             EntityHandle *vs = new EntityHandle[nvF];
             EntityHandle vmax = 0;
             int lv = -1;
             for (int k = 0; k<nvF; k++)
               {
                 int id = lConnMap3D[index].hf2v[i][k];
                 vs[k] = conn[id];
                 if (vmax <= conn[id])
                   {
                     vmax = conn[id];
                     lv = k;
                   }
               }

             int * next = new int[nvF];
             int * prev = new int[nvF];
             error = local_maps_2d(nvF, next, prev);
             if (MB_SUCCESS != error) return error;

             int v = vmax - *_verts.begin();
             EntityHandle v1 = vs[prev[lv]];
             EntityHandle v2 = vs[next[lv]];

             for (int ind = is_index[v]; ind <= is_index[v+1]-1; ind++)
               {
                 if ((v2oe_v1[ind] == v1)&&(v2oe_v2[ind] == v2))
                   {

                     EntityHandle cur_cid = v2hf_map_cid[ind];
                     int cur_lfid = v2hf_map_lfid[ind];

                     sibhfs_cid_ptr[nfpc*(*cid-start_cell)+i] = cur_cid;
                     sibhfs_lfid_ptr[nfpc*(*cid-start_cell)+i] = cur_lfid;

                     sibhfs_cid_ptr[nfpc*(cur_cid-start_cell)+cur_lfid] = *cid;
                     sibhfs_lfid_ptr[nfpc*(cur_cid-start_cell)+cur_lfid] = i;

                   }
               }

             delete [] next;
             delete [] prev;
             delete [] vs;
           }
       }


     delete [] is_index;
     delete [] v2oe_v1;
     delete [] v2oe_v2;
     delete [] v2hf_map_cid;
     delete [] v2hf_map_lfid;

     std::cout<<"Finished creating sibling half-faces map"<<std::endl;

     return MB_SUCCESS;

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode HalfFacetRep::determine_incident_halffaces( Range &cells)
  {
    ErrorCode error;

    int index = get_index_from_type(*cells.begin());
    int nvpc = lConnMap3D[index].num_verts_in_cell;
    int nfpc = lConnMap3D[index].num_faces_in_cell;
  
    std::vector<EntityHandle> conn(nvpc);
    std::vector<EntityHandle> sib_cids(nfpc);

    for (Range::iterator t_it = cells.begin(); t_it != cells.end(); ++t_it){      
      EntityHandle tet = *t_it;
      conn.clear();
      error = mb->get_connectivity(&tet, 1, conn);     
      if (MB_SUCCESS != error) return error;

      sib_cids.clear(); 
      error = mb->tag_get_data(sibhfs_cid, &tet, 1, &sib_cids[0]);
      if (MB_SUCCESS != error) return error;
      
      for(int i=0; i<nvpc; ++i){
	EntityHandle v = conn[i], vcid;	
	error = mb->tag_get_data(v2hf_cid, &v, 1, &vcid);
	if (MB_SUCCESS != error) return error;

	int nhf_pv = lConnMap3D[index].v2hf_num[i];

	for (int j=0; j < nhf_pv; ++j){
	  int ind = lConnMap3D[index].v2hf[i][j];
	  if (vcid==0){
	    error = mb->tag_set_data(v2hf_cid, &v, 1, &tet);
	    if (MB_SUCCESS != error) return error;
	    error = mb->tag_set_data(v2hf_lfid, &v, 1, &ind);
	    if (MB_SUCCESS != error) return error;
	    break;
	  }
	  else if ((vcid!=0) && (sib_cids[ind]==0)){
	    error = mb->tag_set_data(v2hf_cid, &v, 1, &tet);
	    if (MB_SUCCESS != error) return error;
	    error = mb->tag_set_data(v2hf_lfid, &v, 1, &ind);
	    if (MB_SUCCESS != error) return error;
	  }
	}
      }
    }
    
    std::cout<<"Finished creating incident half-faces"<<std::endl;
    return MB_SUCCESS;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode  HalfFacetRep::determine_border_vertices( Range &cells, Tag isborder)
  {
    ErrorCode error;
    EntityHandle start_cell = *cells.begin();
    
    int index = get_index_from_type(start_cell);
    int nvpc = lConnMap3D[index].num_verts_in_cell;
    int nfpc = lConnMap3D[index].num_faces_in_cell;

    int val= 1;
    std::vector<EntityHandle> conn(nvpc);
    std::vector<EntityHandle> sib_cids(nfpc);

    for(Range::iterator t= cells.begin(); t !=cells.end(); ++t){
      conn.clear();     
      error = mb->get_connectivity(&*t, 1, conn);
      if (MB_SUCCESS != error) return error;

      sib_cids.clear();
      error = mb->tag_get_data(sibhfs_cid, &*t, 1, &sib_cids[0]);   
      if (MB_SUCCESS != error) return error;
      
      for (int i = 0; i < nfpc; ++i){	
	if (sib_cids[i]==0){
	  int nvF = lConnMap3D[index].hf2v_num[i];
	  
	  for (int j = 0; j < nvF; ++j){
	    int ind = lConnMap3D[index].hf2v[i][j];
	    error = mb->tag_set_data(isborder, &conn[ind], 1, &val);	    
	    if (MB_SUCCESS != error) return error;
	  }
	}	
      }
    }    

    return MB_SUCCESS;
  }

 //////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode HalfFacetRep::get_up_adjacencies_edg_3d( EntityHandle eid,
                                                     std::vector<EntityHandle> &adjents,
                                                     bool local_id,
                                                     std::vector<int> * leids)
  {
    ErrorCode error; 
    
    EntityHandle cid=0;
    int leid=0;
    //Find one incident cell
    bool found = find_matching_implicit_edge_in_cell(eid, &cid, &leid);

    //Find all incident cells     
    if (found)
      {
        error =get_up_adjacencies_edg_3d(cid, leid, adjents, local_id, leids);
        if (MB_SUCCESS != error) return error;
      }

    return MB_SUCCESS;
  }

  //////////////////////////////////////////////////////////////
  ErrorCode HalfFacetRep::get_up_adjacencies_edg_3d( EntityHandle cid,
                                                     int leid,
                                                     std::vector<EntityHandle> &adjents,
                                                     bool local_id,
                                                     std::vector<int> * leids)
  {
    ErrorCode error;

    int index = get_index_from_type(cid);
    int nvpc = lConnMap3D[index].num_verts_in_cell;
    int nfpc = lConnMap3D[index].num_faces_in_cell;

    adjents.push_back(cid);
    if (local_id)
      leids->push_back(leid);

    std::vector<EntityHandle> conn(nvpc);
    error =mb->get_connectivity(&cid, 1, conn);
    if (MB_SUCCESS != error) return error;

    // Get the end vertices of the edge <cid,leid>
    int id = lConnMap3D[index].e2v[leid][0];
    EntityHandle vert0 = conn[id];
    id = lConnMap3D[index].e2v[leid][1];
    EntityHandle vert1 = conn[id];

    //Loop over the two incident half-faces on this edge
    for(int i=0; i<2; i++)
      {
        EntityHandle cur_cell = cid;
        int cur_leid = leid;

        int lface = i;

	while(true){
	  int lfid = lConnMap3D[index].e2hf[cur_leid][lface];
	  
	  std::vector<EntityHandle> sibcids(nfpc);
	  std::vector<int> siblfids(nfpc);
	  error = mb->tag_get_data(sibhfs_cid, &cur_cell,1, &sibcids[0]);
	  if (MB_SUCCESS != error) return error;
	  error = mb->tag_get_data(sibhfs_lfid, &cur_cell, 1, &siblfids[0]);
	  if (MB_SUCCESS != error) return error;
	  
	  cur_cell = sibcids[lfid];
	  lfid = siblfids[lfid];
	  
	  //Check if loop reached starting cell or a boundary
	  if ((cur_cell == cid) || ( cur_cell==0))
	    break;
	  
	  std::vector<EntityHandle> sib_conn(nvpc);
	  error =mb->get_connectivity(&cur_cell, 1, sib_conn);
	  if (MB_SUCCESS != error) return error;

	  //Find the local edge id wrt to sibhf
	  int nv_curF = lConnMap3D[index].hf2v_num[lfid];
	  int lv0 = -1, lv1 = -1, idx = -1 ;
	  for (int j = 0; j<nv_curF; j++)
	    {
	      idx = lConnMap3D[index].hf2v[lfid][j];
	      if (vert0 == sib_conn[idx])
		lv0 = idx;
	      if (vert1 == sib_conn[idx])
		lv1 = idx;
	    }

	  assert((lv0 >= 0) && (lv1 >= 0));
	  cur_leid = lConnMap3D[index].lookup_leids[lv0][lv1];

	  int chk_lfid = lConnMap3D[index].e2hf[cur_leid][0];

	  if (lfid == chk_lfid)
	    lface = 1;
	  else
	    lface = 0;

	  //Memory allocation
	  if (adjents.capacity()<=adjents.size()){
	    if (100<adjents.size()){
		std::cout<<"Too many incident cells"<<std::endl;
		return MB_FAILURE;
	    }
	    adjents.reserve(20+adjents.size());
	    if (local_id)
	      leids->reserve(20+leids->size());
	  }

	  //insert new incident tetrahedra and local edge ids
	  adjents.push_back(cur_cell);
	  if (local_id)
	    leids->push_back(cur_leid);
	}

	//Loop back
	if (cur_cell != 0)
	  break;
      }
    return MB_SUCCESS;
  }
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode  HalfFacetRep::get_up_adjacencies_face_3d( EntityHandle fid,
                                                       std::vector<EntityHandle> &adjents,
                                                       bool local_id,
                                                       std::vector<int> * lfids)
  {
    ErrorCode error;

    EntityHandle cid = 0;
    int lid = 0;
    bool found = find_matching_halfface(fid, &cid, &lid);

    if (found){
      error = get_up_adjacencies_face_3d(cid, lid, adjents, local_id, lfids);
      if (MB_SUCCESS != error) return error;
      }

    return MB_SUCCESS;
  }
  ///////////////////////////////////////////
  ErrorCode HalfFacetRep::get_up_adjacencies_face_3d( EntityHandle cid,
                                                      int lfid,
                                                      std::vector<EntityHandle> &adjents,
                                                      bool local_id,
                                                      std::vector<int> * lfids)
  {
    ErrorCode error;

    EntityHandle start_cell = *_cells.begin();
    int index = get_index_from_type(start_cell);
    int nfpc = lConnMap3D[index].num_faces_in_cell;

    std::vector<EntityHandle> sibcids(nfpc);
    std::vector<int> siblids(nfpc);
    error = mb->tag_get_data(sibhfs_cid, &cid, 1, &sibcids[0]);
    if (MB_SUCCESS != error) return error;
    error = mb->tag_get_data(sibhfs_lfid, &cid, 1, &siblids[0]);
    if (MB_SUCCESS != error) return error;

    EntityHandle sibcid = sibcids[lfid];
    int siblid = siblids[lfid];
    if (sibcid > 0)
      {
        adjents.push_back(cid);
        adjents.push_back(sibcid);
        if (local_id)
          {
            lfids->push_back(lfid);
            lfids->push_back(siblid);
          }
      }
    else if (sibcid == 0)
      {
        adjents.push_back(cid);
        if (local_id)
          lfids->push_back(lfid);
      }

    return MB_SUCCESS;
  }
  /////////////////////////////////////////////////////////////////
 bool HalfFacetRep::find_matching_implicit_edge_in_cell( EntityHandle eid,
                                                         EntityHandle *cid,
                                                         int *leid)
  {
    ErrorCode error;

    // Find the edge vertices
    std::vector<EntityHandle> econn(2);
    error = mb->get_connectivity(&eid, 1, econn);
    if (MB_SUCCESS != error) return error;

    EntityHandle v_start = econn[0], v_end = econn[1];
    
    // Find an incident cell to v_start
    EntityHandle cell2origin, cell2end;
    error =mb->tag_get_data(v2hf_cid, &v_start, 1, &cell2origin);
    if (MB_SUCCESS != error) return error;
    error =mb->tag_get_data(v2hf_cid, &v_end, 1, &cell2end);
    if (MB_SUCCESS != error) return error;
    
    bool found = false;
    if (cell2origin == 0|| cell2end == 0){
      return found;
    }
    
    EntityHandle start_cell = *_cells.begin();
    int index = get_index_from_type(start_cell);
    int nvpc = lConnMap3D[index].num_verts_in_cell;
    int nfpc = lConnMap3D[index].num_faces_in_cell;

    // Create search queue to add cells
    EntityHandle cellq[MAXSIZE];
    for (int i = 0; i< MAXSIZE; i++)
      {
        cellq[i] = 0;
      }
    cellq[0] = cell2origin;
    cellq[1] = cell2end;

    int qsize = 2, num_qvals = 0;

    while (num_qvals < qsize){
        EntityHandle cell_id = cellq[num_qvals];
        num_qvals += 1;

        std::vector<EntityHandle> conn(nvpc);
        error =mb->get_connectivity(&cell_id, 1, conn);
        if (MB_SUCCESS != error) return error;

        int lv0 = -1, lv1 = -1, lv;

        //locate v_origin in poped out tet, check if v_end is in
        for (int i = 0; i<nvpc; i++){
            if (v_start == conn[i]){
                lv0 = i;
                lv = lv0;
              }
            else if (v_end == conn[i]){
                lv1 = i;
                lv = lv1;
              }
          }

        if ((lv0 >= 0) && (lv1 >= 0))
          {
            found = true;
            *cid = cell_id;
            *leid = lConnMap3D[index].lookup_leids[lv0][lv1];
            break;
          }

        //push back new found unchecked incident tets of v_start
        std::vector<EntityHandle> ngb_cids(nfpc);
        error =mb->tag_get_data(sibhfs_cid,&cell_id,1,&ngb_cids[0]);
        if (MB_SUCCESS != error) return error;

        int nhf_thisv = lConnMap3D[index].v2hf_num[lv];

        for (int i = 0; i < nhf_thisv; i++){
            int ind = lConnMap3D[index].v2hf[lv][i];
            EntityHandle ngb = ngb_cids[ind];

            if (ngb){
                bool found_ent = find_match_in_array(ngb, &cellq[0], qsize-1);

                if (!found_ent)
                  {
                    cellq[qsize] = ngb;
                    qsize += 1;
                  }
              }
          }
      }
    
    return found;
 }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool HalfFacetRep::find_matching_halfface(EntityHandle fid, EntityHandle *cid, int *lid)
  {
    ErrorCode error;
    EntityHandle start_cell = *_cells.begin();
    int index = get_index_from_type(start_cell);
    int nvpc = lConnMap3D[index].num_verts_in_cell;
    int nfpc = lConnMap3D[index].num_faces_in_cell;
    int nvF = local_maps_2d(fid);

    std::vector<EntityHandle> fid_verts(nvF);
    error = mb->get_connectivity(&fid, 1, fid_verts);
    if (MB_SUCCESS != error) return error;

    EntityHandle cur_cid;
    error = mb->tag_get_data(v2hf_cid, &fid_verts[0], 1, &cur_cid);
    if (MB_SUCCESS != error) return error;

    bool found = false;

    if (cur_cid != 0){
        EntityHandle Stkcells[MAXSIZE], trackcells[MAXSIZE];
        for (int i = 0; i< MAXSIZE; i++)
          {
            Stkcells[i] = 0;
            trackcells[i] = 0;
          }
        int Stksize = 0, count = -1;
        Stkcells[0] = cur_cid;

        while (Stksize >= 0 ){
            cur_cid = Stkcells[Stksize];
            Stksize -= 1 ;
            count += 1;
            trackcells[count] = cur_cid;

            std::vector<EntityHandle> conn;
            error = mb->get_connectivity(&cur_cid, 1, conn);
            if (MB_SUCCESS != error) return error;

	    // Local id of fid_verts[0] in the cell
	    int lv0;
	    for (int i = 0; i< nvpc; ++i){
		if (conn[i] == fid_verts[0])
		  {
		    lv0 = i;
		  }
	      };

	    int nhf_thisv = lConnMap3D[index].v2hf_num[lv0];

            // Search each half-face to match input face
            for(int i = 0; i < nhf_thisv; ++i){
                int lfid = lConnMap3D[index].v2hf[lv0][i];
                int nv_curF = lConnMap3D[index].hf2v_num[lfid];
                if (nv_curF != nvF)
                  continue;

                // Connectivity of the current half-face

                std::vector<EntityHandle> vthisface(nvF);
                for(int l = 0; l < nvF; ++l){
                    int ind = lConnMap3D[index].hf2v[lfid][l];
                    vthisface[l] = conn[ind];
                  };


                // Match this half-face with input fid
                int direct,offset;
                bool they_match = CN::ConnectivityMatch(&vthisface[0],&fid_verts[0],nvF,direct,offset);

                if (they_match){
                    found = true;
                    cid[0] = cur_cid;
                    lid[0] = lfid;

                    break;
                  }
              }

	    // Add other cells that are incident on fid_verts[0]
	    std::vector<EntityHandle> sib_cids(nfpc);
	    error = mb->tag_get_data(sibhfs_cid, &cur_cid, 1, &sib_cids[0]);
	    if (MB_SUCCESS != error) return error;

	    // Add new cells into the stack
	    EntityHandle ngb;
	    for (int i = 0; i < nhf_thisv; ++i){
		int ind = lConnMap3D[index].v2hf[lv0][i];
		ngb = sib_cids[ind];

		if (ngb) {

		    bool found_ent = find_match_in_array(ngb, trackcells, count);

		    if (!found_ent){
			Stksize += 1;
			Stkcells[Stksize] = ngb;
		      }
		  }
	      }
	  }
      }
    return found;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ErrorCode HalfFacetRep::get_neighbor_adjacencies_3d( EntityHandle cid, std::vector<EntityHandle> &adjents)
  {
    ErrorCode error;
    int index = get_index_from_type(cid);
    int nfpc = lConnMap3D[index].num_faces_in_cell;
    
    if (cid != 0 ){

      std::vector<EntityHandle> sibcids(nfpc);
      std::vector<int> siblfids(nfpc);
      error = mb->tag_get_data(sibhfs_cid, &cid, 1, &sibcids[0]);
      if (MB_SUCCESS != error) return error;      
      error = mb->tag_get_data(sibhfs_lfid, &cid, 1, &siblfids[0]);
      if (MB_SUCCESS != error) return error;

      for (int lfid = 0; lfid < nfpc; ++lfid){      
	if (sibcids[lfid] != 0) 
	  adjents.push_back(sibcids[lfid]);
      }    
    }

    return MB_SUCCESS; 
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  bool HalfFacetRep::find_match_in_array(EntityHandle ent, EntityHandle *ent_list, int count, bool get_index, int *index)
  {
    bool found = false;
   // std::cout<<"count = "<<count<<std::endl;
    for (int i = 0; i<= count; i++)
      {
     //   std::cout<<"For i = "<<i<<": ent-ent_list[i] = "<<(int)(ent-ent_list[i])<<std::endl;
        if (!((int)(ent - ent_list[i])))
          {
            found = true;
            if (get_index)
              *index = i;
            break;
          }
      }

    return found;
  }

} // namespace moab

