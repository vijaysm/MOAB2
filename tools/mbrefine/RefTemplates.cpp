#ifdef WIN32
#pragma warning (disable : 4786)
#endif 

#include "moab/NestedRefine.hpp"
#include "moab/HalfFacetRep.hpp"
#include "moab/ReadUtilIface.hpp"
#include <iostream>
#include <assert.h>
#include <vector>


namespace moab{


  ///////////////////////////////////////////////
  /*ErrorCode UniRefine::initialize(){
    
    error = ahf.initialize();
    return MB_SUCCESS;
  }
  
    ErrorCode UniRefine::deinitialize(){
    ErrorCode error;
    error = ahf->deinitialize();
    }*/


  ////////////////////////////////////////
  int UniRefine::get_index_from_type(EntityHandle ent){
   
    EntityType type = mb->type_from_handle(ent);
    int index;
    if (type == MBEDGE)
      index = 0;
    else if (type == MBTRI)
      index = 1;
    else if (type == MBQUAD)
      index = 2;
    return index;   
  }

  const UniRefine::LevelTemplates UniRefine::refPatterns[3][MAX_LEVEL] =
  {
    
    //EDGE
    {{1,0,1,2,{2,2},{},{{0,2},{2,1}},{}},
     {3,0,3,4,{2,4},{},{{0,2},{2,3},{3,4},{4,1}},{}},
     {5,0,5,6,{2,6},{},{{0,2},{2,3},{3,4},{4,5},{5,6},{6,1}},{}}},

    //TRI
    {{1,0,3,4,{3,5},{},{{0,3,4},{3,5,4},{3,1,5},{4,5,2}},{{3},{5},{4},{}}},

     {3,3,12,16,{3,14},{},{{0,3,6},{3,7,6},{3,4,7},{4,8,7},{4,5,8},{5,9,8},{5,1,9},{6,7,10},{7,11,10},{7,8,11},{8,12,11},{8,9,12},{10,11,13},{11,14,13},{11,12,14},{13,14,2}},{{3,4,5},{9,12,14},{13,10,6},{7,8,11}}},

     {5,10,25,36,{3,27},{},{{0,3,8},{3,9,8},{3,4,9},{4,10,9},{4,5,10},{5,11,10},{5,6,11},{6,12,11},{6,7,12},{7,13,12},{7,1,13},{8,9,14},{9,15,14},{9,10,15},{10,16,15},{10,11,16},{11,17,16},{11,12,17},{12,18,17},{12,13,18},{14,15,19},{15,20,19},{15,16,20},{16,21,20},{16,17,21},{17,22,21},{17,18,22},{19,20,23},{20,24,23},{20,21,24},{21,25,24},{21,22,25},{23,24,26},{24,27,26},{24,25,27},{26,27,2}},{{3,4,5,6,7},{13,18,22,25,27},{26,23,19,14,8},{9,10,11,12,15,16,17,20,21,24}}}},

    //QUAD
    {{1,1,5,4,{4,8},{},{{0,4,5,6},{4,1,7,6},{5,6,3,8},{6,7,8,2}},{{4},{7},{8},{5},{6}}},

     {3,9,21,16,{4,24},{},{{0,4,8,7},{4,5,9,8},{5,6,10,9},{6,1,11,10},{7,8,13,12},{8,9,14,13},{9,10,15,14},{10,11,16,15},{12,13,18,17},{13,14,19,18},{14,15,20,19},{15,16,21,20},{17,18,22,3},{18,19,23,22},{19,20,24,23},{20,21,2,24}},{{4,5,6},{11,16,21},{24,23,22},{17,12,7},{8,9,10,13,14,15,18,19,20}}},

     {5,25,45,36,{4,48},{},{{0,4,9,10},{4,5,11,10},{5,6,12,11},{6,7,13,12},{7,8,14,13},{8,1,15,14},{9,10,17,16},{10,11,18,17},{11,12,19,18},{12,13,20,19},{13,14,21,20},{14,15,22,21},{16,17,24,23},{17,18,25,24},{18,19,26,25},{19,20,27,26},{20,21,28,27},{21,22,29,28},{23,24,31,30},{24,25,32,31},{25,26,33,32},{26,27,34,33},{27,28,35,34},{28,29,36,35},{30,31,38,37},{31,32,39,38},{32,33,40,39},{33,34,41,40},{34,35,42,41},{35,36,43,42},{37,38,44,3},{38,39,45,44},{39,40,46,45},{40,41,47,46},{41,42,48,47},{42,43,2,48}},{{4,5,6,7,8},{15,22,29,36,43},{48,47,46,45,44},{37,30,23,16,9},{10,11,12,13,14,17,18,19,20,21,24,25,26,27,28,31,32,33,34,35,38,39,40,41,42}}}}

    //
  };
  


  /////////////////////////////////////////////////////////
  ErrorCode UniRefine::get_total_new_verts_and_ents(EntityHandle edg, int nedges, int level, METHOD method, int *newverts, int *newsubents)
  {
    int index = get_index_from_type(edg);
    if (method == MULTI_LEVEL)
      {
        int nve = refPatterns[index][level].num_new_verts_per_edge;
        int nents = refPatterns[index][level].num_new_ents;
	
	newverts[0] = nve*nedges;
	newsubents[0] = nents*nedges;
      }
    else if (method == ALL_LEVELS)
      {
	for (int i = 0; i < level; i++)
	  {
	    int nve = refPatterns[index][i].num_new_verts_per_edge;
	    newverts[0] += nve*nedges;
	    
	    int nents = refPatterns[index][i].num_new_ents;
	    newsubents[0] += nents*nedges; 
	  };
      } 
    return MB_SUCCESS;
  }
  /////////////////////////////////////////////////////////
    ErrorCode UniRefine::get_total_new_verts_and_ents(EntityHandle face, int nedges, int nfaces, int level, METHOD method, int *newverts, int *newsubents)
  {

    int index = get_index_from_type(face);
    if (method == MULTI_LEVEL)
      {
        int nve = refPatterns[index][level].num_new_verts_per_edge;
        int nvf = refPatterns[index][level].num_new_verts_per_face;
        int nents = refPatterns[index][level].num_new_ents;

	newverts[0] = nve*nedges + nvf*nfaces;
	newsubents[0] = nents*nfaces;
      }
    else if (method == ALL_LEVELS)
      {
        for (int i = 0; i < level; i++)
          {
            int nve = refPatterns[index][i].num_new_verts_per_edge;
            int nvf = refPatterns[index][i].num_new_verts_per_face;
            newverts[0] += nve*nedges + nvf*nfaces;

	    int nents = refPatterns[index][i].num_new_ents;
	    newsubents[0] += nents*nfaces;
	  };
      }
    return MB_SUCCESS;
  }
    /////////////////////////////////////////////////////////
    ErrorCode UniRefine::get_total_new_verts_and_ents(EntityHandle cell, int nedges, int nfaces, int ncells, int level, METHOD method, int *newverts, int *newsubents)
  {

    int index = get_index_from_type(cell);
    if (method == MULTI_LEVEL)
      {
        int nve = refPatterns[index][level].num_new_verts_per_edge;
        int nvf = refPatterns[index][level].num_new_verts_per_face;
        int nvc = refPatterns[index][level],num_new_verts_per_cell;
        int nents = refPatterns[index][level].num_new_ents;

	newverts[0] = nve*nedges + nvf*nfaces + nvc*ncells;
	newsubents[0] = nents*ncells;
      }
    else if (method == ALL_LEVELS)
      {
        for (int i = 0; i < level; i++)
          {
            int nve = refPatterns[index][i].num_new_verts_per_edge;
            int nvf = refPatterns[index][i].num_new_verts_per_face;
            int nvc = refPatterns[index][i].num_new_verts_per_cell;
            newverts[0] += nve*nedges + nvf*nfaces+ nvc*ncells;

	    int nents = refPatterns[index][i].num_new_ents;
	    newsubents[0] += nents*ncells;
	  };
      }
    return MB_SUCCESS;
  }

  //////////////////////////////////////////////////////////
  ErrorCode UniRefine::uniform_refinement( Range &ents, int level, METHOD method){
    ErrorCode error;

    Range verts = ents.subset_by_dimension(0);
    Range edges = ents.subset_by_dimension(1);
    Range faces = ents.subset_by_dimension(2);
    Range cells = ents.subset_by_dimension(3);

    MESHTYPE mesh_type = get_mesh_type(verts.size(), edges.size(), faces.size(), cells.size());
    //if parallel, get skin entities for this part mesh and add the explicit skin entities to appropriate range.

    if (cells.size()==0)
      {
        error = uniform_refinement_mixed_2d(verts, edges, faces, level, method);
      }
    else
      {
        error = uniform_refinement_mixed_3d(verts, edges, faces, cells, level, method);
      }

  }

  ////////////////////////////////////////////////////////
  ErrorCode UniRefine::uniform_refinement_mixed_2d(Range &verts, Range &edges, Range &faces, int level, METHOD method)
  {
    ErrorCode error;
    EntityHandle first_edge = *edges.begin();
    int nExpEdges = edges.size();
    EntityHandle first_face = *faces.begin();
    int nFaces = faces.size();
    
    level = level -1 ; // The index for levels in one less.

    //Step 1: Get total no. of new vertices and new entities that have to be created for each entity type.
    int newvertsE = 0 , newsubEdges = 0 ;
    error = get_total_new_verts_and_ents(first_edge, nExpEdges, level, method, &newvertsE, &newsubEdges);
    if (MB_SUCCESS != error) return error;

    int newvertsF = 0, newsubFaces = 0;
    int totalEdges = find_total_edges_2d(faces, sibhes_fid, sibhes_leid);
    error = get_total_new_verts_and_ents(first_face, totalEdges, nFaces, level, method, &newvertsF, &newsubFaces);
 
    // Step 2: Pre-allocate memory and create handles for new nodes and new entities for each entity type.   
    ReadUtilIface *read_iface;
    error = mb->query_interface(read_iface);

    //Vertices
    std::vector<double *> coords;
    EntityHandle start_nverts;
    error = read_iface->get_node_coords(3, newvertsF, 0, start_nverts, coords);

    // Edges
    EntityHandle start_nedges, *econnect;
    error = read_iface->get_element_connect(newsubEdges, 2, MBEDGE, 0, start_nedges, econnect);

    //Faces
    EntityHandle start_nfaces, *fconnect;
    EntityType type = mb->type_from_handle(first_face);
    int nvpf = local_maps_2d(first_face);
    error = read_iface->get_element_connect(newsubFaces, nvpf, type, 0, start_nfaces, fconnect); 
       
    if (method == MULTI_LEVEL){
      EntityHandle vert_Hbnds[2], edge_Hbnds[2], face_Hbnds[2];
      verts_Hbnds[0] = start_nverts; edge_Hbnds[0] = start_nedges; face_Hbnds[0] = start_faces;
      error = uniform_refinement_mixed_single_level_2d(edges, faces, level, vert_Hbnds, coords, edge_Hbnds, econnect, face_Hbnds, fconnect);
      }
    else if (method == ALL_LEVELS)
      {
        EntityHandle *vert_Hbnds = new EntityHandle(level*2);
        EntityHandle *edge_Hbnds = new EntityHandle(level*2);
        EntityHandle *face_Hbnds = new EntityHandle(level*2);
        EntityHandle vert_bnds[2], edge_bnds[2], face_bnds[2];

        vert_Hbnds[0] = start_nverts; edge_Hbnds[0] = start_nedges; face_Hbnds[0] = start_nfaces;

	for (int l = 0; l <= level; l++){
	    vert_bnds[0] = vert_Hbnds[2*l]; edge_bnds[0] = edge_Hbnds[2*l]; face_bnds[0] = face_Hbnds[2*l];

	    error = uniform_refinement_mixed_single_level_2d(edges, faces, level, vert_bnds, coords, edge_bnds, econnect, face_bnds, fconnect);

	    vert_Hbnds[2*l+1] = vert_bnds[1]; edge_Hbnds[2*l+1] = edge_bnds[1]; face_Hbnds[2*l+1] = face_bnds[1];
	    vert_Hbnds[2*(l+1)] = vert_bnds[1]+1; edge_Hbnds[2*(l+1)] = edge_bnds[1]+1; face_Hbnds[2*(l+1)] = face_bnds[1]+1;
	  }
      }
  }

  ////////////////////////////////////////////////////////////
  ErrorCode UniRefine::uniform_refinement_mixed_single_level_2d(Range &edges, Range &faces, int level, EntityHandle vert_bnds[2], std::vector<double *> coords, EntityHandle edge_bnds[2], EntityHandle *econnect, EntityHandle face_bnds[2], EntityHandle *fconnect )
  {
    ErrorCode error;
    EntityHandle first_edge = *edges.begin();
    EntityHandle first_face = *faces.begin();
    int indE = get_index_from_type(first_edge);
    int indF = get_index_from_type(first_face);
    EntityHandle start_nverts = verts_bnds[0];

    // Use a list of vertices to avoid introducing duplicates.
    // This list would be used for parallel communication later.

    int nFaces = faces.size();
    int newvF = refPatterns[indF][level].total_new_verts;
    std::vector<EntityHandle> trackvertsF(newvF*nFaces);
    for (int i=0; i< newvF*nFaces; i++)
      trackverts[i] = 0;

    int nv = refPatterns[indE][level].num_new_verts_per_edge;
    std::vector<EntityHandle> vbuf_Edg(nv+2);
    int vcount = 0, ecount = 0, fcount = 0;

    //Loop over edges
    for (Range::iterator eid = edges.begin(); eid != edges.end(); eid++){
      const EntityHandle *conn;
      int num_nodes = 0;
      error = mb->get_connectivity(*eid, conn, num_nodes);  
      if (MB_SUCCESS != error) return error;

      // Fill up the vbuf for this edge
      vbuf_Edg[0] = conn[0];
      vbuf_Edg[1] = conn[1];
      for (int i=0; i<nv; i++){
        vbuf_Edg[i+2] = start_nverts+vcount;
        vcount += 1;
      };  
      
      // Create connectivities of new sub-edges
      int nents = refPatterns[indE][level].num_new_ents;
      for (int i =0; i< nents; i++){
        int id1 = refPatterns[indE][level].new_entsConn[i][0];
        int id2 = refPatterns[indE][level].new_entsConn[i][1];

	econnect[2*ecount] = vbuf_Edg[id1];
	econnect[2*ecount+1] = vbuf_Edg[id2];
	ecount += 1;
      };

      // Find incident faces on the edges, add the vertex handles
      std::vector<EntityHandle> fids;
      std::vector<int> leids;
      error = get_upward_incidences_2d(*eid, faces, fids, leids);

      for (int i= 0; i< (int)fids.size(); i++){
	EntityHandle curfid = fids[i];
	int curlid = leids[i];
	for (int j = 0; j<nv; j++){
	  int id = refPatterns[indF][level].vert_identify[curlid][j];
	  if (!trackvertsF[newvF*(curfid-first_face)+id-nvpf])
	      trackvertsF[newvF*(curfid-first_face)+id-nvpf] = vbuf_Edg[j+2];
	  }
	}
      };

    // The entity handle of the last subedge.
    edge_bnds[1] = edge_bnds[0]+ecount-1;

    int nvpf = local_maps_2d(first_face);
    nv = refPatterns[indF][level].total_new_verts;
    std::vector<EntityHandle> vbuf_Face(nv+nvpf);

    // Loop over the faces
    for (Range::iterator fid = faces.begin(); fid != faces.end(); fid++){
      const EntityHandle *conn;
      int num_nodes = 0;
      error = mb->get_connectivity(*fid, conn, num_nodes);  
      
      // Get a list of vertices
      for (int i = 0; i < nvpf; i++)
        vbuf_Face[i] = conn[i];
      for (int j = 0; j < nv; j++){
	if (trackvertsF[newvF*(*fid-first_face)+j]!=0)
	  vbuf_Face[j+nvpf] = trackvertsF[newvF*(*fid-first_face)+j];
	else
	  {
	  vbuf_Face[j+nvpf] = start_nverts+vcount;
	  trackvertsF[newvF*(*fid-first_face)+j] = start_nverts+vcount;
	  vcount += 1;
	  };
      };

      // Get connectivities of new sub-entities
      int nents = refPatterns[indF][level].num_new_ents;
      for (int i = 0; i < nents; i++){
	for (int j = 0; j < nvpf; j++){
	  int id = refPatterns[indF][level].new_entsConn[i][j];
	  fconnect[nvpf*fcount+j] = vbuf_Face[id];
	  fcount += 1;
	};
      };

      // Compute the vertex coordinates
      for (int i = 0; i < nv; i++){
          //double xi = refPatterns[indF][level].vert_params[i][0];
          //double eta = refPatterns[indF][level].vert_params[i][1];
          // Call the function to high-order point projection
          EntityHandle vert = trackvertsF(newvF*(*fid-first_face)+i);
        coords[3*(vert-start_nverts)] = 0;
        coords[3*(vert-start_nverts)+1] = 0;
        coords[3*(vert-start_nverts)+2] = 0;
      };

      face_bnds[1] = face_bnds[0]+fcount-1;
      vert_bnds[1] = vert_bnds[0]+vcount-1 ;
    }; 
  }

  ////////////////////////////////////////////////////////
  ErrorCode UniRefine::uniform_refinement_mixed_3d(Range &verts, Range &edges, Range &faces, Range &cells, int level, METHOD method)
  {
    ErrorCode error;
    EntityHandle first_edge = *edges.begin();
    int nExpEdges = edges.size();
    EntityHandle first_face = *faces.begin();
    int nExpFaces = faces.size();
    EntityHandle first_cell = *cells.begin();
    int nCells = cells.size();

    level = level -1 ; // The index for levels in one less.

    //Step 1: Get total no. of new vertices and new entities that have to be created for each entity type.
    int newvertsE = 0 , newsubEdges = 0 ;
    error = get_total_new_verts_and_ents(first_edge, nExpEdges, level, method, &newvertsE, &newsubEdges);
    if (MB_SUCCESS != error) return error;

    int newvertsF = 0, newsubFaces = 0;
    int impEdgesS = find_total_edges_2d(faces);
    error = get_total_new_verts_and_ents(first_face, impEdgesS, nExpFaces, level, method, &newvertsF, &newsubFaces);
    if (MB_SUCCESS != error) return error;

    int newvertsC = 0; newsubCells = 0;
    int impEdgesV = find_total_edges_3d(cells);
    int impFacesV = find_total_faces_3d(cells);
    error = get_total_new_verts_and_ents(first_cell,impEdgesV, impFacesV, nCells, level, method, &newvertsC, &newsubCells);
    if (MB_SUCCESS != error) return error;

    // Step 2: Pre-allocate memory and create handles for new nodes and new entities for each entity type.
    ReadUtilIface *read_iface;
    error = mb->query_interface(read_iface);

    //Vertices
    std::vector<double *> coords;
    EntityHandle start_nverts;
    error = read_iface->get_node_coords(3, newvertsF, 0, start_nverts, coords);

    // Edges
    EntityHandle start_nedges, *econnect;
    error = read_iface->get_element_connect(newsubEdges, 2, MBEDGE, 0, start_nedges, econnect);
    if (MB_SUCCESS != error) return error;

    //Faces
    EntityHandle start_nfaces, *fconnect;
    EntityType ftype = mb->type_from_handle(first_face);
    int nvpf = local_maps_2d(first_face);
    error = read_iface->get_element_connect(newsubFaces, nvpf, ftype, 0, start_nfaces, fconnect);
    if (MB_SUCCESS != error) return error;

    //Cells
    EntityHandle start_ncells, *cconnect;
    EntityType ctype = mb->type_from_handle(first_cell);
    int index = get_index_from_type(first_cell);
    int nvpc = lConnMap3D[index].num_verts_in_cell;
    error = read_iface->get_element_connect(newsubCells, nvpc, ctype, 0, start_ncells, cconnect);
    if (MB_SUCCESS != error) return error;

    // Step 3: Call the refinement routine
    if (method == MULTI_LEVEL){
      EntityHandle vert_Hbnds[2], edge_Hbnds[2], face_Hbnds[2], cell_Hbnds[2];
      verts_Hbnds[0] = start_nverts; edge_Hbnds[0] = start_nedges; face_Hbnds[0] = start_nfaces; cell_Hbnds[0] = start_ncells;
      error = uniform_refinement_mixed_single_level_3d(edges, faces, cells, level, vert_Hbnds, coords, edge_Hbnds, econnect, face_Hbnds, fconnect, cell_Hbnds, cconnect);
      }
    else if (method == ALL_LEVELS)
      {
        EntityHandle *vert_Hbnds = new EntityHandle((level+1)*2);
        EntityHandle *edge_Hbnds = new EntityHandle((level+1)*2);
        EntityHandle *face_Hbnds = new EntityHandle((level+1)*2);
        EntityHandle *cell_Hbnds = new EntityHandle((level+1)*2);
        EntityHandle vert_bnds[2], edge_bnds[2], face_bnds[2], cell_bnds[2];

        vert_Hbnds[0] = start_nverts; edge_Hbnds[0] = start_nedges; face_Hbnds[0] = start_nfaces; cell_Hbnds[0] = start_ncells;

	for (int l = 0; l <= level; l++){
	    vert_bnds[0] = vert_Hbnds[2*l]; edge_bnds[0] = edge_Hbnds[2*l];
	    face_bnds[0] = face_Hbnds[2*l]; cell_bnds[0] = cell_Hbnds[2*l];

	    error = uniform_refinement_mixed_single_level_3d(edges, faces, cells, level, vert_bnds, coords, edge_bnds, econnect, face_bnds, fconnect, cell_bnds, cconnect);

	    vert_Hbnds[2*l+1] = vert_bnds[1]; edge_Hbnds[2*l+1] = edge_bnds[1];
	    face_Hbnds[2*l+1] = face_bnds[1]; cell_Hbnds[2*l+1] = cell_bnds[1];
	    vert_Hbnds[2*(l+1)] = vert_bnds[1]+1; edge_Hbnds[2*(l+1)] = edge_bnds[1]+1;
	    face_Hbnds[2*(l+1)] = face_bnds[1]+1; cell_Hbnds[2*(l+1)] = cell_bnds[1]+1;
	  }
      }
    // Write out or create new Ranges
  }

  ////////////////////////////////////////////////////////////
  ErrorCode UniRefine::uniform_refinement_mixed_single_level_3d(Range &edges, Range &faces, Range &cells, int level, EntityHandle vert_bnds[2], std::vector<double *> coords, EntityHandle edge_bnds[2], EntityHandle *econnect, EntityHandle face_bnds[2], EntityHandle *fconnect , EntityHandle cell_bnds[2], EntityHandle *cconnect)
  {
    ErrorCode error;
    EntityHandle first_edge = *edges.begin();
    EntityHandle first_face = *faces.begin();
    EntityHandle first_cell = *cells.begin();
    int indE = get_index_from_type(first_edge);
    int indF = get_index_from_type(first_face);
    int indC = get_index_from_type(first_cell);

    EntityHandle start_nverts = verts_bnds[0];

    // Use a list of vertices to avoid introducing duplicates.
    // This list would be used for parallel communication later.

    int nCells = cells.size();
    int newvC = refPatterns[indC][level].total_new_verts;
    std::vector<EntityHandle> trackvertsC(newvC*nCells);
    for (int i = 0; i< newvC*nCells; i++)
      trackvertsC[i] = 0;

    int nFaces = faces.size();
    int newvF = refPatterns[indF][level].total_new_verts;
    std::vector<EntityHandle> trackvertsF(newvF*nFaces);
    for (int i = 0; i< newvF*nFaces; i++)
      trackvertsF[i] = 0;

    // Start refinement
    int vcount = 0, ecount = 0, fcount = 0, ccount = 0;

    // First Loop: Over all explicit edges
    int nvE = refPatterns[indE][level].total_new_verts;
    std::vector<EntityHandle> vbuf_Edg(nvE+2);

    for (Range::iterator eid = edges.begin(); eid != edges.end(); eid++){
      const EntityHandle *conn;
      int num_nodes = 0;
      error = mb->get_connectivity(*eid, conn, num_nodes);
      if (MB_SUCCESS != error) return error;

      // Get a list of vertex handles
      vbuf_Edg[0] = conn[0];
      vbuf_Edg[1] = conn[1];
      for (int i=0; i<nvE; i++){
        vbuf_Edg[i+2] = start_nverts+vcount;
        vcount += 1;
      };

      // Create connectivities of new sub-edges
      int nents = refPatterns[indE][level].num_new_ents;
      for (int i =0; i< nents; i++){
        int id1 = refPatterns[indE][level].new_entsConn[i][0];
        int id2 = refPatterns[indE][level].new_entsConn[i][1];

	econnect[2*ecount] = vbuf_Edg[id1];
	econnect[2*ecount+1] = vbuf_Edg[id2];
	ecount += 1;
      };

      // Find incident faces on the edges.
      std::vector<EntityHandle> fids;
      std::vector<int> leids;
      error = get_upward_incidences_2d(*eid, faces, fids, leids);
      if (MB_SUCCESS != error) return error;

      for (int i= 0; i< (int)fids.size(); i++){
        EntityHandle curfid = fids[i];
        int curlid = leids[i];
        for (int j = 0; j<nvE; j++){
          int id = refPatterns[indF][level].vert_identify[curlid][j];
          if (!trackvertsF[newvF*(curfid-first_face)+id-nvpf])
              trackvertsF[newvF*(curfid-first_face)+id-nvpf] = vbuf_Edg[j+2];
        };
      };

      // Find incident cells on the edges
      std::vector<EntityHandle> cids;
      std::vector<int> leids;
      error = get_upward_incidences_3d(*eid, cells, cids, leids);
      if (MB_SUCCESS != error) return error;

      for (int i= 0; i< (int)cids.size(); i++){
        EntityHandle curcid = cids[i];
        int curlid = leids[i];
        for (int j = 0; j<nvE; j++){
          int id = refPatterns[indC][level].vert_identify[curlid][j];
          if (trackvertsC[newvC*(curcid-first_cell)+id-nvpc] != 0)
              trackvertsC[newvC*(curcid-first_cell)+id-nvpc] = vbuf_Edg[j+2];
        };
      };

    };

    // The entity handle of the last subedge.
    edge_bnds[1] = edge_bnds[0]+ecount-1;

    // Second Loop: Over all explicit faces
    int nvpf = local_maps_2d(first_face);
    int nvF = refPatterns[indF][level].total_new_verts;
    std::vector<EntityHandle> vbuf_Face(nvF+nvpf);

    for (Range::iterator fid = faces.begin(); fid != faces.end(); fid++){
      const EntityHandle *conn;
      int num_nodes = 0;
      error = mb->get_connectivity(*fid, conn, num_nodes);

      // Get a list of vertices
      for (int i = 0; i < nvpf; i++)
        vbuf_Face[i] = conn[i];
      for (int j = 0; j < nvF; j++){
        if (trackvertsF[newvF*(*fid-first_face)+j]!=0)
          vbuf_Face[j+nvpf] = trackvertsF[newvF*(*fid-first_face)+j];
        else
          {
          vbuf_Face[j+nvpf] = start_nverts + vcount;
          trackvertsF[newvF*(*fid-first_face)+j] = start_nverts+vcount;
          vcount += 1;
          }
      }

      // Get connectivities of new sub-entities
      int nents = refPatterns[indF][level].num_new_ents;
      for (int i = 0; i < nents; i++){
          for (int j = 0; j < nvpf; j++){
              int id = refPatterns[indF][level].new_entsConn[i][j];
              fconnect[nvpf*fcount+j] = vbuf_Face[id];
            }
        fcount += 1;
      }

      // Find incident cells on the face
      std::vector<EntityHandle> cids;
      std::vector<int> lfids;
      error = get_upward_incidences_3d(*fid, cells, cids, lfids);
      if (MB_SUCCESS != error) return error;

      for (int i= 0; i< (int)cids.size(); i++){
        EntityHandle curcid = cids[i];
        int curlid = lfids[i];
        for (int j = 0; j<nvF; j++){
          int id = refPatterns[indC][level].vert_identify[curlid][j];
          if (trackvertsC[newvC*(curcid-first_cell)+id-nvpc] != 0)
              trackvertsC[newvC*(curcid-first_cell)+id-nvpc] = vbuf_Face[j+nvpf];
        };
      };


      // Compute the vertex coordinates
      for (int i = 0; i < nv; i++){
          //double xi = refPatterns[indF][level].vert_params[i][0];
          //double eta = refPatterns[indF][level].vert_params[i][1];
          // Call the function to high-order point projection
          EntityHandle vert = trackvertsF(newvF*(*fid-first_face)+i);
          coords[3*(vert-start_nverts)] = 0;
          coords[3*(vert-start_nverts)+1] = 0;
          coords[3*(vert-start_nverts)+2] = 0;
      };

      face_bnds[1] = face_bnds[0] + fcount;

      // Third Loop: Over all cells
      int index = get_index_from_type(first_cell);
      int nvpc = lConnMap3D[index].num_verts_in_cell;
      int nvC = refPatterns[indC][level].total_new_verts;
      std::vector<EntityHandle> vbuf_Cell(nvC+nvpc);

      for (Range::iterator cid = cells.begin(); cid != cells.end(); cid++){
        const EntityHandle *conn;
        int num_nodes = 0;
        error = mb->get_connectivity(*cid, conn, num_nodes);

        // Fill out vbuf for this cell
        for (int i = 0; i < nvpc; i++)
          vbuf_Face[i] = conn[i];
        for (int j = 0; j < nvC; j++){
          if (trackvertsC[newvC*(*cid-first_cell)+j]!=0)
            vbuf_Cell[j+nvpc] = trackvertsC[newvC*(*cid-first_cell)+j];
          else
            {
            vbuf_Cell[j+nvpc] = start_nverts + vcount;
            trackvertsC[newvC*(*cid-first_cell)+j] = start_nverts+vcount;
            vcount += 1;
            }
        }

        // Get connectivities of new sub-entities
        int nents = refPatterns[indC][level].num_new_ents;
        for (int i = 0; i < nents; i++){
          for (int j = 0; j < nvpc; j++){
            int id = refPatterns[indC][level].new_entsConn[i][j];
            fconnect[nvpc*ccount+j] = vbuf_Cell[id];
            ccount += 1;
          };
        };
  /*
        // Compute the vertex coordinates
        for (int i = 0; i < nv; i++){
            //double xi = refPatterns[indF][level].vert_params[i][0];
            //double eta = refPatterns[indF][level].vert_params[i][1];
            // Call the function to high-order point projection
            EntityHandle vert = trackvertsC(newvF*(*fid-first_face)+i);
          coords[3*(vert-start_nverts)] = 0;
          coords[3*(vert-start_nverts)+1] = 0;
          coords[3*(vert-start_nverts)+2] = 0;*/
        };

        cell_bnds[1] = cell_bnds[0] + ccount;
        vert_bnds[1] = vert_bnds[0] + vcount;

    };
  }

}//namesapce moab
  
