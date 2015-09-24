
#include <vector>
#include <utility>
#include <iostream>
using namespace std;
using namespace moab;

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"

Interface* gMB = 0;

class EntityCount
{
public:
  unsigned int node;
  unsigned int edge;
  unsigned int quad;
  unsigned int tri;
  unsigned int hex;
  unsigned int tet;

  EntityCount();

  ErrorCode get_counts();
  ErrorCode create_adjacencies(Range &entities, int adj_dim);
  void copy_counts(EntityCount &count);
  void print();
};

EntityCount::EntityCount()
{
  node = 0;
  edge = 0;
  quad = 0;
  tri = 0;
  hex = 0;
  tet = 0;
}

void EntityCount::copy_counts(EntityCount &count)
{
  node = count.node;
  edge = count.edge;
  quad = count.quad;
  tri = count.tri;
  hex = count.hex;
  tet = count.tet;
}

ErrorCode EntityCount::get_counts()
{
  Range entities;
  int do_create = edge == 0;

  if(gMB->get_entities_by_type(0, MBVERTEX, entities) != MB_SUCCESS)
    return MB_FAILURE;
  node = entities.size();

  entities.clear();
  if(gMB->get_entities_by_type(0, MBHEX, entities) != MB_SUCCESS)
    return MB_FAILURE;
  hex = entities.size();
  if(hex > 0 && do_create)
  {
    if(create_adjacencies(entities, 2) != MB_SUCCESS)
      return MB_FAILURE;
    if(create_adjacencies(entities, 1) != MB_SUCCESS)
      return MB_FAILURE;
  }

  entities.clear();
  if(gMB->get_entities_by_type(0, MBQUAD, entities) != MB_SUCCESS)
    return MB_FAILURE;
  quad = entities.size();
  if(quad > 0 && do_create)
  {

    if(create_adjacencies(entities, 1) != MB_SUCCESS)
      return MB_FAILURE;
  }

  entities.clear();
  if(gMB->get_entities_by_type(0, MBTET, entities) != MB_SUCCESS)
    return MB_FAILURE;
  tet = entities.size();
  if(tet > 0 && do_create)
  {
    if(create_adjacencies(entities, 2) != MB_SUCCESS)
      return MB_FAILURE;
    if(create_adjacencies(entities, 1) != MB_SUCCESS)
      return MB_FAILURE;
  }

  entities.clear();
  if(gMB->get_entities_by_type(0, MBTRI, entities) != MB_SUCCESS)
    return MB_FAILURE;
  tri = entities.size();
  if(tri > 0 && do_create)
  {
    if(create_adjacencies(entities, 1) != MB_SUCCESS)
      return MB_FAILURE;
  }

  entities.clear();
  if(gMB->get_entities_by_type(0, MBEDGE, entities) != MB_SUCCESS)
    return MB_FAILURE;
  edge = entities.size();

  return MB_SUCCESS;
}

ErrorCode EntityCount::create_adjacencies(Range &entities, int adj_dim)
{
  ErrorCode result;
  Range::iterator iter;
  std::vector<EntityHandle> adjacencies;
  
  for (iter = entities.begin(); iter != entities.end(); ++iter)
  {
    result = gMB->get_adjacencies(&*iter, 1, adj_dim, true, adjacencies);
    if(result != MB_SUCCESS)
      break;
  }

  return result;
}

void EntityCount::print()
{
  std::cout << "   Vertices: " << node << std::endl;
  cout << "   Edges: " << edge << endl;
  if(quad > 0)
    cout << "   Quad Elements: " << quad << endl;
  if(hex > 0)
    cout << "   Hex Elements: " << hex << endl;
  if(tri > 0)
    cout << "   Tri Elements: " << tri << endl;
  if(tet > 0)
    cout << "   Tet Elements: " << tet << endl;
}

bool points_are_coincident(const double *first, const double *second)
{
  double diff[3];
  diff[0] = first[0] - second[0];
  diff[1] = first[1] - second[1];
  diff[2] = first[2] - second[2];
  double length = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];

  if(fabs(length) < .001)
    return true;

  return false;
}

// dumb n^2 coincident node algorithm
ErrorCode find_coincident_nodes(Range vertices,
  std::vector< std::pair<EntityHandle,EntityHandle> > &coin_nodes)
{
  double first_coords[3], second_coords[3];
  Range::iterator iter, jter;
  std::pair<EntityHandle, EntityHandle> coincident_pair;
  ErrorCode result;

  for (iter = vertices.begin(); iter != vertices.end(); ++iter)
  {
    result = gMB->get_coords(&*iter, 1, first_coords);
    if (result != MB_SUCCESS)
      return result;

    for (jter = iter; jter != vertices.end(); ++jter)
    {
      if (*iter != *jter)
      {
        result = gMB->get_coords(&*jter, 1, second_coords);
        if (result != MB_SUCCESS)
          return result;

        if(points_are_coincident(first_coords, second_coords))
        {
          coincident_pair.first  = *iter;
          coincident_pair.second = *jter;
          coin_nodes.push_back(coincident_pair);
        }
      }
    }
  }
  return MB_SUCCESS;
}

ErrorCode find_coincident_edges(Range entities,
  std::vector< std::pair<EntityHandle,EntityHandle> > &coin_edges)
{
  double coords1[3], coords2[3], coords3[3];
  Range::iterator iter, jter;
  std::vector<EntityHandle> conn(2);
  std::pair<EntityHandle, EntityHandle> coincident_pair;

  for(iter = entities.begin(); iter != entities.end(); ++iter)
  {
    if(gMB->get_connectivity(&*iter, 1, conn) != MB_SUCCESS)
      return MB_FAILURE;

    // Get the coordinates for the edge endpoints.
    if(gMB->get_coords(&conn[0], 1, coords1) != MB_SUCCESS)
      return MB_FAILURE;

    if(gMB->get_coords(&conn[1], 1, coords2) != MB_SUCCESS)
      return MB_FAILURE;

    for(jter = iter; jter != entities.end(); ++jter)
    {
      if(*iter != *jter)
      {
        // Edges should be the same sense to merge.
        if(gMB->get_connectivity(&*jter, 1, conn) != MB_SUCCESS)
          return MB_FAILURE;

        if(gMB->get_coords(&conn[0], 1, coords3) != MB_SUCCESS)
          return MB_FAILURE;

        if(points_are_coincident(coords1, coords3))
        {
          if(gMB->get_coords(&conn[1], 1, coords3) != MB_SUCCESS)
            return MB_FAILURE;

          if(points_are_coincident(coords2, coords3))
          {
            coincident_pair.first  = *iter;
            coincident_pair.second = *jter;
            coin_edges.push_back(coincident_pair);
          }
        }
      }
    }
  }

  return MB_SUCCESS;
}

ErrorCode find_coincident_elements(Range entities, int num_nodes,
  std::vector< std::pair<EntityHandle,EntityHandle> > &coin)
{
  double coords1[8][3], coords2[8][3];
  Range::iterator iter, jter;
  std::vector<EntityHandle> conn(8);
  std::pair<EntityHandle, EntityHandle> coincident_pair;
  int i = 0, j = 0, ii = 0;

  for(iter = entities.begin(); iter != entities.end(); ++iter)
  {
    // Get the coordinates for the element corners.
    if(gMB->get_connectivity(&*iter, 1, conn) != MB_SUCCESS)
      return MB_FAILURE;
    for(ii=0; ii<num_nodes; ii++)
    {
      if(gMB->get_coords(&conn[ii], 1, coords1[ii]) != MB_SUCCESS)
        return MB_FAILURE;
    }

    for(jter = iter; jter != entities.end(); ++jter)
    {
      if(*iter != *jter)
      {
        // Elements should be the same sense to merge.
        if(gMB->get_connectivity(&*jter, 1, conn) != MB_SUCCESS)
          return MB_FAILURE;

        if(gMB->get_coords(&conn[0], 1, coords2[0]) != MB_SUCCESS)
          return MB_FAILURE;

        // Find if first node is coincident before testing the rest.
        for(i=0; i<num_nodes; i++)
        {
          if(points_are_coincident(coords1[i], coords2[0]))
            break;
        }

        if(i < num_nodes)
        {
          for(ii=1; ii<num_nodes; ii++)
          {
            if(gMB->get_coords(&conn[ii], 1, coords2[ii]) != MB_SUCCESS)
              return MB_FAILURE;
          }

          for(j=1; j<num_nodes; j++)
          {
            if(!points_are_coincident(coords1[j], coords2[(j+i)%num_nodes]))
              break;
          }

          if(j == num_nodes)

          {
            coincident_pair.first  = *iter;
            coincident_pair.second = *jter;
            coin.push_back(coincident_pair);
          }
        }
      }
    }
  }

  return MB_SUCCESS;
}

ErrorCode coincident_counts(EntityCount &curr_count,
                               EntityCount &diff_count)
{
  Range entities;
  std::vector< std::pair<EntityHandle,EntityHandle> > coincident;


  if(curr_count.node > 0)
  {
    if(gMB->get_entities_by_type(0, MBVERTEX, entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_nodes(entities, coincident);
    diff_count.node = coincident.size();
    entities.clear();
    coincident.clear();
  }
  if(curr_count.edge > 0)
  {
    if(gMB->get_entities_by_type(0, MBEDGE, entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_edges(entities, coincident);
    diff_count.edge = coincident.size();
    entities.clear();
    coincident.clear();
  }
  if(curr_count.quad > 0)
  {
    if(gMB->get_entities_by_type(0, MBQUAD, entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 4, coincident);
    diff_count.quad = coincident.size();
    entities.clear();
    coincident.clear();
  }
  if(curr_count.tri > 0)
  {
    if(gMB->get_entities_by_type(0, MBTRI, entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 3, coincident);
    diff_count.tri = coincident.size();
    entities.clear();
    coincident.clear();
  }
  if(curr_count.tet > 0)
  {
    if(gMB->get_entities_by_type(0, MBTET, entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 4, coincident);
    diff_count.tet = coincident.size();
    entities.clear();
    coincident.clear();
  }
  if(curr_count.hex > 0)
  {
    if(gMB->get_entities_by_type(0, MBHEX, entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 8, coincident);
    diff_count.hex = coincident.size();
    entities.clear();
    coincident.clear();
  }

  return MB_SUCCESS;
}

ErrorCode merge_top_down(EntityCount &init_count,
                            EntityCount &curr_count)
{
  Range entities;
  EntityCount diff_count;
  std::vector< std::pair<EntityHandle,EntityHandle> > coincident;

  // Find how many objects of each type need to be merged.
  if(coincident_counts(curr_count, diff_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Find the top level object to merge.
  if(diff_count.hex > 0)
  {
    if(gMB->get_entities_by_type(0, MBHEX,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 8, coincident);
  }
  else if(diff_count.quad > 0)
  {
    if(gMB->get_entities_by_type(0, MBQUAD,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 4, coincident);
  }
  else if(diff_count.tet > 0)
  {
    if(gMB->get_entities_by_type(0, MBTET,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 4, coincident);
  }
  else if(diff_count.tri > 0)

  {
    if(gMB->get_entities_by_type(0, MBTRI,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 3, coincident);
  }
  else if(diff_count.edge > 0)
  {
    if(gMB->get_entities_by_type(0, MBEDGE,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_edges(entities, coincident);
  }

  std::vector< std::pair<EntityHandle,EntityHandle> >::iterator iter;
  for(iter=coincident.begin(); iter != coincident.end(); ++iter)
    gMB->merge_entities((*iter).first, (*iter).second, false, true);

  // Get the new entity totals.
  curr_count.get_counts();

  // Make sure we didn't merge anything.
  if(init_count.node != curr_count.node ||
     init_count.edge != curr_count.edge ||
     init_count.quad != curr_count.quad ||
     init_count.tri != curr_count.tri ||
     init_count.hex != curr_count.hex ||
     init_count.tet != curr_count.tet)
  {
    cout << "***ERROR: Merged top down when not using auto merge.***" << endl;
    return MB_FAILURE;

  }

  return MB_SUCCESS;
}

ErrorCode merge_nodes(EntityCount &init_count, EntityCount &curr_count)
{
  cout << "Merging Coincident Nodes:" << endl;

  // Get the list of vertices from the database.
  Range vertices;
  ErrorCode result = gMB->get_entities_by_type(0, MBVERTEX,  vertices);
  if(result != MB_SUCCESS)
    return result;

  // find the coincident node pairs
  std::vector< std::pair<EntityHandle,EntityHandle> > coincident_nodes;
  find_coincident_nodes(vertices, coincident_nodes);

  // merge the coincident nodes
  std::vector< std::pair<EntityHandle,EntityHandle> >::iterator iter;
  for(iter=coincident_nodes.begin(); iter != coincident_nodes.end(); ++iter)
  {
    cout << "   Coincident nodes: " << (*iter).first << "-"
         << (*iter).second << endl;
    result = gMB->merge_entities((*iter).first, (*iter).second, false, true);
    if(result != MB_SUCCESS)
      return result;
  }

  // Get the reduced list of vertices.
  curr_count.get_counts();



  // Make sure the coincident nodes were all merged.
  if(init_count.node - curr_count.node != coincident_nodes.size())
  {
    cout << "***ERROR: Not all coincident nodes were merged.***" << endl;
    return MB_FAILURE;
  }

  init_count.node = curr_count.node;

  // Make sure we didn't merge anything else.
  if(init_count.edge != curr_count.edge ||
     init_count.quad != curr_count.quad ||
     init_count.tri != curr_count.tri ||
     init_count.hex != curr_count.hex ||
     init_count.tet != curr_count.tet)
  {
    cout << "***ERROR: Merged other objects when merging nodes.***" << endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode merge_edges(EntityCount &init_count, EntityCount &curr_count)
{
  cout << "Merging Coincident Edges:" << endl;

  // Get the list of entities from the database.
  Range entities;
  ErrorCode result = gMB->get_entities_by_type(0, MBEDGE,  entities);
  if(result != MB_SUCCESS)
    return result;

  // find the coincident edge pairs
  std::vector< std::pair<EntityHandle,EntityHandle> > coincident_edges;
  find_coincident_edges(entities, coincident_edges);

  // merge the coincident edges
  unsigned long id1, id2;
  std::vector< std::pair<EntityHandle,EntityHandle> >::iterator iter;
  for(iter=coincident_edges.begin(); iter != coincident_edges.end(); ++iter)
  {
    id1 = gMB->id_from_handle((*iter).first);
    id2 = gMB->id_from_handle((*iter).second);
    cout << "   Coincident edges: " << id1 << "-" << id2 << endl;
    result = gMB->merge_entities((*iter).first, (*iter).second, false, true);
    if(result != MB_SUCCESS)
      return result;
  }

  // Get the reduced list of edges.
  curr_count.get_counts();

  // Make sure the coincident edges were all merged.
  if(init_count.edge - curr_count.edge != coincident_edges.size())
  {
    cout << "***ERROR: Not all coincident edges were merged.***" << endl;
    return MB_FAILURE;
  }

  init_count.edge = curr_count.edge;

  // Make sure we didn't merge anything else.
  if(init_count.node != curr_count.node ||
     init_count.quad != curr_count.quad ||
     init_count.tri != curr_count.tri ||
     init_count.hex != curr_count.hex ||
     init_count.tet != curr_count.tet)
  {
    cout << "***ERROR: Merged other objects when merging edges.***" << endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode merge_2D_elem(EntityCount &init_count, EntityCount &curr_count)
{
  cout << "Merging Coincident 2D Elements:" << endl;

  // Get the list of tris from the database.
  Range entities;
  ErrorCode result = gMB->get_entities_by_type(0, MBTRI,  entities);
  if(result != MB_SUCCESS)
    return result;

  // find the coincident tri pairs
  std::vector< std::pair<EntityHandle,EntityHandle> > coincident;
  find_coincident_elements(entities, 3, coincident);

  // merge the coincident tris
  unsigned long id1, id2;
  unsigned int tri_diff = coincident.size();
  std::vector< std::pair<EntityHandle,EntityHandle> >::iterator iter;
  for(iter=coincident.begin(); iter != coincident.end(); ++iter)
  {
    id1 = gMB->id_from_handle((*iter).first);
    id2 = gMB->id_from_handle((*iter).second);
    cout << "   Coincident tris: " << id1 << "-" << id2 << endl;
    result = gMB->merge_entities((*iter).first, (*iter).second, false, true);

    if(result != MB_SUCCESS)
      return result;
  }

  // Get the list of quads from the database.
  entities.clear();
  result = gMB->get_entities_by_type(0, MBQUAD,  entities);
  if(result != MB_SUCCESS)
    return result;

  // find the coincident tri pairs
  coincident.clear();
  find_coincident_elements(entities, 4, coincident);

  // merge the coincident tris
  unsigned int quad_diff = coincident.size();
  for(iter=coincident.begin(); iter != coincident.end(); ++iter)
  {
    id1 = gMB->id_from_handle((*iter).first);
    id2 = gMB->id_from_handle((*iter).second);
    cout << "   Coincident quads: " << id1 << "-" << id2 << endl;
    result = gMB->merge_entities((*iter).first, (*iter).second, false, true);
    if(result != MB_SUCCESS)
      return result;
  }

  // Get the reduced list of faces.
  curr_count.get_counts();

  // Make sure the coincident faces were all merged.
  if(init_count.tri - curr_count.tri != tri_diff)
  {
    cout << "***ERROR: Not all coincident tris were merged.***" << endl;
    return MB_FAILURE;
  }
  if(init_count.quad - curr_count.quad != quad_diff)
  {
    cout << "***ERROR: Not all coincident quads were merged.***" << endl;
    return MB_FAILURE;
  }

  init_count.tri = curr_count.tri;
  init_count.quad = curr_count.quad;

  // Make sure we didn't merge anything else.
  if(init_count.node != curr_count.node ||
     init_count.edge != curr_count.edge ||
     init_count.hex != curr_count.hex ||
     init_count.tet != curr_count.tet)
  {
    cout << "***ERROR: Merged other objects when merging faces.***" << endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode merge_3D_elem(EntityCount &init_count, EntityCount &curr_count)
{
  cout << "Merging Coincident 3D Elements:" << endl;


  // Get the list of tets from the database.
  Range entities;
  ErrorCode result = gMB->get_entities_by_type(0, MBTET,  entities);
  if(result != MB_SUCCESS)

    return result;

  // find the coincident tet pairs
  std::vector< std::pair<EntityHandle,EntityHandle> > coincident;
  find_coincident_elements(entities, 4, coincident);

  // merge the coincident tets
  unsigned long id1, id2;
  unsigned int tet_diff = coincident.size();
  std::vector< std::pair<EntityHandle,EntityHandle> >::iterator iter;
  for(iter=coincident.begin(); iter != coincident.end(); ++iter)
  {
    id1 = gMB->id_from_handle((*iter).first);
    id2 = gMB->id_from_handle((*iter).second);
    cout << "   Coincident tets: " << id1 << "-" << id2 << endl;
    result = gMB->merge_entities((*iter).first, (*iter).second, false, true);
    if(result != MB_SUCCESS)
      return result;
  }


  // Get the list of hexs from the database.
  entities.clear();
  result = gMB->get_entities_by_type(0, MBHEX,  entities);
  if(result != MB_SUCCESS)
    return result;

  // find the coincident hex pairs
  coincident.clear();
  find_coincident_elements(entities, 8, coincident);

  // merge the coincident tris
  unsigned int hex_diff = coincident.size();
  for(iter=coincident.begin(); iter != coincident.end(); ++iter)
  {
    id1 = gMB->id_from_handle((*iter).first);
    id2 = gMB->id_from_handle((*iter).second);
    cout << "   Coincident hexs: " << id1 << "-" << id2 << endl;
    result = gMB->merge_entities((*iter).first, (*iter).second, false, true);
    if(result != MB_SUCCESS)
      return result;
  }

  // Get the reduced list of elements.
  curr_count.get_counts();

  // Make sure the coincident elements were all merged.
  if(init_count.tet - curr_count.tet != tet_diff)
  {
    cout << "***ERROR: Not all coincident tets were merged.***" << endl;
    return MB_FAILURE;
  }
  if(init_count.hex - curr_count.hex != hex_diff)
  {
    cout << "***ERROR: Not all coincident hexs were merged.***" << endl;
    return MB_FAILURE;
  }

  init_count.tet = curr_count.tet;
  init_count.hex = curr_count.hex;

  // Make sure we didn't merge anything else.
  if(init_count.node != curr_count.node ||
     init_count.edge != curr_count.edge ||
     init_count.quad != curr_count.quad ||
     init_count.tri != curr_count.tri)
  {
    cout << "***ERROR: Merged other objects when merging elements.***" << endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode read_file(std::string &file_name, EntityCount &counts)
{
  // Make sure the database is empty.
  gMB->delete_mesh();

  // Read the model from the file.
  if(gMB->load_mesh(file_name.c_str(), 0) != MB_SUCCESS)
  {
    cout << "***ERROR: Unable to load mesh file.***" << endl;
    return MB_FAILURE;
  }

  // Get the number of each entity types in the mesh.
  if(counts.get_counts() != MB_SUCCESS)
  {
    cout << "***ERROR: Unable to get entity list counts.***" << endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode write_file(std::string &file_name)
{
  //get the block tag
  ErrorCode result;
  Range block_range;
  Tag block_tag;

  if(gMB->tag_get_handle("MATERIAL_SET", block_tag) == MB_SUCCESS)
  {
    //get all the blocks
    result = gMB->get_entities_by_type_and_tag(0, MBENTITYSET, &block_tag, 0, 1, block_range);
    if(result != MB_SUCCESS)
      return result;
  }

  //transfer range contents into vectors
  std::vector<EntityHandle> output_list;

  Range::iterator range_iter, end_iter;
  range_iter = block_range.begin();
  end_iter = block_range.end();

  for(; range_iter != end_iter; ++range_iter)
  {
    int id;
    result = gMB->tag_get_handle("MATERIAL_SET",  block_tag);
    if(result != MB_SUCCESS)
      return result;

    result = gMB->tag_get_data(block_tag, &*range_iter, 1, &id);
    if(result != MB_SUCCESS)
      return result;

    //if(id != 2)
    output_list.push_back( *range_iter );
  }

  // write the file
  static std::string mrg_splice(".mrg");
  file_name.insert(file_name.size()-2, mrg_splice);
  result = gMB->write_mesh(file_name.c_str(), &output_list[0], output_list.size());

  return result;
}

ErrorCode process_td_auto_merge(std::string &file_name)
{
  EntityCount init_count;
  EntityCount curr_count;
  EntityCount diff_count;

  // Read in the mesh and get the number of entities of each type.
  if(read_file(file_name, init_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Copy the initial counts into the current count object.
  curr_count.copy_counts(init_count);

  // Print out the list of initial objects.
  cout << "Initial Entities:" << endl;
  curr_count.print();

  // Try auto merging from the top down.
  Range entities;
  std::vector< std::pair<EntityHandle,EntityHandle> > coincident;

  // Find how many objects of each type need to be merged.
  if(coincident_counts(curr_count, diff_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Find the top level object to merge.
  if(diff_count.hex > 0)
  {
    if(gMB->get_entities_by_type(0, MBHEX,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 8, coincident);
  }
  else if(diff_count.quad > 0)
  {
    if(gMB->get_entities_by_type(0, MBQUAD,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 4, coincident);
  }
  else if(diff_count.tet > 0)
  {
    if(gMB->get_entities_by_type(0, MBTET,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 4, coincident);
  }
  else if(diff_count.tri > 0)
  {
    if(gMB->get_entities_by_type(0, MBTRI,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_elements(entities, 3, coincident);
  }
  else if(diff_count.edge > 0)
  {
    if(gMB->get_entities_by_type(0, MBEDGE,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_edges(entities, coincident);
  }
  else if(diff_count.node > 0)
  {
    if(gMB->get_entities_by_type(0, MBVERTEX,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_nodes(entities, coincident);
  }

  cout << "Merging coincident entities(top down)..." << endl;
  std::vector< std::pair<EntityHandle,EntityHandle> >::iterator iter;
  for(iter=coincident.begin(); iter != coincident.end(); ++iter)
    gMB->merge_entities((*iter).first, (*iter).second, true, true);

  // Get the new entity totals.
  if(curr_count.get_counts() != MB_SUCCESS)
  {
    cout << "***ERROR: Unable to get entity list counts.***" << endl;
    return MB_FAILURE;
  }

  // Make sure we merged everything.
  if(init_count.node - curr_count.node != diff_count.node ||
     init_count.edge - curr_count.edge != diff_count.edge ||
     init_count.quad - curr_count.quad != diff_count.quad ||
     init_count.tri - curr_count.tri != diff_count.tri ||
     init_count.hex - curr_count.hex != diff_count.hex ||
     init_count.tet - curr_count.tet != diff_count.tet)
  {
    cout << "***ERROR: Not all coincident objects merged automatically.***"
         << endl;
    curr_count.print();
    return MB_FAILURE;
  }

  // Print out the final list of objects.
  cout << "Final Entities:" << endl;
  curr_count.print();

  return MB_SUCCESS;
}

ErrorCode process_mo_auto_merge(std::string &file_name)
{
  EntityCount init_count;
  EntityCount curr_count;
  EntityCount diff_count;

  // Read in the mesh and get the number of entities of each type.
  if(read_file(file_name, init_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Copy the initial counts into the current count object.
  curr_count.copy_counts(init_count);

  // Print out the list of initial objects.
  cout << "Initial Entities:" << endl;
  curr_count.print();

  // Try auto merging from the middle out.
  Range entities;
  std::vector< std::pair<EntityHandle,EntityHandle> > coincident;

  // Find how many objects of each type need to be merged.
  if(coincident_counts(curr_count, diff_count) != MB_SUCCESS)
    return MB_FAILURE;



  // Start merge at the edge level if present.
  if(diff_count.edge > 0)
  {
    if(gMB->get_entities_by_type(0, MBEDGE,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_edges(entities, coincident);
  }
  else
  {
    if(gMB->get_entities_by_type(0, MBVERTEX,  entities) != MB_SUCCESS)
      return MB_FAILURE;
    find_coincident_nodes(entities, coincident);
  }

  cout << "Merging coincident entities(middle out)..." << endl;
  std::vector< std::pair<EntityHandle,EntityHandle> >::iterator iter;
  for(iter=coincident.begin(); iter != coincident.end(); ++iter)

    gMB->merge_entities((*iter).first, (*iter).second, true,true);

  // Get the new entity totals.
  if(curr_count.get_counts() != MB_SUCCESS)
  {
    cout << "***ERROR: Unable to get entity list counts.***" << endl;
    return MB_FAILURE;
  }

  // Make sure we merged everything.
  if(init_count.node - curr_count.node != diff_count.node ||
     init_count.edge - curr_count.edge != diff_count.edge ||
     init_count.quad - curr_count.quad != diff_count.quad ||
     init_count.tri - curr_count.tri != diff_count.tri ||
     init_count.hex - curr_count.hex != diff_count.hex ||
     init_count.tet - curr_count.tet != diff_count.tet)
  {
    cout << "***ERROR: Not all coincident objects merged automatically.***"
         << endl;
    curr_count.print();
    return MB_FAILURE;
  }

  // Print out the final list of objects.
  cout << "Final Entities:" << endl;
  curr_count.print();

  return MB_SUCCESS;
}

ErrorCode process_bu_auto_merge(std::string &file_name)
{
  EntityCount init_count;
  EntityCount curr_count;
  EntityCount diff_count;

  // Read in the mesh and get the number of entities of each type.
  if(read_file(file_name, init_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Copy the initial counts into the current count object.
  curr_count.copy_counts(init_count);

  // Print out the list of initial objects.
  cout << "Initial Entities:" << endl;
  curr_count.print();

  // Try auto merging from the bottom up.
  Range entities;
  std::vector< std::pair<EntityHandle,EntityHandle> > coincident;

  // Find how many objects of each type need to be merged.
  if(coincident_counts(curr_count, diff_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Start merging from the nodes and go up.
  if(gMB->get_entities_by_type(0, MBVERTEX,  entities) != MB_SUCCESS)
    return MB_FAILURE;
  find_coincident_nodes(entities, coincident);

  cout << "Merging coincident entities(bottom up)..." << endl;
  std::vector< std::pair<EntityHandle,EntityHandle> >::iterator iter;
  for(iter=coincident.begin(); iter != coincident.end(); ++iter)
    gMB->merge_entities((*iter).first, (*iter).second, true, true);

  // Get the new entity totals.
  if(curr_count.get_counts() != MB_SUCCESS)
  {
    cout << "***ERROR: Unable to get entity list counts.***" << endl;
    return MB_FAILURE;
  }

  // Make sure we merged everything.
  if(init_count.node - curr_count.node != diff_count.node ||
     init_count.edge - curr_count.edge != diff_count.edge ||
     init_count.quad - curr_count.quad != diff_count.quad ||
     init_count.tri - curr_count.tri != diff_count.tri ||
     init_count.hex - curr_count.hex != diff_count.hex ||
     init_count.tet - curr_count.tet != diff_count.tet)
  {
    cout << "***ERROR: Not all coincident objects merged automatically.***"
         << endl;
    curr_count.print();

    return MB_FAILURE;
  }

  // Print out the final list of objects.
  cout << "Final Entities:" << endl;
  curr_count.print();

  return MB_SUCCESS;
}

ErrorCode process_merge(std::string &file_name)
{
  EntityCount init_count;
  EntityCount curr_count;

  // Read in the mesh and get the number of entities of each type.
  if(read_file(file_name, init_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Copy the initial counts into the current count object.
  curr_count.copy_counts(init_count);

  // Print out the list of initial objects.
  cout << "Initial Entities:" << endl;
  curr_count.print();

  // Try to merge elements before nodes (should fail).
  if(merge_top_down(init_count, curr_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Merge the nodes.
  if(merge_nodes(init_count, curr_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Next, merge the edges.
  if(merge_edges(init_count, curr_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Now, merge the 2D elements.
  if(merge_2D_elem(init_count, curr_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Finally, merge the 3D elements.
  if(merge_3D_elem(init_count, curr_count) != MB_SUCCESS)
    return MB_FAILURE;

  // Print out the final list of objects.
  if(curr_count.get_counts() != MB_SUCCESS)
  {
    cout << "***ERROR: Unable to get entity list counts.***" << endl;
    return MB_FAILURE;
  }
  cout << "Final Entities:" << endl;
  curr_count.print();

  // write the file out (modifies name)
  return write_file(file_name);
}


int main()
{
  ErrorCode result;
  std::string test_files[] = {std::string("test/2barcase1.g"),
                              std::string("test/2barcase2.g"),
                              std::string("test/2hexcase1.g"),
                              std::string("test/2hexcase2.g"),
                              std::string("test/2hexcase3.g"),
                              std::string("test/2hexcase4.g"),
                              std::string("test/2hexcase5.g"),
                              std::string("test/2quadcase1.g"),
                              std::string("test/2quadcase2.g"),
                              std::string("test/2quadcase3.g"),
                              std::string("test/2quadcase4.g"),
                              std::string("test/2tetcase1.g"),
                              std::string("test/2tetcase2.g"),
                              std::string("test/2tetcase3.g"),
                              std::string("test/2tetcase4.g"),
                              std::string("test/2tricase1.g"),
                              std::string("test/2tricase2.g"),
                              std::string("test/2tricase3.g")};

  // Create the MB database instance.
  gMB = new Core();

  // Loop through the list of test files.
  unsigned int i;
  cout << "---Starting Top Down Auto Merge Tests---" << endl << endl;
  for(i=0; i<(sizeof(test_files)/sizeof(std::string)); i++)
  {
    cout << "---Testing:\"" << test_files[i] << "\"---" << endl;
    result = process_td_auto_merge(test_files[i]);
    if(result == MB_SUCCESS)
      cout << "---Success---";
    else
      cout << "---Failure---";
    cout << endl << endl;
  }

  cout << "---Starting Bottom Up Auto Merge Tests---" << endl << endl;
  for(i=0; i<(sizeof(test_files)/sizeof(std::string)); i++)
  {
    cout << "---Testing:\"" << test_files[i] << "\"---" << endl;
    result = process_bu_auto_merge(test_files[i]);
    if(result == MB_SUCCESS)

      cout << "---Success---";
    else
      cout << "---Failure---";
    cout << endl << endl;
  }

  cout << "---Starting Middle Out Auto Merge Tests---" << endl << endl;
  for(i=0; i<(sizeof(test_files)/sizeof(std::string)); i++)
  {
    cout << "---Testing:\"" << test_files[i] << "\"---" << endl;
    result = process_mo_auto_merge(test_files[i]);
    if(result == MB_SUCCESS)
      cout << "---Success---";
    else
      cout << "---Failure---";

    cout << endl << endl;
  }

  cout << "---Starting Merge Tests---" << endl << endl;
  for(i=0; i<(sizeof(test_files)/sizeof(std::string)); i++)
  {
    cout << "---Testing:\"" << test_files[i] << "\"---" << endl;
    result = process_merge(test_files[i]);
    if(result == MB_SUCCESS)
      cout << "---Success---";
    else
      cout << "---Failure---";
    cout << endl << endl;
  }

  // Destroy the MB database instance.
  delete gMB;
  gMB = NULL;
}
