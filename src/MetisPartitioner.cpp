// SpaFEDTe, a Template based C++ library for creating 
// Discontinuous Finite Element Spaces,
// Copyright (C) 2012 Lorenzo Alessio Botti

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */ 
/* License as published by the Free Software Foundation either */ 
/* version 3.0 of the License, or (at your option) any later version. */

/* This software is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this software; if not, a copy of the full */
/* GNU Lesser General Public License can be found at */
/* http://www.gnu.org/licenses/ */

// This implementation is mostly borrowed from the mbzoltan MOAB partitioning tool

#include <iostream>
#include <assert.h>
#include <sstream>

#include "MetisMOABPartitioner.hpp"

using namespace moab;

const bool debug = false;

MetisMOABPartitioner::MetisMOABPartitioner( Interface *impl, 
                                            const bool use_coords,
                                            int argc, 
                                            char **argv) 
                                           : mbImpl(impl), 
                                             useCoords(use_coords),
                                             argcArg(argc), 
                                             argvArg(argv)
{
  mbpc = ParallelComm::get_pcomm(mbImpl, 0);
  if (!mbpc)
    mbpc = new ParallelComm( impl, MPI_COMM_WORLD, 0 );
}

MetisMOABPartitioner::~MetisMOABPartitioner() 
{
  ;
}

ErrorCode MetisMOABPartitioner::partition_mesh_geom(const int nparts,
                                                    const char *method,
                                                    const int part_dim,
                                                    const bool write_as_sets,
                                                    const bool write_as_tags,
					            const bool partition_tagged_sets,
					            const bool partition_tagged_ents,
					            const char *aggregating_tag)
{
    // should only be called in serial
  if (mbpc->proc_config().proc_size() != 1) {
    std::cout << "MetisMOABPartitioner::partition_mesh_geom must be called in serial." 
              << std::endl;
    return MB_FAILURE;
  }
  
  if (NULL != method && strcmp(method, "ML_RB") && strcmp(method, "ML_KWAY"))
  {
    std::cout << "ERROR node " << mbpc->proc_config().proc_rank() << ": Method must be "
              << "ML_RB or ML_KWAY"
              << std::endl;
    return MB_FAILURE;
  }
  
  std::vector<double> pts; // x[0], y[0], z[0], ... from MOAB
  std::vector<int> ids; // point ids from MOAB
  std::vector<int> adjs, length, parts;
  Range elems;
  // Get a mesh from MOAB and diide it across processors.

  ErrorCode result;
  std::cout << "Assembling graph..." << std::endl;
  if (!partition_tagged_sets && !partition_tagged_ents)
  {
    result = assemble_graph(part_dim, pts, ids, adjs, length, elems);
  }
  else if (partition_tagged_sets) 
  {
    result = assemble_taggedsets_graph(part_dim, pts, ids, adjs, length, elems, &(*aggregating_tag)); 
  }
  else if (partition_tagged_ents) 
  {
    result = assemble_taggedents_graph(part_dim, pts, ids, adjs, length, elems, &(*aggregating_tag)); 
  }
  if (MB_SUCCESS != result) return result;
  
  std::cout << "Computing partition using " << method 
            <<" method for " << nparts << " processors..." << std::endl;

  int nelems = length.size()-1;
  int *assign_parts;
  assign_parts = (int *)malloc(sizeof(int) * nelems);
  int nconstraints = 1;
  int edgeCut = 0;
  int nOfPartitions = nparts;
  int metis_RESULT;

  if (strcmp(method, "ML_KWAY") == 0)
  {
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_CONTIG] = 1;  
    metis_RESULT = METIS_PartGraphKway(&nelems, &nconstraints, &length[0], &adjs[0], NULL, NULL, NULL, &nOfPartitions, NULL, NULL, options, &edgeCut, assign_parts);
  }
  else if (strcmp(method, "ML_RB") == 0)
  {
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // CUT 
    options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW; // GROW or RANDOM
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM; // RM or SHEM
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // FM
    options[METIS_OPTION_NCUTS] = 10; // Number of different partitionings to compute, then chooses the best one, default = 1
    options[METIS_OPTION_NITER] = 10;  // Number of refinements steps, default = 10
    options[METIS_OPTION_UFACTOR] = 30; // Imabalance, default = 1
    options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;
    metis_RESULT = METIS_PartGraphRecursive(&nelems, &nconstraints, &length[0], &adjs[0], NULL, NULL, NULL, &nOfPartitions, NULL, NULL, options, &edgeCut, assign_parts);
  }
    // assign global node ids, starting from one! TODO
  result = mbpc->assign_global_ids(0, 0, 1); 

  if (metis_RESULT != METIS_OK)
    return MB_FAILURE;
  
  // take results & write onto MOAB partition sets
  std::cout << "Saving partition information to MOAB..." << std::endl;

  if (partition_tagged_sets || partition_tagged_ents)
    result = write_aggregationtag_partition(nparts, elems, assign_parts,
                                            write_as_sets, write_as_tags);
  else
    result = write_partition(nparts, elems, assign_parts,
                             write_as_sets, write_as_tags);

  free(assign_parts);
  if (MB_SUCCESS != result) return result;

  return MB_SUCCESS;
}

ErrorCode MetisMOABPartitioner::assemble_taggedents_graph(const int dimension,
                                                          std::vector<double> &coords,
                                                          std::vector<int> &moab_ids,
                                                          std::vector<int> &adjacencies, 
                                                          std::vector<int> &length,
                                                          Range &elems,
					                  const char *aggregating_tag)
{
  Tag partSetTag;
  ErrorCode result = mbImpl->tag_get_handle(aggregating_tag, 1, MB_TYPE_INTEGER, partSetTag);
  if (MB_SUCCESS != result) return result;
 
  Range allSubElems;
  result = mbImpl->get_entities_by_dimension(0, dimension, allSubElems);
  if (MB_SUCCESS != result || allSubElems.empty()) return result;
  int partSet;
  std::map<int, Range> aggloElems;
  for (Range::iterator rit = allSubElems.begin(); rit != allSubElems.end(); rit++) 
  {
    EntityHandle entity = *rit;
    result = mbImpl->tag_get_data(partSetTag,&entity,1,&partSet);
    if (MB_SUCCESS != result) return result;
    if (partSet >= 0)
      aggloElems[partSet].insert(entity);
  }
  // clear aggregating tag data
  TagType type;
  result = mbImpl->tag_get_type(partSetTag, type);
  if (type == MB_TAG_DENSE)
  {
    // clear tag on ents and sets
    result = mbImpl->tag_delete(partSetTag); 
    if (MB_SUCCESS != result) return result;
  }
  if (type == MB_TAG_SPARSE)
  {
    // clear tag on ents
    result = mbImpl->tag_delete_data(partSetTag, allSubElems); 
    if (MB_SUCCESS != result) return result;
    // clear tag on sets
    result = mbImpl->get_entities_by_type_and_tag(0 , MBENTITYSET, &partSetTag, 0, 1, elems);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->tag_delete_data(partSetTag, elems); 
    if (MB_SUCCESS != result) return result;
    elems.clear();
  }
  result = mbImpl->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
                                  partSetTag, MB_TAG_SPARSE|MB_TAG_CREAT); 
  if (MB_SUCCESS != result) return result;
  
  for (std::map<int, Range>::iterator mit = aggloElems.begin(); mit != aggloElems.end(); mit++) 
  {
    EntityHandle new_set;
    result = mbImpl->create_meshset(MESHSET_SET, new_set);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->add_entities(new_set, mit->second);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->tag_set_data (partSetTag, &new_set, 1, &mit->first);
    if (MB_SUCCESS != result) return result;
  }

  result = assemble_taggedsets_graph(dimension, coords, moab_ids, adjacencies, length, elems, &(*aggregating_tag));
  return MB_SUCCESS;
}

ErrorCode MetisMOABPartitioner::assemble_taggedsets_graph(const int dimension,
                                                          std::vector<double> &coords,
                                                          std::vector<int> &moab_ids,
                                                          std::vector<int> &adjacencies, 
                                                          std::vector<int> &length,
                                                          Range &elems,
					                  const char *aggregating_tag)
{
  length.push_back(0);
    // assemble a graph with vertices equal to elements of specified dimension, edges
    // signified by list of other elements to which an element is connected

  // get the tagged elements 
  Tag partSetTag;
  ErrorCode result = mbImpl->tag_get_handle(aggregating_tag, 1, MB_TYPE_INTEGER, partSetTag);
  //ErrorCode result = mbImpl->tag_get_handle("PARALLEL_PARTITION_SET", 1, MB_TYPE_INTEGER, partSetTag);
  if (MB_SUCCESS != result) return result;

  result = mbImpl->get_entities_by_type_and_tag(0 , MBENTITYSET, &partSetTag, 0, 1, elems);
  if (MB_SUCCESS != result || elems.empty()) return result;

  //assign globla ids to elem sets based on aggregating_tag data 
  Tag gid_tag;
  int zero1 = -1;
  result = mbImpl->tag_get_handle("GLOBAL_ID_AGGLO", 1, MB_TYPE_INTEGER, gid_tag, MB_TAG_SPARSE|MB_TAG_CREAT, &zero1);
  if (MB_SUCCESS != result) return result;
  for (Range::iterator rit = elems.begin(); rit != elems.end(); rit++) 
  {
    int partSet;
    result = mbImpl->tag_get_data(partSetTag,&(*rit),1,&partSet);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->tag_set_data(gid_tag, &(*rit), 1, &partSet);
    if (MB_SUCCESS != result) return result;
  }
  // clear aggregating tag data
  TagType type;
  result = mbImpl->tag_get_type(partSetTag, type);
  if (type == MB_TAG_DENSE)
  {
    result = mbImpl->tag_delete(partSetTag); 
    if (MB_SUCCESS != result) return result;
  }
  if (type == MB_TAG_SPARSE)
  {
    result = mbImpl->tag_delete_data(partSetTag, elems); 
    if (MB_SUCCESS != result) return result;
  }
  
  // assemble the graph, using Skinner to get d-1 dimensional neighbors and then intersecting to get adjacencies
  std::vector<Range> skin_subFaces(elems.size());
  unsigned int i = 0;
  for (Range::iterator rit = elems.begin(); rit != elems.end(); rit++) 
  {
    Range part_ents;
    result = mbImpl->get_entities_by_handle(*rit, part_ents, false);
    if (mbImpl->dimension_from_handle(*part_ents.rbegin()) != mbImpl->dimension_from_handle(*part_ents.begin())) 
    {
      Range::iterator lower = part_ents.lower_bound(CN::TypeDimensionMap[0].first),
      upper = part_ents.upper_bound(CN::TypeDimensionMap[dimension-1].second);
      part_ents.erase(lower, upper);
    }
    Skinner skinner(mbImpl);
    result = skinner.find_skin(part_ents, false, skin_subFaces[i], NULL, false, true, false);
    if (MB_SUCCESS != result) return result;
    i++;
  }
  std::vector<EntityHandle> adjs;
  std::vector<int> neighbors;
  double avg_position[3];
  int moab_id;
  MeshTopoUtil mtu(mbImpl);
  for (unsigned int k = 0; k < i; k++)
  {
      // get bridge adjacencies for element k
    adjs.clear();
    for (unsigned int t = 0; t < i; t++)
    {
      if (t != k)
      {
        Range subFaces = intersect(skin_subFaces[k],skin_subFaces[t]);
        if (subFaces.size() > 0)
  	  adjs.push_back(elems[t]);
      }
    }
    if (!adjs.empty()) 
    {
      neighbors.resize(adjs.size());
      result = mbImpl->tag_get_data(gid_tag, &adjs[0], adjs.size(), &neighbors[0]); 
    }
      // copy those into adjacencies vector
    length.push_back(length.back()+(int)adjs.size());
    std::copy(neighbors.begin(), neighbors.end(), std::back_inserter(adjacencies));
      // get the graph vertex id for this element
    const EntityHandle& setk = elems[k];
    result = mbImpl->tag_get_data(gid_tag, &setk, 1, &moab_id); 
    moab_ids.push_back(moab_id);
      // get average position of vertices
    Range part_ents;
    result = mbImpl->get_entities_by_handle(elems[k], part_ents, false);
    result = mtu.get_average_position(part_ents, avg_position); 
    std::copy(avg_position, avg_position+3, std::back_inserter(coords));
  }
  for (unsigned int k = 0; k < i; k++)
  {
    for (unsigned int t = 0; t < k; t++)
    {
      Range subFaces = intersect(skin_subFaces[k],skin_subFaces[t]);
      if (subFaces.size() > 0)
        mbImpl->delete_entities(subFaces);
    }
  }

  if (debug) {
    std::cout << "Length vector: " << std::endl;
    std::copy(length.begin(), length.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Adjacencies vector: " << std::endl;
    std::copy(adjacencies.begin(), adjacencies.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Moab_ids vector: " << std::endl;
    std::copy(moab_ids.begin(), moab_ids.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Coords vector: " << std::endl;
    std::copy(coords.begin(), coords.end(), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;
  }
  return MB_SUCCESS;
}

ErrorCode MetisMOABPartitioner::assemble_graph(const int dimension,
                                               std::vector<double> &coords,
                                               std::vector<int> &moab_ids,
                                               std::vector<int> &adjacencies, 
                                               std::vector<int> &length,
                                               Range &elems) 
{
  length.push_back(0);
    // assemble a graph with vertices equal to elements of specified dimension, edges
    // signified by list of other elements to which an element is connected

    // get the elements of that dimension
  ErrorCode result = mbImpl->get_entities_by_dimension(0, dimension, elems);
  if (MB_SUCCESS != result || elems.empty()) return result;
  
    // assign global ids
  result = mbpc->assign_global_ids(0, dimension, 0); 

    // now assemble the graph, calling MeshTopoUtil to get bridge adjacencies through d-1 dimensional
    // neighbors
  MeshTopoUtil mtu(mbImpl);
  Range adjs;
    // can use a fixed-size array 'cuz the number of lower-dimensional neighbors is limited
    // by MBCN
  int neighbors[5*MAX_SUB_ENTITIES];
  double avg_position[3];
  int moab_id;
  
    // get the global id tag hanlde
  Tag gid;
  result = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
                                  gid, MB_TAG_DENSE|MB_TAG_CREAT);
  if (MB_SUCCESS != result) return result;
  
  for (Range::iterator rit = elems.begin(); rit != elems.end(); rit++) {

      // get bridge adjacencies
    adjs.clear();
    result = mtu.get_bridge_adjacencies(*rit, (dimension > 0 ? dimension-1 : 3), 
                                        dimension, adjs); 
    
    
      // get the graph vertex ids of those
    if (!adjs.empty()) {
      assert(adjs.size() < 5*MAX_SUB_ENTITIES);
      result = mbImpl->tag_get_data(gid, adjs, neighbors); 
    }

      // copy those into adjacencies vector
    length.push_back(length.back()+(int)adjs.size());
    std::copy(neighbors, neighbors+adjs.size(), std::back_inserter(adjacencies));

      // get average position of vertices
    result = mtu.get_average_position(*rit, avg_position); 
    
      // get the graph vertex id for this element
    result = mbImpl->tag_get_data(gid, &(*rit), 1, &moab_id); 

      // copy those into coords vector
    moab_ids.push_back(moab_id);
    std::copy(avg_position, avg_position+3, std::back_inserter(coords));
  }

  if (debug) {
    std::cout << "Length vector: " << std::endl;
    std::copy(length.begin(), length.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Adjacencies vector: " << std::endl;
    std::copy(adjacencies.begin(), adjacencies.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Moab_ids vector: " << std::endl;
    std::copy(moab_ids.begin(), moab_ids.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Coords vector: " << std::endl;
    std::copy(coords.begin(), coords.end(), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;
  }

  return MB_SUCCESS;
}

ErrorCode MetisMOABPartitioner::write_aggregationtag_partition(const int nparts,
                                                               Range &elems, 
                                                               const int *assignment,
                                                               const bool write_as_sets,
                                                               const bool write_as_tags)
{
  ErrorCode result;

    // get the partition set tag
  Tag part_set_tag;
  int i;
  result = mbImpl->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
                                  part_set_tag, MB_TAG_SPARSE|MB_TAG_CREAT); 
  if (MB_SUCCESS != result) return result;
    // get any sets already with this tag, and clear them
  Range tagged_sets;
  result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &part_set_tag, NULL, 1,
                                                tagged_sets, Interface::UNION); 
  if (!tagged_sets.empty()) {
    result = mbImpl->clear_meshset(tagged_sets); 
    if (!write_as_sets) {
      result = mbImpl->tag_delete_data(part_set_tag, tagged_sets); 
    }
  }
 
  if (write_as_sets) {
      // first, create partition sets and store in vector
    partSets.clear();
  
    if (nparts > (int) tagged_sets.size()) {
        // too few partition sets - create missing ones
      int num_new = nparts - tagged_sets.size();
      for (i = 0; i < num_new; i++) {
        EntityHandle new_set;
        result = mbImpl->create_meshset(MESHSET_SET, new_set); 
        tagged_sets.insert(new_set);
      }
    }
    else if (nparts < (int) tagged_sets.size()) {
        // too many partition sets - delete extras
      int num_del = tagged_sets.size() - nparts;
      for (i = 0; i < num_del; i++) {
        EntityHandle old_set = tagged_sets.pop_back();
        result = mbImpl->delete_entities(&old_set, 1); 
      }
    }
  
      // assign partition sets to vector
    partSets.swap(tagged_sets);
  
      // write a tag to those sets denoting they're partition sets, with a value of the
      // proc number
    int *dum_ids = new int[nparts];
    for (i = 0; i < nparts; i++) dum_ids[i] = i;
  
    result = mbImpl->tag_set_data(part_set_tag, partSets, dum_ids); 

      // assign entities to the relevant sets
    std::vector<EntityHandle> tmp_part_sets;
    std::copy(partSets.begin(), partSets.end(), std::back_inserter(tmp_part_sets));
    Range::iterator rit;
    for (i = 0, rit = elems.begin(); rit != elems.end(); rit++, i++) {
      result = mbImpl->add_entities(tmp_part_sets[assignment[i]], &(*rit), 1); 
    }

      // check for empty sets, warn if there are any
    Range empty_sets;
    for (Range::iterator rit = partSets.begin(); rit != partSets.end(); rit++) {
      int num_ents = 0;
      result = mbImpl->get_number_entities_by_handle(*rit, num_ents);
      if (MB_SUCCESS != result || !num_ents) empty_sets.insert(*rit);
    }
    if (!empty_sets.empty()) {
      std::cout << "WARNING: " << empty_sets.size() << " empty sets in partition: ";
      for (Range::iterator rit = empty_sets.begin(); rit != empty_sets.end(); rit++)
        std::cout << *rit << " ";
      std::cout << std::endl;
    }
  }

  if (write_as_tags) {
    Tag gid_tag;
    result = mbImpl->tag_get_handle("GLOBAL_ID_AGGLO", 1, MB_TYPE_INTEGER, gid_tag, MB_TAG_SPARSE); 
    if (MB_SUCCESS != result) return result;
  
      // allocate integer-size partitions
    unsigned int i = 0;
    int gid;
    for (Range::iterator rit = elems.begin(); rit != elems.end(); rit++) 
    {
      result = mbImpl->tag_get_data(gid_tag, &(*rit), 1, &gid);
      Range part_ents;
//      std::cout<<"part ents "<<part_ents.size()<<std::endl;
      result = mbImpl->get_entities_by_handle(*rit, part_ents, false);
      if (MB_SUCCESS != result) return result;
      for (Range::iterator eit = part_ents.begin(); eit != part_ents.end(); eit++) 
      {
        result = mbImpl->tag_set_data(part_set_tag, &(*eit), 1, &assignment[i]);
        if (MB_SUCCESS != result) return result;
        result = mbImpl->tag_set_data(gid_tag, &(*eit), 1, &gid);
        if (MB_SUCCESS != result) return result;
      }
      i++;
    }
  }
  return MB_SUCCESS;
}

ErrorCode MetisMOABPartitioner::write_partition(const int nparts,
                                                Range &elems, 
                                                const int *assignment,
                                                const bool write_as_sets,
                                                const bool write_as_tags) 
{
  ErrorCode result;

    // get the partition set tag
  Tag part_set_tag;
  int dum_id = -1, i;
  result = mbImpl->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
                                  part_set_tag, MB_TAG_SPARSE|MB_TAG_CREAT, &dum_id); 
  
    // get any sets already with this tag, and clear them
  Range tagged_sets;
  result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &part_set_tag, NULL, 1,
                                                tagged_sets, Interface::UNION); 
  if (!tagged_sets.empty()) {
    result = mbImpl->clear_meshset(tagged_sets); 
    if (!write_as_sets) {
      result = mbImpl->tag_delete_data(part_set_tag, tagged_sets); 
    }
  }

  if (write_as_sets) {
      // first, create partition sets and store in vector
    partSets.clear();
  
    if (nparts > (int) tagged_sets.size()) {
        // too few partition sets - create missing ones
      int num_new = nparts - tagged_sets.size();
      for (i = 0; i < num_new; i++) {
        EntityHandle new_set;
        result = mbImpl->create_meshset(MESHSET_SET, new_set); 
        tagged_sets.insert(new_set);
      }
    }
    else if (nparts < (int) tagged_sets.size()) {
        // too many partition sets - delete extras
      int num_del = tagged_sets.size() - nparts;
      for (i = 0; i < num_del; i++) {
        EntityHandle old_set = tagged_sets.pop_back();
        result = mbImpl->delete_entities(&old_set, 1); 
      }
    }
  
      // assign partition sets to vector
    partSets.swap(tagged_sets);
  
      // write a tag to those sets denoting they're partition sets, with a value of the
      // proc number
    int *dum_ids = new int[nparts];
    for (i = 0; i < nparts; i++) dum_ids[i] = i;
  
    result = mbImpl->tag_set_data(part_set_tag, partSets, dum_ids); 
    delete dum_ids;

      // assign entities to the relevant sets
    std::vector<EntityHandle> tmp_part_sets;
    std::copy(partSets.begin(), partSets.end(), std::back_inserter(tmp_part_sets));
    Range::iterator rit;
    for (i = 0, rit = elems.begin(); rit != elems.end(); rit++, i++) {
      result = mbImpl->add_entities(tmp_part_sets[assignment[i]], &(*rit), 1); 
    }

      // check for empty sets, warn if there are any
    Range empty_sets;
    for (Range::iterator rit = partSets.begin(); rit != partSets.end(); rit++) {
      int num_ents = 0;
      result = mbImpl->get_number_entities_by_handle(*rit, num_ents);
      if (MB_SUCCESS != result || !num_ents) empty_sets.insert(*rit);
    }
    if (!empty_sets.empty()) {
      std::cout << "WARNING: " << empty_sets.size() << " empty sets in partition: ";
      for (Range::iterator rit = empty_sets.begin(); rit != empty_sets.end(); rit++)
        std::cout << *rit << " ";
      std::cout << std::endl;
    }
  }
  
  if (write_as_tags) {
      // allocate integer-size partitions
    result = mbImpl->tag_set_data(part_set_tag, elems, assignment); 
  }
  
  return MB_SUCCESS;
}

