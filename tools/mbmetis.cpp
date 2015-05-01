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
// Checkout mbzoltan for a more comprehensive set of partitioning options

#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ReorderTool.hpp"

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <list>
#include <time.h>
#include "moab/MetisPartitioner.hpp"

const char DEFAULT_TAGGEDSETS_TAG[] = "PARALLEL_PARTITION";
const char ALTERNATIVE_METIS_METHOD[] = "ML_RB";
const char DEFAULT_METIS_METHOD[] = "ML_KWAY";

using namespace moab;

const char BRIEF_DESC[] = 
 "Use Metis to partition MOAB meshes for use on parallel computers";
std::ostringstream LONG_DESC;

int main( int argc, char* argv[] )
{
  int err = MPI_Init( &argc, &argv );
  if (err) {
    std::cerr << "MPI_Init failed.  Aborting." << std::endl;
    return 3;
  }

  Core moab;
  Interface& mb = moab;
  std::vector<int> set_l;
  
  LONG_DESC << "This utility invokes the MetisPartitioner component of MOAB/CGM"
    "to partition a mesh/geometry." << std::endl;
  
  ProgOptions opts(LONG_DESC.str(), BRIEF_DESC);

  int part_dim = 3;
  opts.addOpt<int>( "dimension", "Specify dimension of entities to partition."
		    "  Default is  largest in file.", 
                    &part_dim, ProgOptions::int_flag );
  
  std::string metis_method;
  opts.addOpt<std::string>( "metis,m", "Specify Metis partition method.  "
                            "ML_RB or ML_KWAY", &metis_method);

  bool write_sets = true, write_tags = false;
  opts.addOpt<void>( "sets,s", "Write partition as tagged sets (Default)", &write_sets);
  opts.addOpt<void>( "tags,t",  "Write partition by tagging entities", &write_tags);

  bool reorder = false;
  opts.addOpt<void>( "reorder,R", "Reorder mesh to group entities by partition", &reorder);

  int num_parts;
  opts.addRequiredArg<int>( "#parts", "Number of parts in partition" );

  std::string input_file, output_file;
  opts.addRequiredArg<std::string>( "input_file", "Mesh/geometry to partition", &input_file );
  opts.addRequiredArg<std::string>( "output_file", "File to which to write partitioned mesh/geometry", &output_file );

  bool print_time = false;
  opts.addOpt<void>( ",T", "Print CPU time for each phase.", &print_time);
  
  bool partition_tagged_sets = false;
  opts.addOpt<void>( "taggedsets,x", "Partition tagged sets.", &partition_tagged_sets);
  
  bool partition_tagged_ents = false;
  opts.addOpt<void>( "taggedents,y", "Partition tagged ents.", &partition_tagged_ents);

  std::string aggregating_tag;
  opts.addOpt<std::string>( "aggregatingtag,a", "Specify aggregating tag to partion tagged sets or tagged entities.", &aggregating_tag);
  
  std::string aggregating_bc_tag;
  opts.addOpt<std::string>( "aggregatingBCtag,B", "Specify boundary id tag name used to group cells with same boundary ids.", &aggregating_bc_tag);

  std::string boundaryIds;
  std::vector<int> BCids;
  opts.addOpt<std::string>( "aggregatingBCids,I", " Specify id or ids of boundaries to be aggregated before partitioning (all elements with same boundary id will be in the same partition). Comma separated e.g. -I 1,2,5 ", &boundaryIds);

  opts.parseCommandLine( argc, argv ); 

  MetisPartitioner *tool = NULL;

  tool = new MetisPartitioner (&mb, false, argc, argv);
  num_parts = opts.getReqArg<int>("#parts");

  if ((aggregating_tag.empty() && partition_tagged_sets) || (aggregating_tag.empty() && partition_tagged_ents))
    aggregating_tag = DEFAULT_TAGGEDSETS_TAG;
  if (!write_sets && !write_tags)
    write_sets = true;

  if (!boundaryIds.empty())
  {
    std::vector<std::string> ids;
    std::stringstream ss(boundaryIds);
    std::string item;
    while (std::getline(ss, item, ',')) {
       ids.push_back(item);
    }
    for (unsigned int i = 0; i < ids.size(); i++)
      BCids.push_back(std::atoi(ids[i].c_str()));
  }

  if (part_dim < 0 || part_dim > 3) {
    std::cerr << part_dim << " : invalid dimension" << std::endl;
    return 1;
  }

  if (metis_method.empty()) {
    metis_method = DEFAULT_METIS_METHOD;
  }
  
  clock_t t = clock();

  const char* options = NULL;
  ErrorCode rval;
  std::cout << "Loading file " << input_file << "..." << std::endl;

  int dotindexinput = input_file.find_last_of("."); 
  std::string extensioninput = input_file.substr(dotindexinput+1); 
  t = clock();
  rval = mb.load_file( input_file.c_str(), 0, options );
  if (MB_SUCCESS != rval) {
    std::cerr << input_file << " : failed to read file." << std::endl;
    std::cerr << "  Error code: " << mb.get_error_string(rval) << " (" << rval << ")" << std::endl;
    std::string errstr;
    mb.get_last_error(errstr);
    if (!errstr.empty())
      std::cerr << "  Error message: " << errstr << std::endl;
    return 2;
  }
  if (print_time)
    std::cout << "Read input file in " << (clock() - t)/(double)CLOCKS_PER_SEC << " seconds" << std::endl;
  
  for (int dim = part_dim; dim >= 0; --dim) {
    int n;
    rval = mb.get_number_entities_by_dimension( 0, dim, n );
    if (MB_SUCCESS == rval && 0 != n) {
      part_dim = dim;
      break;
    }
  }
  if (part_dim < 0) {
    std::cerr << input_file << " : file does not contain any mesh entities" << std::endl;
    return 2;
  }
 
  ReorderTool reorder_tool(&moab);
  t = clock();
  rval = tool->partition_mesh_geom( num_parts, metis_method.c_str(), part_dim,
                                    write_sets, write_tags, 
				    partition_tagged_sets, partition_tagged_ents, 
        			    aggregating_tag.c_str());
  if (MB_SUCCESS != rval) {
    std::cerr << "Partitioner failed!" << std::endl;
    std::cerr << "  Error code: " << mb.get_error_string(rval) << " (" << rval << ")" << std::endl;
    std::string errstr;
    mb.get_last_error(errstr);
    if (!errstr.empty())
      std::cerr << "  Error message: " << errstr << std::endl;
    return 3;
  }
  if (print_time) 
    std::cout << "Generated " << num_parts << " part partitioning in "
                << (clock() - t)/(double)CLOCKS_PER_SEC << " seconds" 
                << std::endl;
  
  if (reorder) {
    std::cout << "Reordering mesh for partition..." << std::endl;

    Tag tag, order;
    rval = mb.tag_get_handle( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, tag );
    if (MB_SUCCESS != rval) {
      std::cerr << "Partitioner did not create PARALLEL_PARTITION tag" << std::endl;
      return 2;
    }
    
    t = clock();
    if (write_sets) {
      Range sets;
      mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, 0, 1, sets );
      rval = reorder_tool.handle_order_from_sets_and_adj( sets, order );
    }
    else {
      rval = reorder_tool.handle_order_from_int_tag( tag, -1, order );
    }
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to calculate reordering!" << std::endl;
      return 2;
    }
    
    rval = reorder_tool.reorder_entities( order );
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to perform reordering!" << std::endl;
      return 2;
    }
    
    mb.tag_delete( order );
    if (print_time) 
      std::cout << "Reorderd mesh in "
                << (clock() - t)/(double)CLOCKS_PER_SEC << " seconds" 
                << std::endl;
  }
  
  int dotindex = output_file.find_last_of("."); 
  std::string extension = output_file.substr(dotindex+1); 
  t = clock();
  rval = mb.write_file( output_file.c_str() );
  if (MB_SUCCESS != rval) {
    std::cerr << output_file << " : failed to write file." << std::endl;
    std::cerr << "  Error code: " << mb.get_error_string(rval) << " (" << rval << ")" << std::endl;
    std::string errstr;
    mb.get_last_error(errstr);
    if (!errstr.empty())
      std::cerr << "  Error message: " << errstr << std::endl;
    return 2;
  }
  
  if (print_time)
    std::cout << "Wrote \"" << output_file << "\" in "
              << (clock() - t)/(double)CLOCKS_PER_SEC << " seconds" 
              << std::endl;
              
  delete tool;
  return 0;    
}
