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

//Point Locater
extern "C" 
{
#include "moab/FindPtFuncs.h"
}

#include "moab/point_locater/point_locater.hpp"
#include "moab/point_locater/io.hpp"
#include "moab/point_locater/tree/element_tree.hpp"
#include "moab/point_locater/tree/bvh_tree.hpp"
#include "moab/point_locater/parametrizer.hpp"

//iMesh
#include "imesh/iMesh_extensions.h"
#include "imesh/MBiMesh.hpp"

//MOAB
#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "MBTagConventions.hpp"

//STL
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace io = moab::point_locator::io;

void print_usage() {
  std::cerr << "Usage: point_search" ;
  std::cerr << "-meshes <source_mesh> <target_mesh> ";
  std::cerr << " -itag <interp_tag> [-gnorm <gnorm_tag>] " ;
  std::cerr << " [-ssnorm <ssnorm_tag> <ssnorm_selection>] [-ropts <roptions>]";
  std::cerr << " [-outfile <out_file> [-wopts <woptions>]]"; 
  std::cerr << " [-dbgout [<dbg_file>]]" << std::endl;
  std::cerr << "    -meshes" << std::endl;
  std::cerr << "        Read in mesh files <source_mesh> and <target_mesh>." 
	    << std::endl;
  std::cerr << "    -itag" << std::endl;
  std::cerr << "        Interpolate tag <interp_tag> from source mesh to target"
	    << " mesh." << std::endl;
  std::cerr << "    -gnorm" << std::endl;
  std::cerr << "        Normalize the value of tag <gnorm_tag> over then entire"
	    << "mesh and save to" << std::endl;
  std::cerr << "        tag \"<gnorm_tag>_normf\" on the mesh set.  Do this for"
	    << " all meshes." << std::endl;
  std::cerr << "    -ssnorm" << std::endl;
  std::cerr << "        Normalize the value of tag <ssnorm_tag> over "
	    << "subsets of a mesh and save to" << std::endl;
  std::cerr << "        tag \"<ssnorm_tag>_normf\" on the Entity Set for each subset.  Subsets are selected" << std::endl;
  std::cerr << "        using criteria in <ssnorm_selection>.  Do this for all meshes." << std::endl;
  std::cerr << "    -ropts" << std::endl;
  std::cerr << "        Read in the mesh files using options in <roptions>." << std::endl;
  std::cerr << "    -outfile" << std::endl;
  std::cerr << "        Write out target mesh to <out_file>." << std::endl;
  std::cerr << "    -wopts" << std::endl;
  std::cerr << "        Write out mesh files using options in <woptions>." << std::endl;
  std::cerr << "    -dbgout" << std::endl;
  std::cerr << "        Write stdout and stderr streams to the file \'<dbg_file>.txt\'." << std::endl;
}

template< typename Pair>
struct Count_predicate: std::unary_function< Pair, bool>{
	bool operator()(const Pair & pair){ return pair.first != 0; }
};

// default types.. whatevs.
typedef moab::common_tree::Box< double> Bounding_box;
typedef std::vector< Bounding_box> Boxes;
typedef moab::EntityHandle Element;
typedef moab::element_utility::Parametrizer Parametrizer;
typedef Range Elements;
typedef moab::Element_tree< Elements, Bounding_box, moab::Core, Parametrizer> 
							Element_tree_;
typedef moab::Bvh_tree< Elements, Bounding_box, moab::Core, Parametrizer> 
							Bvh_tree_;
typedef moab::ParallelComm Communicator; 
typedef moab::Point_search< Element_tree_, std::vector< Bounding_box> > 
								Element_locater;
typedef moab::Point_search< Bvh_tree_, std::vector< Bounding_box> > 
								Bvh_locater;

template< typename Communicator, typename Moab, 
	  typename Root, typename String, typename String_>
void read_mesh( Communicator & comm, 
		Moab & moab, Root root, 
		String filename, String_ & options){
  moab::ErrorCode result = moab::MB_SUCCESS;
  std::stringstream ss;
  ss << options << ";PARALLEL_COMM=" << comm.get_id();
  std::string read_options( ss.str());
  std::cout << read_options << std::endl;
  result = moab.load_file( filename, 
			   root, 
			   read_options.c_str() );
  if(result != moab::MB_SUCCESS){ 
	std::string error;
	moab.get_last_error( error);
	std::cerr << error << std::endl;
	MPI_Abort(MPI_COMM_WORLD, result);
	std::exit( -1); 
  }
}

int main(int argc, char* argv[]){
 	MPI_Init(&argc, &argv);
 	io::File_options<> options;

	moab::ErrorCode result = moab::MB_SUCCESS;
	bool help = false;
	result = io::get_file_options(argc, argv, options);

	if (result != moab::MB_SUCCESS || help) {
 	  print_usage();
 	  MPI_Finalize();
 	  return 1;
 	}

	int nprocs, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	moab::Core moab;
	Communicator source_comm( &moab, MPI_COMM_WORLD);
	Communicator target_comm( &moab, MPI_COMM_WORLD);
	std::vector< moab::EntityHandle > roots( options.meshFiles.size());
	moab.create_meshset( moab::MESHSET_SET, roots[ 0]);
	moab.create_meshset( moab::MESHSET_SET, roots[ 1]);
	read_mesh( source_comm, moab, &roots[ 0], 
		   options.meshFiles[ 0].c_str(), options.readOpts);
	read_mesh( target_comm, moab, &roots[ 1], 
		   options.meshFiles[ 1].c_str(), options.readOpts);
	
	Range source_element_range, target_element_range; 
	Range target_vertex_handles, non_owned_vertices;
	source_comm.get_part_entities( source_element_range, 3);
	target_comm.get_part_entities( target_element_range, 3);
	moab.get_adjacencies( target_element_range, 0, false, 
			      target_vertex_handles, 
			      moab::Interface::UNION);
	target_comm.get_pstatus_entities( 0, PSTATUS_NOT_OWNED, 
					  non_owned_vertices);
	moab::subtract( target_vertex_handles, non_owned_vertices);
	moab::common_tree::Box< double> bounding_box;	
	Parametrizer p;
	//double construction_start = MPI_Wtime();
	//Element_tree_ tree_1(source_element_range, moab, bounding_box, p);
	//double construction = MPI_Wtime() - construction_start;
	//std::cout << "element construction time: " 
	//<< construction << std::endl;
	double construction_start = MPI_Wtime();
	Bvh_tree_ tree_2(source_element_range, moab, bounding_box, p);
	double construction = MPI_Wtime() - construction_start;
	std::cout << "bvh construction time: " << construction << std::endl;
	Boxes boxes;
	Bvh_locater locator( tree_2, boxes);
	typedef std::vector< double> Point;
	typedef std::pair< moab::EntityHandle, Point > Result;
	typedef Count_predicate< Result> Predicate;
	std::vector< Result > located_elements;
	std::vector< Point > target_vertices;
	target_vertices.reserve( target_vertex_handles.size());
	std::size_t index = 0;
	for( Range::iterator i = target_vertex_handles.begin();
			     i != target_vertex_handles.end();
					++i, ++index){
		std::vector< double> coordinate( 3, 0);
		moab.get_coords( &*i, 3, &coordinate[ 0]);
		target_vertices.push_back( coordinate);
	}
	const double tol=1e-6;
	std::cout << "point location beginning" << std::endl;
	double locate_start = MPI_Wtime();
	locator.locate_points( target_vertices, located_elements, tol); 
	double locate_points = MPI_Wtime() - locate_start;
	std::cout << "point_location: " << locate_points << std::endl;
	std::size_t unfound_points = std::count_if( located_elements.begin(),
						    located_elements.end(), 
						     Predicate());
	typedef  ct::_Element_data< const Bounding_box, double > Element_data;
	typedef std::tr1::unordered_map< Element, 
					 Element_data> Entity_map;
	if( unfound_points > 0){
		std::cout << "unlocated points: " <<  unfound_points
		<< "/" << located_elements.size() << std::endl;
		locator.bruteforce_locate_points( target_vertices, 
						  located_elements, tol);
	}else{
		std::cout << "all points located!" << std::endl;
	}
}
