#include "MBZoltan.hpp"
#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ReorderTool.hpp"

#ifdef CGM
#include "InitCGMA.hpp"
#include "CubitCompat.hpp"
#endif

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <list>
#include <time.h>

using namespace moab;

const char DEFAULT_ZOLTAN_METHOD[] = "RCB";
const char ZOLTAN_PARMETIS_METHOD[] = "PARMETIS";
const char ZOLTAN_OCTPART_METHOD[] = "OCTPART";

const char BRIEF_DESC[] = 
 "Use Zoltan to partition MOAB meshes for use on parallel computers";
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
  
  LONG_DESC << "This utility invokes the MBZoltan componemnt of MOAB/CGM"
               "to partition a mesh/geometry." << std::endl
            << "If no partitioning method is specified, the default is "
               "the Zoltan \"" << DEFAULT_ZOLTAN_METHOD << "\" method" << std::endl;
  
  ProgOptions opts(LONG_DESC.str(), BRIEF_DESC);
  opts.addOpt<int>( "dimension", "Specify dimension of entities to partition."
                                 "  Default is  largest in file.", 
                                 0, ProgOptions::int_flag );
  opts.addOpt<std::string>( "zoltan,z", "Specify Zoltan partition method.  "
                                         "One of RR, RCB, RIB, HFSC, PHG, " 
                                         " or Hypergraph (PHG and Hypergraph "
                                         "are synonymous)." );
  opts.addOpt<std::string>( "parmetis,p", "Specify Parmetis partition method.");
  opts.addOpt<std::string>( "octpart,o", "Specify OctPart partition method.");
  opts.addOpt<void>( "sets,s", "Write partition as tagged sets (Default)" );
  opts.addOpt<void>( "tags,t",  "Write partition by tagging entities");
  opts.addOpt<double>( "imbalance,i",  "Imbalance tolerance (used in PHG/Hypergraph method)");
  opts.addOpt<int>( "power,m", "Generate multiple partitions, in powers of 2, up to 2^(pow)");
  opts.addOpt<void>( "reorder,R", "Reorder mesh to group entities by partition");
  opts.addOpt<double>( "geom,g", "Specify if partition geometry and mesh size.");
  opts.addOpt<void>( "surf,f", "Specify if partition geometry surface.");
  opts.addOpt<void>( "ghost,h", "Specify if partition ghost geometry body.");
  opts.addOpt<int>( "vertex_w,v", "Number of weights associated with a graph vertex.");
  opts.addOpt<int>( "edge_w,e", "Number of weights associated with an edge.");
  opts.addRequiredArg<int>( "#parts", "Number of parts in partition" );
  opts.addRequiredArg<std::string>( "input_file", "Mesh/geometry to partition" );
  opts.addRequiredArg<std::string>( "output_file", "File to which to write partitioned mesh/geometry" );
  opts.addOpt<void>( ",T", "Print CPU time for each phase.");
  opts.parseCommandLine( argc, argv ); 

  MBZoltan *tool = NULL;
  double part_geom_mesh_size = -1.0;
  bool part_surf = false;
  bool ghost = false;
  int part_dim = -1;
  long num_parts;
  bool write_sets = true, write_tags = false;
  double imbal_tol = 1.10;
  int power = -1;
  int obj_weight = -1;
  int edge_weight = -1;
  
  bool print_time = opts.getOpt<void>(",T",0);
  
  // check if partition geometry, if it is, should get mesh size for the geometry
  if (opts.getOpt( "geom", &part_geom_mesh_size )) {
    if (part_geom_mesh_size < 0.0) {
      std::cerr << part_geom_mesh_size << ": invalid geometry partition mesh size." << std::endl;
      return 1;
    }
  }

  if (part_geom_mesh_size < 0.) { // partition mesh
    tool = new MBZoltan (&mb, false, argc, argv);
  }
  else { // partition geometry
#ifdef CGM
    CubitStatus status = InitCGMA::initialize_cgma();
    if (CUBIT_SUCCESS != status) {
      std::cerr << "CGM couldn't be initialized." << std::endl;
      return 1;
    }
    GeometryQueryTool *gti = GeometryQueryTool::instance();
    tool = new MBZoltan (&mb, false, argc, argv, gti);
#else
    std::cerr << "CGM should be configured to partition geometry." << std::endl;
    return 1;
#endif
  }

  part_surf = opts.numOptSet("surf") > 0;
  ghost = opts.numOptSet("ghost") > 0;

  std::string zoltan_method, other_method;
  if (!opts.getOpt("zoltan", &zoltan_method))
    zoltan_method = DEFAULT_ZOLTAN_METHOD;
  if (opts.getOpt("parmetis", &other_method))
    zoltan_method = ZOLTAN_PARMETIS_METHOD;
  if (opts.getOpt("octpart", &other_method))
    zoltan_method = ZOLTAN_OCTPART_METHOD;
  
  write_sets = opts.getOpt<void>("sets",0);
  write_tags = opts.getOpt<void>("tags",0);
  if (!write_sets && !write_tags)
    write_sets = true;

  if (!opts.getOpt("power",&power)) {
    num_parts = opts.getReqArg<int>("#parts");
    power = 1;
  }
  else if (power < 1 || power > 18) {
    std::cerr << power << ": invalid power for multiple patitions. Expected value in [1,18]" << std::endl;
    return 1;
  }
  else {
    num_parts = 2;
  }

  if (!opts.getOpt("vertex_w", &obj_weight)) obj_weight = 0;
  if (!opts.getOpt("edge_w", &edge_weight)) edge_weight = 0;

  if (opts.getOpt( "imbalance", &imbal_tol )) {
    if (imbal_tol < 0.0) {
      std::cerr << imbal_tol << ": invalid imbalance tolerance" << std::endl;
      return 1;
    }
  }
  
  std::string input_file = opts.getReqArg<std::string>("input_file");
  std::string output_file = opts.getReqArg<std::string>("output_file");
  clock_t t = clock();

  const char* options = NULL;
  if (part_geom_mesh_size > 0.) options = "FACET_DISTANCE_TOLERANCE=0.1";
  ErrorCode rval = mb.load_file( input_file.c_str(), 0, options );
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
  
  if (opts.getOpt("dimension", &part_dim)) {
    if (part_dim < 0 || part_dim > 3) {
      std::cerr << part_dim << " : invalid dimension" << std::endl;
      return 1;
    }
  } 
  else {
    for (int dim = 3; dim >= 0; --dim) {
      int n;
      rval = mb.get_number_entities_by_dimension( 0, dim, n );
      if (MB_SUCCESS == rval && 0 != n) {
        part_dim = dim;
        break;
      }
    }
  }
  if (part_dim < 0) {
    std::cerr << input_file << " : file does not contain any mesh entities" << std::endl;
    return 2;
  }
  
  ReorderTool reorder(&moab);
  for (int p = 0; p < power; p++) {
    t = clock();
    rval = tool->partition_mesh_geom( part_geom_mesh_size, num_parts, zoltan_method.c_str(), other_method.c_str(),
                                      imbal_tol, write_sets, write_tags, part_dim, obj_weight, edge_weight,
                                      part_surf, ghost );
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
    
    if (opts.getOpt<void>("reorder",0) && part_geom_mesh_size < 0.) {
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
        rval = reorder.handle_order_from_sets_and_adj( sets, order );
      }
      else {
        rval = reorder.handle_order_from_int_tag( tag, -1, order );
      }
      if (MB_SUCCESS != rval) {
        std::cerr << "Failed to calculate reordering!" << std::endl;
        return 2;
      }
      
      rval = reorder.reorder_entities( order );
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
  
    std::ostringstream tmp_output_file;
    
    if (power > 1) {
        // append num_parts to output filename
      std::string::size_type idx = output_file.find_last_of( "." );
      if (idx == std::string::npos) {
        tmp_output_file << output_file << "_" << num_parts;
        if (part_geom_mesh_size < 0.) tmp_output_file << ".h5m";
        else {
          std::cerr << "output file type is not specified." << std::endl;
          return 1;
        }
      }
      else {
        tmp_output_file << output_file.substr(0, idx) << "_" << num_parts
                        << output_file.substr(idx);
      }
    }
    else
      tmp_output_file << output_file;

    t = clock();
    if (part_geom_mesh_size < 0.) {
      rval = mb.write_file( tmp_output_file.str().c_str() );
      if (MB_SUCCESS != rval) {
        std::cerr << tmp_output_file << " : failed to write file." << std::endl;
        std::cerr << "  Error code: " << mb.get_error_string(rval) << " (" << rval << ")" << std::endl;
        std::string errstr;
        mb.get_last_error(errstr);
        if (!errstr.empty())
          std::cerr << "  Error message: " << errstr << std::endl;
        return 2;
      }
    }
    else {
#ifdef CGM
      std::string::size_type idx = output_file.find_last_of( "." );
      int c_size = output_file.length() - idx;
      const char* file_type = NULL;
      if (output_file.compare(idx, c_size, ".occ") == 0
          || output_file.compare(idx, c_size, ".OCC") == 0) file_type = "OCC";
      else if (output_file.compare(idx, c_size, ".sab") == 0) file_type = "ACIS_SAB";
      else if (output_file.compare(idx, c_size, ".sat") == 0) file_type = "ACIS_SAT";
      else {
        std::cerr << "File type for " << output_file.c_str() << " not supported." << std::endl;
        return 1;
      }

      int junk;
      DLIList<RefEntity*> ref_entity_list;
      CubitStatus status = CubitCompat_export_solid_model(ref_entity_list,
                           tmp_output_file.str().c_str(),
                           file_type, junk,
                           CubitString(__FILE__));
      if (CUBIT_SUCCESS != status) {
        std::cerr << "CGM couldn't export models." << std::endl;
        return 1;
      }
#endif
    }
    
    if (print_time)
      std::cout << "Wrote \"" << tmp_output_file.str() << "\" in "
                << (clock() - t)/(double)CLOCKS_PER_SEC << " seconds" 
                << std::endl;
                
    num_parts *= 2;
  }

  delete tool;
  
  return 0;    
}
