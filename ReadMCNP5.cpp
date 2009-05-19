/**
 * \class ReadMCNP5
 * \brief  Read output from MCNP5
 * \author Brandon Smith
 **/
#include "ReadMCNP5.hpp"
#include "MBInterface.hpp"
#include "MBReadUtilIface.hpp"
#include "MBInternals.hpp" // for MB_START_ID
#include "MBRange.hpp"
#include "FileOptions.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "assert.h"
#include "math.h"


// this setup is mostly copied from ReadSms.cpp
MBReaderIface* ReadMCNP5::factory( MBInterface* iface ) { 
  return new ReadMCNP5( iface );
}

ReadMCNP5::ReadMCNP5(MBInterface* impl)
  : MBI(impl) {
    assert( NULL!=impl);
    void *ptr = 0;
    MBI->query_interface("MBReadUtilIface", &ptr);
    assert( NULL!=ptr );
    readMeshIface = reinterpret_cast<MBReadUtilIface*>(ptr);
}

ReadMCNP5::~ReadMCNP5() {
  if (readMeshIface) {
    MBI->release_interface("MBReadUtilIface", readMeshIface);
    readMeshIface = 0;
  }
}

MBErrorCode ReadMCNP5::load_file(const char* fname, 
                                 MBEntityHandle& input_meshset, 
                                 const FileOptions& options,
                                 const char* set_tag_name,       // not used
				 const int* material_set_list,   // not used
                                 const int num_material_sets ) { // not used

  std::cout << "begin MCNP5 reader" << std::endl;
  
  // no support for reading a subset of the file
  if (set_tag_name) {
    readMeshIface->report_error( "Reading subset of files not supported for meshtal." );
    return MB_UNSUPPORTED_OPERATION;
  }

  bool debug=true;
  MBErrorCode result;
  std::fstream file;
  file.open( fname, std::fstream::in );
  char line[10000];

  // Options specify if this file will be averaged with a file already existing in
  // the interface.
  bool average;
  if ( MB_SUCCESS ==options.get_null_option("AVERAGE_TALLY") )
    average = true;
  else
    average = false;

  // Create tags
  MBTag date_and_time_tag,  
        title_tag,           
	nps_tag,  
	tally_number_tag,    
	tally_comment_tag, 
        tally_particle_tag, 
	tally_coord_sys_tag, 
	tally_tag, 
	error_tag;

  result = create_tags( date_and_time_tag,   
                        title_tag,         
			nps_tag,            
                        tally_number_tag,    
			tally_comment_tag, 
			tally_particle_tag, 
			tally_coord_sys_tag, 
			tally_tag,         
			error_tag );
    assert(MB_SUCCESS == result);

  // ******************************************************************
  // This info exists only at the top of each meshtal file
  // ******************************************************************

  // define characteristics of the entire file
  char date_and_time[100] = "";
  char title[100] = "";
  // this file's number of particles
  unsigned long long int nps;
  // sum of this file's and existing file's nps for averaging
  unsigned long long int new_nps;

  // read the file header
  result = read_file_header( file, 
                             debug, 
                             date_and_time, 
			     title,
			     nps );
    assert(MB_SUCCESS == result);

  // blank line
  file.getline(line, 10000);

  // Everything stored in the file being read will be in the output_meshset.
  // if this is being saved in MOAB, set header tags
  MBEntityHandle output_meshset;
  if (!average) {
    result = MBI->create_meshset( MESHSET_SET, output_meshset );
      assert(MB_SUCCESS == result);
    result = set_header_tags( output_meshset, 
                              date_and_time,
                              title,
		              nps,
			      date_and_time_tag,
                              title_tag,
		              nps_tag );
      assert(MB_SUCCESS == result);
  }

  // ******************************************************************
  // This info is repeated for each tally in the meshtal file.
  // ******************************************************************

  // If averaging, nps will hold the sum of particles simulated in both tallies.
  while( !file.eof() ) {

    // define characteristics of this tally
    unsigned int        tally_number;
    char                tally_comment[100] = "";
    particle            tally_particle;
    coordinate_system   tally_coord_sys;
    std::vector<double> planes[3]; 
    int                 n_chopped_x2_planes;
    
    // read tally header
    result = read_tally_header( file, 
                                debug, 
                                tally_number, 
				tally_comment,
                                tally_particle );
      assert(MB_SUCCESS == result);
    
    // blank line
    file.getline(line, 10000);

    // read mesh planes
    result = read_mesh_planes( file, 
                               debug, 
			       planes, 
			       tally_coord_sys );
      assert(MB_SUCCESS == result);

    // get energy boundaries
    file.getline(line, 10000);
    std::string a = line;
    if (debug) std::cout << "Energy bin boundaries:=| " << a << std::endl;

    // blank
    file.getline(line, 10000);

    // column headers
    file.getline(line, 10000);

    // If using cylidrical mesh, it may be necessary to chop off the last theta element.
    // We chop off the last theta plane because the elements will be wrong and skew up
    // the tree building code. This is
    // because the hex elements are a linear approximation to the cylindrical elements.
    // Chopping off the last plane is problem-dependent, and due to MCNP5's mandate 
    // that the cylidrical mesh must span 360 degrees.
    if ( CYLINDRICAL==tally_coord_sys &&
         MB_SUCCESS ==options.get_null_option("REMOVE_LAST_CYLIDRICAL_PLANE") ) {
      planes[2].pop_back();
      n_chopped_x2_planes = 1;
      if (debug) std::cout << "remove last cylindrical plane option found" << std::endl;
    } else {
      n_chopped_x2_planes = 0;
    }
    
    // read the values and errors of each element from the file.
    // Do not read values that are chopped off.
    unsigned int n_elements = (planes[0].size()-1) * (planes[1].size()-1) * (planes[2].size()-1);
    double *values, *errors;
    values = new double [n_elements];
    errors = new double [n_elements];
    result = read_element_values_and_errors( file, 
                                             debug, 
					     planes,
					     n_chopped_x2_planes,
					     values, 
					     errors );
      assert(MB_SUCCESS == result);
    
    // blank line
    file.getline(line, 10000);
   
    // ****************************************************************
    // This tally has been read. If it is not being averaged, build tags,
    // vertices
    // and elements. If it is being averaged, average the data with a
    // tally already existing in the MOAB instance.
    // ****************************************************************
    if (!average) {
      MBEntityHandle tally_meshset;
      result = MBI->create_meshset(MESHSET_SET, tally_meshset);
        assert(MB_SUCCESS == result);
      if (MB_SUCCESS != result) return result;
      
      // set tags on the tally
      result = set_tally_tags( tally_meshset,
                               tally_number,
			       tally_comment,
			       tally_particle,
			       tally_coord_sys,
                               tally_number_tag, 
			       tally_comment_tag,
                               tally_particle_tag,
                               tally_coord_sys_tag );
        assert(MB_SUCCESS == result);

      // The only info needed to build elements is the mesh plane boundaries.
      // Build vertices...
      //MBRange vert_handles;
      //result = create_vertices( planes, 
      //                          vert_handles, 
      //        	          tally_coord_sys,
      //                          tally_meshset );
      MBEntityHandle start_vert = 0;
      result = create_vertices( planes, 
                                start_vert, 
                                tally_coord_sys,
                                tally_meshset );
        assert(MB_SUCCESS == result); 
      
      // Build elements and tag them with tally values and errors, then add
      // them to the tally_meshset.
      result = create_elements( debug, 
                                planes, 
	  		        n_chopped_x2_planes,
                                //vert_handles, 
                                start_vert, 
			        values, 
			        errors, 
			        tally_tag, 
			        error_tag, 
			        tally_meshset );
      assert(MB_SUCCESS == result); 
      
      // add this tally's meshset to the output meshset
      std::cout << "not averaging tally" << std::endl;
      result = MBI->add_entities( output_meshset, &tally_meshset, 1);
        if (MB_SUCCESS != result) std::cout << "error=" << result << std::endl;
        assert(MB_SUCCESS == result);

    // average the tally values, then delete stuff that was created
    } else {
      std::cout << "averaging tally" << std::endl;
      result = average_with_existing_tally( new_nps,
                                            nps,
					    tally_number,
				            tally_number_tag,
				            nps_tag,
				            tally_tag,
				            error_tag,
					    values,
					    errors,
					    n_elements );
        assert( MB_SUCCESS == result );
    }

    // clean up
    delete[] values;
    delete[] errors;
  }

  // If we are averaging, delete the remainder of this file's information.
  // Add the new nps to the existing file's nps if we are averaging.
  // This is calculated during every tally averaging but only used after the last one.
  if (average) {
    MBRange matching_nps_sets;
    result = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &nps_tag, 
                                                0, 1, matching_nps_sets );
      assert(MB_SUCCESS == result);
    std::cout << "number of matching nps  meshsets=" << matching_nps_sets.size() << std::endl;
    assert( 1 == matching_nps_sets.size() );
    result = MBI->tag_set_data( nps_tag, matching_nps_sets, &new_nps );
      if (MB_SUCCESS != result) std::cout<< "result=" << result << std::endl;
      assert(MB_SUCCESS == result);

  // If this file is not being averaged, add it to the meshset that was input after we create it.
  } else {
    result = MBI->create_meshset( MESHSET_SET, input_meshset );
      assert(MB_SUCCESS == result);
    result = MBI->add_entities( input_meshset, &output_meshset, 1 );
      assert(MB_SUCCESS == result);
  }

  file.close();
  return MB_SUCCESS;
}

// Create tags needed for this reader
MBErrorCode ReadMCNP5::create_tags( MBTag &date_and_time_tag,
                                    MBTag &title_tag,
                                    MBTag &nps_tag,
				    MBTag &tally_number_tag,
				    MBTag &tally_comment_tag,
				    MBTag &tally_particle_tag,
				    MBTag &tally_coord_sys_tag,
				    MBTag &tally_tag,
				    MBTag &error_tag ) {
  MBErrorCode result;
  result = MBI->tag_create("DATE_AND_TIME_TAG", sizeof(char[100]), MB_TAG_SPARSE, 
                           MB_TYPE_OPAQUE, date_and_time_tag, 0);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  result = MBI->tag_create("TITLE_TAG", sizeof(char[100]), MB_TAG_SPARSE, 
                           MB_TYPE_OPAQUE, title_tag, 0);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  result = MBI->tag_create("NPS_TAG", sizeof(unsigned long long int), 
                           MB_TAG_SPARSE, MB_TYPE_OPAQUE, nps_tag, 0);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  result = MBI->tag_create("TALLY_NUMBER_TAG", sizeof(int), MB_TAG_SPARSE, 
                           MB_TYPE_INTEGER, tally_number_tag, 0);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  result = MBI->tag_create("TALLY_COMMENT_TAG", sizeof(char[100]), MB_TAG_SPARSE,
                           MB_TYPE_OPAQUE, tally_comment_tag, 0);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  result = MBI->tag_create("TALLY_PARTICLE_TAG", sizeof(particle), MB_TAG_SPARSE,
                           MB_TYPE_OPAQUE, tally_particle_tag, 0);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  result = MBI->tag_create("TALLY_COORD_SYS_TAG", sizeof(coordinate_system), MB_TAG_SPARSE,
                           MB_TYPE_OPAQUE, tally_coord_sys_tag, 0);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  result = MBI->tag_create("TALLY_TAG", sizeof(double), MB_TAG_DENSE, 
                           MB_TYPE_DOUBLE, tally_tag, 0);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  result = MBI->tag_create("ERROR_TAG", sizeof(double), MB_TAG_DENSE, 
                           MB_TYPE_DOUBLE, error_tag, 0);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  return MB_SUCCESS;
}

MBErrorCode ReadMCNP5::read_file_header( std::fstream           &file,
                                         bool                   debug,
                                         char                   date_and_time[100],  
				         char                   title[100], 
				         unsigned long long int &nps ) {

  // get simulation date and time
  // mcnp   version 5     ld=11242008  probid =  03/23/09 13:38:56
  char line[100];
  file.getline(line, 100);
  date_and_time = line;
  if (debug) std::cout << "date_and_time=| " << date_and_time << std::endl;

  // get simulation title
  // iter Module 4                                                                   
  file.getline(line, 100);
  title = line;
  if (debug) std::cout << "title=| " << title  << std::endl;

  // get number of histories
  // Number of histories used for normalizing tallies =      50000000.00
  file.getline(line, 100);
  std::string a = line;
  unsigned int b = a.find("Number of histories used for normalizing tallies =");
  if (std::string::npos!=b) {
    //std::cout << b << std::endl;
    //std::string nps_string = a.substr(b+52,100);
    //std::cout << nps_string << std::endl;
    std::istringstream nps_ss( 
      a.substr(b+sizeof("Number of histories used for normalizing tallies = "),100) );
    nps_ss >> nps;
    if (debug) std::cout << "nps=| " << nps << std::endl;
  } else return MB_FAILURE;

  return MB_SUCCESS;
}


MBErrorCode ReadMCNP5::set_header_tags( MBEntityHandle         output_meshset, 
                                        char                   date_and_time[100],
                                        char                   title[100],
				        unsigned long long int nps,
				        MBTag                  data_and_time_tag,
                                        MBTag                  title_tag,
				        MBTag                  nps_tag ) {
  MBErrorCode result;
  result = MBI->tag_set_data( data_and_time_tag, &output_meshset, 1, &date_and_time);
    assert(MB_SUCCESS == result);
  result = MBI->tag_set_data( title_tag, &output_meshset, 1, &title);
    assert(MB_SUCCESS == result);
  result = MBI->tag_set_data( nps_tag, &output_meshset, 1, &nps);
  if(MB_SUCCESS != result) std::cout << "result=" << result << std::endl;
    assert(MB_SUCCESS == result);
  return MB_SUCCESS;
}

MBErrorCode ReadMCNP5::read_tally_header( std::fstream   &file,
                                          bool           debug,
					  unsigned int   &tally_number,
					  char           tally_comment[100],
					  particle       &tally_particle ) {

  // get tally number
  // Mesh Tally Number 104
  MBEntityHandle result;
  char line[100];
  file.getline(line, 100);
  std::string a = line;
  unsigned int b = a.find("Mesh Tally Number");
  if (std::string::npos != b) {
    //std::istringstream tally_number_ss( a.substr(b+18,100) );
    std::istringstream tally_number_ss( a.substr(b+sizeof("Mesh Tally Number "),100) );
    tally_number_ss >> tally_number;
    if (debug) std::cout << "tally_number=| " << tally_number << std::endl;
  } else return MB_FAILURE;
     
  // next get the tally comment (optional) and particle type
  // 3mm neutron heating in Be (W/cc)     
  // This is a neutron mesh tally.
  // std::string tally_comment;
  
  // get tally particle
  file.getline(line, 100);
  a = line;
  result = get_tally_particle(a, debug, tally_particle);
  if (MB_FAILURE == result) {
    // If this line does not specify the particle type, then it is a tally comment.
    // Get the comment, then get the particle type from the next line.
    tally_comment = line;
    file.getline(line, 100);
    a = line;
    result = get_tally_particle(a, debug, tally_particle);
      assert(MB_SUCCESS==result);
  }
  if (debug) std::cout << "tally_comment=| " << tally_comment << std::endl;
  return MB_SUCCESS;
}

MBErrorCode ReadMCNP5::get_tally_particle( std::string    a,
                                           bool           debug,
					   particle       &tally_particle ) {
  
  if        (std::string::npos != a.find("This is a neutron mesh tally.")) {
    tally_particle = NEUTRON;
  } else if (std::string::npos != a.find("This is a photon mesh tally.") ) {
    tally_particle = PHOTON;
  } else if (std::string::npos != a.find("This is an electron mesh tally.") ) {
    tally_particle = ELECTRON;
  } else return MB_FAILURE;

  if (debug) std::cout << "tally_particle=| " << tally_particle << std::endl;
  return MB_SUCCESS;
}

MBErrorCode ReadMCNP5::read_mesh_planes( std::fstream         &file, 
                                         bool                 debug, 
					 std::vector<double>  planes[3], 
					 coordinate_system    &coord_sys ) {

  // Tally bin boundaries:
  MBErrorCode result;
  char line[10000];
  file.getline(line, 10000);
  std::string a = line;
  if (std::string::npos == a.find("Tally bin boundaries:"))
    return MB_FAILURE;
 
  // decide what coordinate system the tally is using
  // First check for Cylindrical coordinates:
  file.getline(line, 10000);
  a = line;
  unsigned int b = a.find("Cylinder origin at");
  if (std::string::npos != b) {
    coord_sys = CYLINDRICAL;
    if (debug) std::cout << "origin, axis, direction=| " << a << std::endl;
    std::istringstream ss(a.substr(b+sizeof("Cylinder origin at"),10000));
    // get origin (not used)
    // Cylinder origin at   0.00E+00  0.00E+00  0.00E+00, axis in  0.000E+00 0.000E+00 1.000E+00 direction
    double origin[3];
    std::cout << "origin=| ";
    for (int i=0; i<3; i++) {
      ss >> origin[i];
      std::cout << origin[i] << " ";
    }
    std::cout << std::endl;
    int length_of_string = 10;
    ss.ignore( length_of_string, ' ');
    ss.ignore( length_of_string, ' ');
    ss.ignore( length_of_string, ' ');
    // get acis (not used)
    double axis[3];
    std::cout << "axis=| ";
    for (int i=0; i<3; i++) {
      ss >> axis[i];
      std::cout << axis[i] << " ";
    }
    std::cout << std::endl;
    file.getline(line, 10000);
    a = line;

    // get r planes
    if (debug) std::cout << "R direction:=| "; // << a << std::endl;
    b = a.find("R direction:");
    if (std::string::npos != b) {
      std::istringstream ss(a.substr(b+sizeof("R direction"),10000));
      result = get_mesh_plane( ss, planes[0] );
       assert(MB_SUCCESS == result);
    } else return MB_FAILURE;

    // get z planes
    file.getline(line, 10000);
    a = line;
    if (debug) std::cout << "Z direction:=| "; // << a << std::endl;
    b = a.find("Z direction:");
    if (std::string::npos != b) {
      std::istringstream ss(a.substr(b+sizeof("Z direction"),10000));
      result = get_mesh_plane( ss, planes[1] );
        assert(MB_SUCCESS == result);
    } else return MB_FAILURE;

    // get theta planes
    file.getline(line, 10000);
    a = line;
    if (debug) std::cout << "Theta direction:=| "; // << a << std::endl;
    b = a.find("Theta direction (revolutions):");
    if (std::string::npos != b) {
      std::istringstream ss(a.substr(b+sizeof("Theta direction (revolutions):"),10000));
      result = get_mesh_plane( ss, planes[2] );
        assert(MB_SUCCESS == result);
    } else return MB_FAILURE;
    
  // Cartesian coordinate system:
  }  else if (std::string::npos != a.find("X direction:") ) {
    coord_sys = CARTESIAN;
    // get x planes
    b = a.find("X direction:");
    if (std::string::npos != b) {
      std::istringstream ss(a.substr(b+sizeof("X direction"),10000));
      result = get_mesh_plane( ss, planes[0] );
        assert(MB_SUCCESS == result);
    } else return MB_FAILURE;

    // get y planes
    file.getline(line, 10000);
    a = line;
    if (debug) std::cout << "Y direction:=| "; // << a << std::endl;
    b = a.find("Y direction:");
    if (std::string::npos != b) {
      std::istringstream ss(a.substr(b+sizeof("Y direction"),10000));
      result = get_mesh_plane( ss, planes[1] );
        assert(MB_SUCCESS == result);
    } else return MB_FAILURE;

    // get z planes
    file.getline(line, 10000);
    a = line;
    if (debug) std::cout << "Z direction:=| "; // << a << std::endl;
    b = a.find("Z direction:");
    if (std::string::npos != b) {
      std::istringstream ss(a.substr(b+sizeof("Z direction"),10000));
      result = get_mesh_plane( ss, planes[2] );
        assert(MB_SUCCESS == result);
    } else return MB_FAILURE;

  // Spherical coordinate system not yet implemented:
  } else return MB_FAILURE;
    
  return MB_SUCCESS;
}

// Given a stringstream, return a vector of values in the string.
MBErrorCode ReadMCNP5::get_mesh_plane( std::istringstream &ss, 
                                       std::vector<double> &plane) {
  double value;
  plane.clear();
  while (!ss.eof()) {
    ss >> value;
    plane.push_back( value );
    std::cout << value << " ";
  }
  std::cout << std::endl;
  return MB_SUCCESS;
}

MBErrorCode ReadMCNP5::read_element_values_and_errors( std::fstream        &file,
                                                       bool                debug,
						       std::vector<double> planes[3],
						       int                 n_chopped_x2_planes,
						       double              values[],
						       double              errors[] ) {

  unsigned int index = 0;
  for (unsigned int i=0; i<planes[0].size()-1; i++) {
    for (unsigned int j=0; j<planes[1].size()-1; j++) {
      for (unsigned int k=0; k<planes[2].size()-1+n_chopped_x2_planes; k++) {
          
	// need to read every line in the file, even if we chop off some elements
        char line[100];
        file.getline(line, 100);
	
	// if this element has been chopped off, skip it
	if (k>=planes[2].size()-1 && k<planes[2].size()-1+n_chopped_x2_planes) continue;
        std::string a=line;
        std::stringstream ss(a);
        double centroid[3];
        ss >> centroid[0];
        ss >> centroid[1];
        ss >> centroid[2];
	
        ss >> values[index];
        ss >> errors[index];
	index++;
      }
    }
  }
  return MB_SUCCESS;
}

MBErrorCode ReadMCNP5::set_tally_tags( MBEntityHandle    tally_meshset,
                                       unsigned int      tally_number,
			               char              tally_comment[100],
			               particle          tally_particle,
			               coordinate_system tally_coord_sys,
                                       MBTag             tally_number_tag, 
			               MBTag             tally_comment_tag,
                                       MBTag             tally_particle_tag,
                                       MBTag             tally_coord_sys_tag ) {
  MBErrorCode result;
  result = MBI->tag_set_data( tally_number_tag,    &tally_meshset, 1, &tally_number);
    assert(MB_SUCCESS == result);
  result = MBI->tag_set_data( tally_comment_tag,   &tally_meshset, 1, &tally_comment);
    assert(MB_SUCCESS == result);
  result = MBI->tag_set_data( tally_particle_tag,  &tally_meshset, 1, &tally_particle);
    assert(MB_SUCCESS == result);
  result = MBI->tag_set_data( tally_coord_sys_tag, &tally_meshset, 1, &tally_coord_sys);
    assert(MB_SUCCESS == result);
  return MB_SUCCESS;
}


//MBErrorCode ReadMCNP5::create_vertices( std::vector<double> planes[3],
//                                          MBRange           &vert_handles,
//					  coordinate_system coord_sys,
//					  MBEntityHandle    tally_meshset) {
MBErrorCode ReadMCNP5::create_vertices( std::vector<double> planes[3],
                                          MBEntityHandle    &start_vert,
					  coordinate_system coord_sys,
					  MBEntityHandle    tally_meshset) {
                                         
  // The only info needed to build elements is the mesh plane boundaries.
  MBErrorCode result;
  unsigned int n_verts = planes[0].size() * planes[1].size() * planes[2].size();
  //int n_verts = planes[0].size() * planes[1].size() * planes[2].size();

  std::vector<double*> coord_arrays;
  int start_id = 10000; // MB_START_ID
  result = readMeshIface->get_node_arrays( 3, n_verts, start_id, start_vert, coord_arrays );
    assert( MB_SUCCESS == result );

  //double *coords;
  //coords = new double [3*n_verts];
  for (unsigned int k=0; k < planes[2].size(); k++) {
    for (unsigned int j=0; j < planes[1].size(); j++) {
      for (unsigned int i=0; i < planes[0].size(); i++) {
        //unsigned int ijk = 3*(k*planes[0].size()*planes[1].size() + j*planes[0].size() + i);
        unsigned int idx = (k*planes[0].size()*planes[1].size() + j*planes[0].size() + i);
        double in[3], out[3];

        in[0] = planes[0][i];
        in[1] = planes[1][j];
        in[2] = planes[2][k];
        result = transform_point_to_cartesian( in, out, coord_sys );
          assert( MB_SUCCESS == result );
        //coords[ ijk   ] = out[0];
        //coords[ ijk+1 ] = out[1];
	//coords[ ijk+2 ] = out[2];
	*(coord_arrays[0]+idx) = out[0];
	*(coord_arrays[1]+idx) = out[1];
	*(coord_arrays[2]+idx) = out[2];
      }
    }
  }
  //result = MBI->create_vertices(coords, n_verts, vert_handles);
  //  assert( MB_SUCCESS == result );
  //delete[] coords;
 
  // add the vertices to the tally_meshset
  //result = MBI->add_entities( tally_meshset, vert_handles );
  //  assert( MB_SUCCESS == result );
  result = MBI->add_entities( tally_meshset, &start_vert, n_verts );
    assert( MB_SUCCESS == result );
  return MB_SUCCESS;
}


MBErrorCode ReadMCNP5::create_elements( bool                debug, 
                                        std::vector<double> planes[3],
					int                 n_chopped_x2_planes,
					//MBRange             vert_handles,
					MBEntityHandle      start_vert,
					double              values[],
					double              errors[],
					MBTag               tally_tag,
					MBTag               error_tag,
					MBEntityHandle      tally_meshset ) {
  MBErrorCode result;
  //MBEntityHandle connect[8], start_vert, index, hex_handle;
  //MBEntityHandle connect[8], index, hex_handle;
  //MBRange element_handles;
  //start_vert = *(vert_handles.begin());
  unsigned int index, idx;

  MBEntityHandle start_hex, *connect;
  unsigned int n_elements = (planes[0].size()-1) * (planes[1].size()-1) * (planes[2].size()-1);
  result = readMeshIface->get_element_array( n_elements, 8, MBHEX, 1, start_hex, connect );
    assert( MB_SUCCESS == result );



  for (unsigned int i=0; i<planes[0].size()-1; i++) {
    for (unsigned int j=0; j<planes[1].size()-1; j++) {
      for (unsigned int k=0; k<planes[2].size()-1+n_chopped_x2_planes; k++) {
          
        idx = i + j*planes[0].size() + k*planes[0].size()*planes[1].size();
        index = start_vert + idx;

        // if this element has been chopped off, skip it
	if (k>=planes[2].size()-1 && k<planes[2].size()-1+n_chopped_x2_planes) continue;

        // set connectivity and create the element
        connect[0] = index;
        connect[1] = index + 1;  
	connect[2] = index + 1 + planes[0].size();	      
	connect[3] = index +     planes[0].size();
	connect[4] = index +                        planes[0].size()*planes[1].size();
	connect[5] = index + 1 +                    planes[0].size()*planes[1].size();    
	connect[6] = index + 1 + planes[0].size() + planes[0].size()*planes[1].size();
	connect[7] = index +     planes[0].size() + planes[0].size()*planes[1].size();
        //result = MBI->create_element(MBHEX, connect, 8, hex_handle);
        //  assert(MB_SUCCESS == result);
	//if (MB_SUCCESS != result) return MB_FAILURE;
        //element_handles.insert(hex_handle);
 
        // assign parsed data to the element
        //result = MBI->tag_set_data(tally_tag, &hex_handle, 1, &values[idx]);
        result = MBI->tag_set_data(tally_tag, &index, 1, &values[idx]);
	  assert(MB_SUCCESS == result);
        //result = MBI->tag_set_data(error_tag, &hex_handle, 1, &errors[idx]);
        result = MBI->tag_set_data(error_tag, &index, 1, &errors[idx]);
          assert(MB_SUCCESS == result);
      }
    }
  }
  
  // add the elements to the tally set
  //result = MBI->add_entities( tally_meshset, element_handles );
  result = MBI->add_entities( tally_meshset, &start_hex, n_elements );
    assert( MB_SUCCESS == result );
  //std::cout << "Read " << element_handles.size() << " elements from tally." << std::endl;
  std::cout << "Read " << n_elements << " elements from tally." << std::endl;
  return MB_SUCCESS;
}

// Average a tally that was recently read in with one that alread exists in
// the interface. Only the existing values will be updated.
MBErrorCode ReadMCNP5::average_with_existing_tally( unsigned long long int &new_nps,
                                                    unsigned long long int nps1,
                                                    unsigned int  tally_number,
						     MBTag        tally_number_tag,
						     MBTag        nps_tag,
						     MBTag        tally_tag,
						     MBTag        error_tag,
						     double       values1[],
						     double       errors1[],
						     unsigned int n_elements ) {
    
  // get the tally number
  MBErrorCode result;
  //MBTag tally_number_tag, nps_tag, tally_tag, error_tag;
  //result = MBI->tag_create("TALLY_NUMBER_TAG", sizeof(int), MB_TAG_DENSE, 
  //                         MB_TYPE_INTEGER, tally_number_tag, 0);
  //  assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  //int tally_number;
  //result = MBI->tag_get_data( tally_number_tag, &tally_meshset, 1, &tally_number);
  //  assert(MB_SUCCESS == result);

  // match the tally number with one from the existing meshtal file
  MBRange matching_tally_number_sets;
  const void* const tally_number_val[] = {&tally_number};
  result = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &tally_number_tag, 
                                              tally_number_val, 1, matching_tally_number_sets );
    assert(MB_SUCCESS == result);
  std::cout << "number of matching meshsets=" << matching_tally_number_sets.size() << std::endl;
    assert( 1 == matching_tally_number_sets.size() );

  // identify which of the meshsets is existing
  MBEntityHandle existing_meshset;
  //if ( tally_meshset == matching_tally_number_sets.front() ) 
    //existing_meshset = matching_tally_number_sets.back();
  //else
    existing_meshset = matching_tally_number_sets.front();

  // get the existing elements from the set
  MBRange existing_elements;
  result = MBI->get_entities_by_type( existing_meshset, MBHEX, existing_elements );
    assert(MB_SUCCESS == result);

  // get this tally's elements
  //MBRange new_elements;
  //result = MBI->get_entities_by_type( tally_meshset, MBHEX, new_elements );
  //  assert(MB_SUCCESS == result);

  // check to make sure they have the same number of elements
  //assert( existing_elements.size() == new_elements.size() );
  assert( existing_elements.size() == n_elements );

  // get the nps of the existing and new tally
  //result = MBI->tag_create("NPS_TAG", sizeof(unsigned long long int), MB_TAG_SPARSE, 
  //                         MB_TYPE_OPAQUE, nps_tag, 0);
  //  assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  //unsigned long long int nps0, nps1;
  unsigned long long int nps0;
  MBRange sets_with_this_tag;
  result = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &nps_tag, 0, 1, sets_with_this_tag);
    assert(MB_SUCCESS == result);
    std::cout << "number of nps sets=" << sets_with_this_tag.size() << std::endl;
    //assert( 2 == sets_with_this_tag.size() );
    assert( 1 == sets_with_this_tag.size() );
  result = MBI->tag_get_data( nps_tag, &sets_with_this_tag.front(), 1, &nps0);
    assert(MB_SUCCESS == result);
  //result = MBI->tag_get_data( nps_tag, &sets_with_this_tag.back(), 1, &nps1);
  //  assert(MB_SUCCESS == result);
  std::cout << "nps0=" << nps0 << " nps1=" << nps1 << std::endl;
  new_nps = nps0 + nps1;

  // get tally values from the existing elements
  //result = MBI->tag_create("TALLY_TAG", sizeof(double), MB_TAG_DENSE, 
  //                         MB_TYPE_DOUBLE, tally_tag, 0);
  //  assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  //result = MBI->tag_create("ERROR_TAG", sizeof(double), MB_TAG_DENSE, 
  //                         MB_TYPE_DOUBLE, error_tag, 0);
  //  assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);
  //double *values0, *values1, *errors0, *errors1;
  //double *values1, *errors1;
  //values0 = new double [existing_elements.size()];
  //errors0 = new double [existing_elements.size()];
  double values0[existing_elements.size()];
  double errors0[existing_elements.size()];
  result = MBI->tag_get_data( tally_tag, existing_elements, &values0 );
    assert(MB_SUCCESS == result);
  //result = MBI->tag_get_data( error_tag, existing_elements, errors0 ); this worked with erros0 = new double [size];
  //  assert(MB_SUCCESS == result);
  result = MBI->tag_get_data( error_tag, existing_elements, &errors0 );
    assert(MB_SUCCESS == result);

  // get the tally values from the new elements
  //values1 = new double [new_elements.size()];
  //errors1 = new double [new_elements.size()];
  //result = MBI->tag_get_data( tally_tag, new_elements, values1 );
  //  assert(MB_SUCCESS == result);
  //result = MBI->tag_get_data( error_tag, new_elements, errors1 );
  //  assert(MB_SUCCESS == result);

  // average the values and errors
  result = average_tally_values( nps0, nps1, values0, values1, 
                                 errors0, errors1, n_elements );
    assert(MB_SUCCESS == result);
  
  // set the averaged information back onto the existing elements
  result = MBI->tag_set_data( tally_tag, existing_elements, &values0 );
    assert(MB_SUCCESS == result);
  result = MBI->tag_set_data( error_tag, existing_elements, &errors0 );
    assert(MB_SUCCESS == result);
  
  // cleanup
  //delete[] values1;
  //delete[] errors1;

  return MB_SUCCESS;
}


MBErrorCode ReadMCNP5::transform_point_to_cartesian(double *in, double *out, 
                                                    coordinate_system coord_sys ) {
                  //coordinate_system coord_sys, double *rmat) {

      //double q[3];

      // Apply the rotation matrix
      //for (unsigned int i=0; i < 3; i++) {
      //  q[i] =  p[0] * rmat[4*i  ] + p[1] * rmat[4*i+1]
      //        + p[2] * rmat[4*i+2] +        rmat[4*i+3];
      //}


  // Transform coordinate system
  switch( coord_sys ) {
    case CARTESIAN :
      out[0] = in[0]; out[1] = in[1]; out[2] = in[2];  // x, y, z
      break;
    case CYLINDRICAL :
      //r[0] = sqrt( q[0]*q[0] + q[1]*q[1] );   // r
      //r[1] = q[2];                            // z
      //r[2] = c2pi * ( atan2( q[1], q[0] ) );  // theta (in rotations)
      //std::cout << "r=" << p[0] << " z=" << p[1] << " theta=" << p[2] << std::endl;
      out[0] = in[0]*cos( 2*pi*in[2] ); // x
      out[1] = in[0]*sin( 2*pi*in[2] ); // y
      out[2] = in[1];                  // z
      break;
    case SPHERICAL :
      return MB_NOT_IMPLEMENTED;
      break;
    default :
      return MB_NOT_IMPLEMENTED;
      break;
  }
  
  return MB_SUCCESS;
}

// Average two tally values and their error. Return average values in the
// place of first tally values.
MBErrorCode ReadMCNP5::average_tally_values(const unsigned long long int nps0, 
                                            const unsigned long long int nps1,
					    double *values0,
					    const double *values1,
					    double *errors0,
					    const double *errors1,
					    const unsigned long int n_values) {
  
  for(unsigned long int i=0; i<n_values; i++) {
    //std::cout << " values0=" << values0[i] << " values1=" << values1[i] 
    //          << " errors0=" << errors0[i] << " errors1=" << errors1[i] << " nps0=" << nps0 << " nps1=" << nps1 << std::endl;
    errors0[i] = sqrt( pow(values0[i]*errors0[i]*nps0,2) + 
                       pow(values1[i]*errors1[i]*nps1,2) ) / 
		 (values0[i]*nps0 + values1[i]*nps1);

    values0[i] = ( values0[i]*nps0 + values1[i]*nps1 ) / (nps0 + nps1);

    //std::cout << " values0=" << values0[i] << " errors0=" << errors0[i] << std::endl;
  }
  // REMEMBER TO UPDATE NPS0 = NPS0 + NPS1 after this
  return MB_SUCCESS;
}
