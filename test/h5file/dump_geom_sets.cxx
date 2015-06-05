// This utility reads a MOAB HDF5 mesh file, and prints out
// the list of all geometry entitysets and the parent and
// child geometry entitysets of each such entityset.

#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <vector>
#include <map>
#include <string>

// Open HDF5 file
hid_t open_file( const char* name )
{
  hid_t handle = H5Fopen( name, H5F_ACC_RDONLY, H5P_DEFAULT );
  if (handle < 0) {
    fprintf(stderr, "Cannot open file: \"%s\"\n", name );
    exit (1);
  }
  return handle;
}

// Read a scalar attribute value from an HDF5 file
long read_scalar_attrib( hid_t file, const char* table, const char* attrib )
{
  hid_t table_id = H5Dopen( file, table );
  if (table_id < 0) {
    fprintf(stderr, "Invalid file.  Data not found: \"%s\"\n", table );
    exit (1);
  }
  
  hid_t attr_id = H5Aopen_name( table_id, attrib );
  H5Dclose( table_id );
  if (attr_id < 0) {
    fprintf(stderr, "Invalid file.  No \"%s\" attrib on \"%s\"\n", attrib, table );
    exit (1);
  }
  
  long value;
  herr_t rval = H5Aread( attr_id, H5T_NATIVE_LONG, &value );
  H5Aclose( attr_id );
  if (rval < 0) {
    fprintf(stderr, "Failed to read \"%s\" attrib on \"%s\"\n", attrib, table );
    exit (2);
  }
  
  return value;
}


// Read a data table from an HDF5 file
void read_data( hid_t file, const char* name, int expected_cols, std::vector<long>& data )
{
  hid_t handle = H5Dopen( file, name );
  if (handle < 0) {
    fprintf(stderr, "Invalid file.  Data not found: \"%s\"\n", name );
    exit (1);
  }

  hid_t space = H5Dget_space( handle );
  if (space < 0) {
    fprintf(stderr, "Internal error accessing: \"%s\"\n", name );
    exit (2);
  }
  
  int ndims = H5Sget_simple_extent_ndims( space );
  if (ndims < 0) {
    fprintf(stderr, "Internal error accessing: \"%s\"\n", name );
    exit (2);
  }
  else if (ndims < 1 || ndims > 2) {
    fprintf(stderr, "\"%s\" is an %d-dimension table.  Corrupt file?", name, ndims );
    exit (2);
  }
  
  hsize_t dims[2] = { 0, 1 };
  H5Sget_simple_extent_dims( space, dims, 0 );
  H5Sclose( space );
  
  if (dims[1] != expected_cols) {
    fprintf(stderr, "Error reading \"%s\": expected %d cols, has %d\n", name, expected_cols, (int)dims[1] );
    exit (2);
  }
  
  data.resize( dims[0] * dims[1] );
  herr_t rval = H5Dread( handle, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0] );
  if (rval < 0) {
    fprintf(stderr, "Error reading data from: \"%s\"\n", name );
    exit (1);
  }
  
  H5Dclose( handle );
}

// Given two vectors of the same length, create a map where the
// key is taken from the first vectgor and the value is taken 
// from the second vector.
void create_map( const std::vector<long>& ents, 
                 const std::vector<long>& vals,
                 std::map<long,long>& map,
                 const char* name )
{
  if (ents.size() != vals.size()) {
    fprintf(stderr, "Invalid data for tag \"%s\": mismatched table lengths.\n", name );
    exit (1);
  }
  
  std::vector<long>::const_iterator e_iter = ents.begin(), v_iter = vals.begin();
  for (; e_iter != ents.end(); ++e_iter, ++v_iter)
    map[*e_iter] = *v_iter;
} 

// Construct a string designator for a geometry entity given
// the handle of that geometry entity and maps from handle to
// dimension and handle to global id.
std::string ent_from_handle( const std::map<long,long>& dimmap,
                             const std::map<long,long>&  idmap,
                             long handle )
{
  std::string result;
  std::map<long,long>::const_iterator d_iter, i_iter;
  d_iter = dimmap.find( handle );
  i_iter =  idmap.find( handle );
  if (d_iter == dimmap.end() || i_iter == idmap.end()) 
    return result;

  switch (d_iter->second) {
    case 0: result += "v"; break;
    case 1: result += "c"; break;
    case 2: result += "s"; break;
    case 3: result += "V"; break;
    default:
      fprintf(stderr,"Invalid value in GEOM_DIMENSION tag data.\n");
      exit (1);
  }
  
  char buffer[128];
  sprintf(buffer, "%ld", i_iter->second );
  result += buffer;
  return result;
}
      
// Construct a string designator for list of geometry entities given
// a list of handles and maps from handle to
// dimension and handle to global id.
std::string ent_list_from_handles( const std::map<long,long>& dimmap,
                                   const std::map<long,long>&  idmap,
                                   const std::vector<long>& vect,
                                   long start, 
                                   long stop )
{
  std::string result;
  if (start >= (long)vect.size() || stop >= (long)vect.size() || stop < start) {
    fprintf(stderr, "Invalid set data.  Corrupt file?\n");
    exit(2);
  }
  std::vector<long>::const_iterator iter = vect.begin() + start+1,
                                     end = vect.begin() + stop+1;
  
  for (; iter != end; ++iter)
  {
    std::string tmp = ent_from_handle( dimmap, idmap, *iter );
    if (!tmp.empty()) {
      result += tmp;
      result += " ";
    }
    else
      result += "? ";
  }
  return result;
}

int main( int argc, char* argv[] )
{
    // Need a file name to read
  if (argc != 2) {
    printf("Usage: %s <file>\n", argv[0] );
    exit(1);
  }
  
    // Read everything we need from the file
  std::vector<long> set_meta, set_child, set_parent, dim_ents, dim_vals, id_ents, id_vals;
  hid_t file = open_file( argv[1] );
  read_data( file, "/tstt/sets/list", 4, set_meta );
  read_data( file, "/tstt/sets/parents", 1, set_parent );
  read_data( file, "/tstt/sets/children", 1, set_child );
  read_data( file, "/tstt/tags/GEOM_DIMENSION/id_list", 1, dim_ents );
  read_data( file, "/tstt/tags/GEOM_DIMENSION/values", 1, dim_vals );
  read_data( file, "/tstt/tags/GLOBAL_ID/id_list", 1, id_ents );
  read_data( file, "/tstt/tags/GLOBAL_ID/values", 1, id_vals );
  const long startid = read_scalar_attrib( file, "/tstt/sets/list", "start_id" );
  H5Fclose( file );
  
    // Construct handle->dimension and handle->global_id maps.
  std::map<long,long> dimmap, idmap;
  create_map( dim_ents, dim_vals, dimmap, "GEOM_DIMENSION" );
  create_map(  id_ents,  id_vals,  idmap, "GLOBAL_ID" ); 
  
    // For each entity set
  long parent_start = -1l, child_start = -1l;
  printf("handle  ID     Children             Parents\n");
  for (unsigned i = 0; i < set_meta.size(); i += 4)
  {
      // Get name
    long handle = startid + i/4;
    std::string name = ent_from_handle( dimmap, idmap, handle );
    if (name.empty()) // not a geometry set
      continue;
     
      // Get parents and children  
    long  child_end = set_meta[i+1];
    long parent_end = set_meta[i+2];
    std::string children = ent_list_from_handles( dimmap, idmap,  set_child,  child_start,  child_end );
    std::string  parents = ent_list_from_handles( dimmap, idmap, set_parent, parent_start, parent_end );
     child_start =  child_end;
    parent_start = parent_end;
    
      // Print
    printf( "%6ld  %-6s %-20s %-20s\n", handle, name.c_str(), children.c_str(), parents.c_str() );
  }
  
  return 0;
}
