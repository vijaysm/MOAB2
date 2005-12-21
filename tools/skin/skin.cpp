#include <iostream>
#include <time.h>
#include <vector>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "MBInterface.hpp"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"

void get_time_mem(double &tot_time, double &tot_mem);

// Different platforms follow different conventions for usage
#ifndef NT
#include <sys/resource.h>
#endif
#ifdef SOLARIS
extern "C" int getrusage(int, struct rusage *);
#ifndef RUSAGE_SELF
#include </usr/ucbinclude/sys/rusage.h>
#endif
#endif

const char FIXED_TAG[] = "fixed"; 

#define CHKERROR( A ) do { if (MB_SUCCESS != (A)) { \
 std::cerr << "Internal error at line " << __LINE__ << std::endl; \
 return 3; } } while(false)

int main( int argc, char* argv[] )
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] 
              << " [-b <block_num> [-b ...] ] [-p] [-s <sideset_num>] [-t] [-w]"
              << " <input_file> [<output_file>]" << std::endl;
    std::cerr << "Options: " << std::endl;
    std::cerr << "-a : Compute skin using vert-elem adjacencies (more memory, less time)." << std::endl;
    std::cerr << "-b <block_num> : Compute skin only for material set/block <block_num>." << std::endl;
    std::cerr << "-p : Print cpu & memory performance." << std::endl;
    std::cerr << "-s <sideset_num> : Put skin in neumann set/sideset <sideset_num>." << std::endl;
    std::cerr << "-t : Set 'FIXED' tag on skin vertices." << std::endl;
    std::cerr << "-w : Write out whole mesh (otherwise just writes skin)." << std::endl;
    
    return 1;
  }

  int i = 1;
  std::vector<int> matsets;
  int neuset_num = -1;
  bool write_tag = false, write_whole_mesh = false;
  bool print_perf = false;
  bool use_vert_elem_adjs = false;
  
  while (i < argc) {
    if (!strcmp(argv[i], "-a")) {
      i++;
      use_vert_elem_adjs = true;
    }
    else if (!strcmp(argv[i], "-b")) {
      i++;
      matsets.push_back(atoi(argv[i]));
      i++;
    }
    else if (!strcmp(argv[i], "-p")) {
      i++;
      print_perf = true;
    }
    else if (!strcmp(argv[i], "-s")) {
      i++;
      neuset_num = atoi(argv[i]);
      i++;
    }
    else if (!strcmp(argv[i], "-t")) {
      i++;
      write_tag = true;
    }
    
    else if (!strcmp(argv[i], "-w")) {
      i++;
      write_whole_mesh = true;
    }
    else {
      break;
    }
  }
  
  const char* input_file = argv[i++];
  const char* output_file = NULL;
  if (i < argc) 
    output_file = argv[i++];
  
  MBErrorCode result;
  MBCore mbimpl;
  MBInterface* iface = &mbimpl;
  
  if (print_perf) {
    double tmp_time1, tmp_mem1;
    get_time_mem(tmp_time1, tmp_mem1);
    std::cout << "Before reading: cpu time = " << tmp_time1 << ", memory = " 
              << tmp_mem1/1.0e6 << "MB." << std::endl;
  }

    // read input file
  result = iface->load_mesh( input_file );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to load \"" << input_file << "\"." << std::endl; 
    return 2;
  }
  std::cerr << "Read \"" << input_file << "\"" << std::endl;
  if (print_perf) {
    double tmp_time2, tmp_mem2;
    get_time_mem(tmp_time2, tmp_mem2);
    std::cout << "After reading: cpu time = " << tmp_time2 << ", memory = " 
              << tmp_mem2/1.0e6 << "MB." << std::endl;
  }
  
    // get entities of largest dimension
  int dim = 4;
  MBRange entities;
  while (entities.empty() && dim > 1)
  {
    dim--;
    result = iface->get_entities_by_dimension( 0, dim, entities );
    CHKERROR(result);
  }

  MBRange skin_ents;
  MBTag matset_tag = 0, neuset_tag = 0;
  result = iface->tag_get_handle(MATERIAL_SET_TAG_NAME, matset_tag);
  result = iface->tag_get_handle(NEUMANN_SET_TAG_NAME, neuset_tag);

  if (matsets.empty()) skin_ents = entities;
  else {
      // get all entities in the specified blocks
    if (0 == matset_tag) {
      std::cerr << "Couldn't find any material sets in this mesh." << std::endl;
      return 1;
    }
    
    for (std::vector<int>::iterator vit = matsets.begin(); vit != matsets.end(); vit++) {
      int this_matset = *vit;
      const void *this_matset_ptr = &this_matset;
      MBRange this_range, ent_range;
      result = iface->get_entities_by_type_and_tag(0, MBENTITYSET, &matset_tag,
                                                    &this_matset_ptr, 1, this_range);
      if (MB_SUCCESS != result) {
        std::cerr << "Trouble getting material set #" << *vit << std::endl;
        return 1;
      }
      else if (this_range.empty()) {
        std::cerr << "Warning: couldn't find material set " << *vit << std::endl;
        continue;
      }
      
      result = iface->get_entities_by_dimension(*this_range.begin(), dim, ent_range, true);
      if (MB_SUCCESS != result) continue;
      skin_ents.merge(ent_range);
    }
  }
  
  if (skin_ents.empty()) {
    std::cerr << "No entities for which to compute skin; exiting." << std::endl;
    return 1;
  }

  if (use_vert_elem_adjs) {
      // make a call which we know will generate vert-elem adjs
    MBRange dum_range;
    result = iface->get_adjacencies(&(*skin_ents.begin()), 1, 1, false,
                                    dum_range);
  }
  
  double tmp_time = 0.0, tmp_mem = 0.0;
  if (print_perf) {
    get_time_mem(tmp_time, tmp_mem);
    std::cout << "Before skinning: cpu time = " << tmp_time << ", memory = " 
              << tmp_mem/1.0e6 << "MB." << std::endl;
  }

    // skin the mesh
  MBRange forward_lower, reverse_lower;
  MBSkinner tool( iface );
  result = tool.find_skin( skin_ents, forward_lower, reverse_lower );
  MBRange boundary;
  boundary.merge( forward_lower );
  boundary.merge( reverse_lower );
  if (MB_SUCCESS != result || boundary.empty())
  {
    std::cerr << "Mesh skinning failed." << std::endl;
    return 3;
  }

  if (write_tag) {
      // get tag handle
    MBTag tag;
    result = iface->tag_get_handle( FIXED_TAG, tag );
    if (result == MB_SUCCESS)
    {
      int size;
      MBDataType type;
      iface->tag_get_size(tag, size);
      iface->tag_get_data_type(tag, type);
    
      if (size != sizeof(int) || type != MB_TYPE_INTEGER)
      {
        std::cerr << '"' << FIXED_TAG << "\" tag defined with incorrect size or type" << std::endl;
        return 3;
      }
    }
    else if (result == MB_TAG_NOT_FOUND)
    {
      int zero = 0;
      result = iface->tag_create( FIXED_TAG, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, tag, &zero );
      CHKERROR(result);
    }
    else
    {
      CHKERROR(result);
    }
  
      // Set tags
    std::vector<int> ones;
    MBRange bverts;
    result = iface->get_adjacencies(boundary, 0, false, bverts, MBInterface::UNION);
    if (MB_SUCCESS != result) {
      std::cerr << "Trouble getting vertices on boundary." << std::endl;
      return 1;
    }
    ones.resize( bverts.size(), 1 );
    result = iface->tag_set_data( tag, bverts, &ones[0] );
    CHKERROR(result);
  }
  
  if (-1 != neuset_num) {
      // create a neumann set with these entities
    if (0 == neuset_tag) {
      result = iface->tag_create("NEUMANN_SET_TAG_NAME", sizeof(int), MB_TAG_SPARSE,
                                  MB_TYPE_INTEGER, neuset_tag, NULL);
      if (MB_SUCCESS != result || 0 == neuset_tag) return 1;
    }
    

      // always create a forward neumann set, assuming we have something in the set
    MBEntityHandle forward_neuset = 0;
    result = iface->create_meshset(MESHSET_SET, forward_neuset);
    if (MB_SUCCESS != result || 0 == forward_neuset) return 1;
    result = iface->tag_set_data(neuset_tag, &forward_neuset, 1, &neuset_num);
    if (MB_SUCCESS != result) return 1;

    if (!forward_lower.empty()) {
      result = iface->add_entities(forward_neuset, forward_lower);
      if (MB_SUCCESS != result) return 1;
    }
    if (!reverse_lower.empty()) {
      MBEntityHandle reverse_neuset = 1;
      result = iface->create_meshset(MESHSET_SET, reverse_neuset);
      if (MB_SUCCESS != result || 0 == forward_neuset) return 1;

      result = iface->add_entities(reverse_neuset, reverse_lower);
      if (MB_SUCCESS != result) return 1;
      MBTag sense_tag;
      result = iface->tag_get_handle("SENSE", sense_tag);
      if (result == MB_TAG_NOT_FOUND) {
        int dum_sense = 0;
        result = iface->tag_create("SENSE", sizeof(int), MB_TAG_SPARSE, sense_tag, &dum_sense);
      }
      if (result != MB_SUCCESS) return 1;
      int sense_val = -1;
      result = iface->tag_set_data(neuset_tag, &reverse_neuset, 1, &sense_val);
      if (MB_SUCCESS != result) return 0;
      result = iface->add_entities(forward_neuset, &reverse_neuset, 1);
      if (MB_SUCCESS != result) return 0;
    }
  }

  if (NULL != output_file && write_whole_mesh) {
    
      // write output file
    result = iface->write_mesh( output_file);
    if (MB_SUCCESS != result)
    { 
      std::cerr << "Failed to write \"" << output_file << "\"." << std::endl; 
      return 2;
    }
    std::cerr << "Wrote \"" << output_file << "\"" << std::endl;
  }
  else if (NULL != output_file) {
      // write only skin; write them as one set
    MBEntityHandle skin_set;
    result = iface->create_meshset(MESHSET_SET, skin_set);
    if (MB_SUCCESS != result) return 1;
    result = iface->add_entities(skin_set, forward_lower);
    if (MB_SUCCESS != result) return 1;
    result = iface->add_entities(skin_set, reverse_lower);
    if (MB_SUCCESS != result) return 1;

    MBRange this_range, ent_range;
    result = iface->get_entities_by_type_and_tag(0, MBENTITYSET, &matset_tag,
                                                  NULL, 0, this_range);
    if (!this_range.empty()) iface->delete_entities(this_range);

    int dum = 10000;
    result = iface->tag_set_data(matset_tag, &skin_set, 1, &dum);
    

    result = iface->write_mesh( output_file, &skin_set, 1);
    if (MB_SUCCESS != result)
    { 
      std::cerr << "Failed to write \"" << output_file << "\"." << std::endl; 
      return 2;
    }
    std::cerr << "Wrote \"" << output_file << "\"" << std::endl;
  }

  if (print_perf) {
    double tot_time, tot_mem;
    get_time_mem(tot_time, tot_mem);
    std::cout << "Total cpu time = " << tot_time << " seconds." << std::endl;
    std::cout << "Total skin cpu time = " << tot_time-tmp_time << " seconds." << std::endl;
    std::cout << "Total memory = " << tot_mem/1.0e6 << " MB." << std::endl;
    std::cout << "Total skin memory = " << (tot_mem-tmp_mem)/1.0e6 << " MB." << std::endl;
    std::cout << "Entities: " << std::endl;
    iface->list_entities(0, 0);
  }
  
  return 0;
}

void get_time_mem(double &tot_time, double &tot_mem) 
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  double utime = (double)r_usage.ru_utime.tv_sec +
    ((double)r_usage.ru_utime.tv_usec/1.e6);
  double stime = (double)r_usage.ru_stime.tv_sec +
    ((double)r_usage.ru_stime.tv_usec/1.e6);
  tot_time = utime + stime;
  tot_mem = 0;
  if (0 != r_usage.ru_maxrss) {
    tot_mem = r_usage.ru_idrss; 
  }
  else {
      // this machine doesn't return rss - try going to /proc
      // print the file name to open
    char file_str[4096], dum_str[4096];
    int file_ptr = -1, file_len;
    file_ptr = open("/proc/self/stat", O_RDONLY);
    file_len = read(file_ptr, file_str, sizeof(file_str)-1);
    if (file_len == 0) return;
    
    close(file_ptr);
    file_str[file_len] = '\0';
      // read the preceeding fields and the ones we really want...
    int dum_int;
    unsigned int dum_uint, vm_size, rss;
    int num_fields = sscanf(file_str, 
                            "%d " // pid
                            "%s " // comm
                            "%c " // state
                            "%d %d %d %d %d " // ppid, pgrp, session, tty, tpgid
                            "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
                            "%d %d %d %d %d %d " // utime, stime, cutime, cstime, counter, priority
                            "%u %u " // timeout, itrealvalue
                            "%d " // starttime
                            "%u %u", // vsize, rss
                            &dum_int, 
                            dum_str, 
                            dum_str, 
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, 
                            &dum_int,
                            &vm_size, &rss);
    if (num_fields == 24)
      tot_mem = ((double)vm_size);
  }
}

  
  
  

