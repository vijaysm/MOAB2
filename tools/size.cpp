/*
 * Copyright (c) 2005 Lawrence Livermore National Laboratory under
 * contract number B545069 with the University of Wisconsin - Madison.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <stdio.h>
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#  include <unistd.h>
#  include <termios.h>
#  include <sys/ioctl.h>
#endif
#include <math.h>
#include <assert.h>
#include <float.h>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"
#include "moab/Interface.hpp"

using namespace moab;

#include "measure.hpp"

void usage( const char* exe )
{
  std::cerr << "Usage: " << exe << " [-g] [-m] [-l] [-ll] <mesh file list>" << std::endl
            << "-g  : print counts by geometric owner" << std::endl
            << "-m  : print counts per block/boundary" << std::endl
            << "-l  : print counts of mesh" << std::endl
            << "-ll : verbose listing of every entity" << std::endl
            ;
  exit(1);
}


Core mb;


struct stat_set 
{
  double sum;
  double sqr;
  double min; 
  double max;
  long count;
  
  inline stat_set() : sum(0), sqr(0), min(HUGE_VAL), max(0), count (0) {}
  
  inline void add( double val )
  {
    if (val < min)
      min = val;
    if (val > max)
      max = val;
    sum += val;
    sqr += val*val;
    ++count;
  }
  
  inline void add( const stat_set& stats )
  {
    if (stats.min < min)
      min = stats.min;
    if (stats.max > max)
      max = stats.max;
    sum += stats.sum;
    sqr += stats.sqr;
    count += stats.count;
  }
  
  inline void clear()
  {
    sum = sqr = max = count = 0;
    min = HUGE_VAL;
  }
};

struct set_stats {
  stat_set stats[MBMAXTYPE];
  stat_set edge_uses;
  size_t nodes;
  
  void add( const set_stats& other )
  {
    for (int i = 0; i < MBMAXTYPE; ++i)
      stats[i].add( other.stats[i] );
    edge_uses.add( other.edge_uses );
  }
  
  void clear()
  {
    for (int i = 0; i < MBMAXTYPE; ++i)
      stats[i].clear();
    edge_uses.clear();
  }
    
};


ErrorCode gather_set_stats( EntityHandle set, set_stats& stats )
{
  ErrorCode rval = MB_SUCCESS;
  
  int count;
  rval = mb.get_number_entities_by_type( set, MBVERTEX, count );
  if (MB_SUCCESS != rval) return rval;
  stats.nodes = count;
  
  int edge_vtx_idx[2];
  std::vector<EntityHandle> conn;
  std::vector<double> coords;
  for (EntityType type = MBEDGE; type < MBENTITYSET; ++type)
  {
    int num_edges = CN::NumSubEntities( type, 1 );
    
    Range range;
    rval = mb.get_entities_by_type( set, type, range, true );
    if (MB_SUCCESS != rval) return rval;
    for (Range::iterator i = range.begin(); i != range.end(); ++i)
    {
      rval = mb.get_connectivity( &*i, 1, conn, true );
      if (MB_SUCCESS != rval) return rval;
      coords.resize( 3*conn.size() );
      rval = mb.get_coords( &conn[0], conn.size(), &coords[0] );
      if (MB_SUCCESS != rval) return rval;
      stats.stats[type].add( measure( type, conn.size(), &coords[0] ) );
      
      if (type != MBEDGE)
      {
        if (type == MBPOLYGON)
          num_edges = conn.size();

        for (int e = 0; e < num_edges; ++e)
        {
          if (type == MBPOLYGON) {
            edge_vtx_idx[0] = e;
            edge_vtx_idx[1] = e+1;
          }
          else
            CN::SubEntityVertexIndices( type, 1, e, edge_vtx_idx );
          stats.edge_uses.add( edge_length( &coords[3*edge_vtx_idx[0]],
                                            &coords[3*edge_vtx_idx[1]] ) );
        }
      }
    }
  }
  return MB_SUCCESS;
}

const char* dashes( unsigned count )
{
  static std::vector<char> dashes;
  dashes.clear();
  dashes.resize( count + 1, '-' );
  dashes[count] = '\0';
  return &dashes[0];
}

void print_stats( set_stats& stats )
{
  const char* edge_use_name = "1D Side";
  const char* vertex_name = "Vertex";
  
  bool have_some = stats.edge_uses.count > 0;
  for (int i = 0; i < MBMAXTYPE; ++i)
    if (stats.stats[i].count > 0)
      have_some = true;
  
  if (!have_some)
  {
    std::cout << "NO MESH" << std::endl;
    return;
  }
  
    // get field widths
  unsigned type_width = std::max( strlen(vertex_name), strlen( edge_use_name ) );
  unsigned count_width = 5;
  unsigned total_width = 5;
  unsigned total_prec = 2;
  unsigned precision = 5;
  int total_log = -10000;
  
  unsigned node_count_width = (unsigned)(ceil(log10(stats.nodes))) + 1;
  if (count_width < node_count_width)
    count_width = node_count_width;
  
  for (EntityType i = MBEDGE; i < MBMAXTYPE; ++i)
  {
    stat_set& s = (i == MBMAXTYPE) ? stats.edge_uses : stats.stats[i];

    if (s.count == 0)
      continue;
    
    unsigned len = strlen(CN::EntityTypeName(i));
    if (len > type_width)
      type_width = len;
    
    unsigned cw = (unsigned)(ceil(log10(s.count))) + 1;
    if (cw > count_width)
      count_width = cw;
    
    int tl = (unsigned)(ceil(log10(fabs(s.sum)))) + 1;
    if (tl > total_log)
      total_log = tl;
  }
  
  if (total_log > (int)total_width)
  {
    total_width = 8;
    total_prec = 2;
  }
  else if (total_log <= -(int)total_width)
  {
    total_width = -total_log + 5;
    total_prec = 2;
  }
  else if (total_log < 1)
  {
    total_width = -total_log + 4;
    total_prec = -total_log + 1;
  }
  else
  {
    total_width += 2;
  }
    
  
  // get terminal width
  unsigned term_width = 80;
#if !defined(_MSC_VER) && !defined(__MINGW32__)
  struct winsize size;
  if ( ioctl( fileno(stdout), TIOCGWINSZ, (char*)&size ) == 0 )
    term_width = size.ws_col;
  if (!term_width) term_width = 80;
#endif
  assert(term_width > 7 + type_width + count_width + total_width);
  
  term_width -= 7; // spaces
  term_width -= type_width;
  term_width -= count_width;
  term_width -= total_width;
  unsigned val_width = term_width / 5;
  if (val_width < 8)
    val_width = 8;
  
  printf( "%*s %*s %*s %*s %*s %*s %*s %*s\n",
          type_width, "type",
          count_width, "count",
          total_width, "total",
          val_width, "minimum",
          val_width, "average",
          val_width, "rms",
          val_width, "maximum",
          val_width, "std.dev." );
  
  printf( "%*s ", type_width, dashes(type_width) );
  printf( "%*s ", count_width, dashes(count_width) );
  printf( "%*s ", total_width, dashes(total_width) );
  printf( "%*s ", val_width, dashes(val_width) );
  printf( "%*s ", val_width, dashes(val_width) );
  printf( "%*s ", val_width, dashes(val_width) );
  printf( "%*s ", val_width, dashes(val_width) );
  printf( "%*s\n", val_width, dashes(val_width) );
  
  for (EntityType i = MBEDGE; i <= MBMAXTYPE; ++i)
  {
    stat_set& s = (i == MBMAXTYPE) ? stats.edge_uses : stats.stats[i];
    
    if (s.count == 0)
      continue;

    double tmp_dbl = s.sqr / s.count - s.sum*s.sum / (double)s.count / (double)s.count;
    if (tmp_dbl < 0.0) {
      if (tmp_dbl < -100.0*DBL_EPSILON)
        std::cout << "WARNING: stat values dubious, s^2 - sig_s = " << tmp_dbl << std::endl;
      tmp_dbl = 0.0;
    }
    
    printf( "%*s %*ld %*.*g %*.*g %*.*g %*.*g %*.*g %*.*g\n",
            type_width, i == MBMAXTYPE ? edge_use_name : CN::EntityTypeName(i),
            count_width, s.count,
            total_width, total_prec, s.sum,
            val_width, precision, s.min,
            val_width, precision, s.sum / s.count,
            val_width, precision, sqrt( s.sqr / s.count ),
            val_width, precision, s.max,
            val_width, precision,  
            sqrt(tmp_dbl)
          );
  }
  printf( "%*s %*lu\n", type_width, vertex_name, count_width, (unsigned long)stats.nodes );
  
  puts("");
}

const char* geom_type_names[] = { "Vertex", "Curve", "Surface", "Volume" } ;
const char* mesh_type_names[] = { "Nodeset", "Sideset", "Block" };
const char* mesh_type_tags[] = { DIRICHLET_SET_TAG_NAME, NEUMANN_SET_TAG_NAME, MATERIAL_SET_TAG_NAME };

int main( int argc, char* argv[] )
{
  bool geom_owners = false;
  bool mesh_owners = false;
  bool just_list = false;
  bool just_list_basic = false;
  std::vector<std::string> file_list;
  set_stats total_stats, file_stats;
  
  for (int i = 1; i < argc; ++i)
  {
    if (!strcmp(argv[i],"-g"))
      geom_owners = true;
    else if (!strcmp(argv[i],"-ll"))
      just_list = true;
    else if (!strcmp(argv[i],"-l"))
      just_list_basic = true;
    else if (!strcmp(argv[i],"-m"))
      mesh_owners = true;
    else if (*argv[i] && *argv[i] != '-')
      file_list.push_back( argv[i] );
    else
    {
      std::cerr << "Invalid option: \"" << argv[i] << '"' << std::endl;
      usage(argv[0]);
    }
  }
  
  if (file_list.empty())
    usage(argv[0]);
  
  for (std::vector<std::string>::iterator f = file_list.begin(); 
       f != file_list.end(); ++f)
  {
    printf("File %s:\n", f->c_str() );
    if (MB_SUCCESS != mb.load_mesh( f->c_str(), 0, 0 ))
    {
      fprintf(stderr, "Error reading file: %s\n", f->c_str() );
      return 1;
    }

    if (MB_SUCCESS != gather_set_stats( 0, file_stats ))
    {
      fprintf(stderr, "Error processing mesh from file: %s\n", f->c_str());
      return 1;
    }
    
    total_stats.add( file_stats );
    print_stats( file_stats );
    file_stats.clear();
    
    if (geom_owners)
    {
      Range entities;
      Tag dim_tag = 0, id_tag = 0;
      ErrorCode rval = mb.tag_get_handle( GEOM_DIMENSION_TAG_NAME, dim_tag );
      if (MB_TAG_NOT_FOUND == rval) 
      {
        fprintf( stderr, "No geometry tag defined.\n" );
      }
      else if (MB_SUCCESS != rval)
      {
        fprintf( stderr, "Error retreiving geometry tag.\n");
        return 2;
      }
      
      rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, id_tag );
      if (MB_TAG_NOT_FOUND == rval) 
      {
        fprintf( stderr, "No ID tag defined.\n" );
      }
      else if (MB_SUCCESS != rval)
      {
        fprintf( stderr, "Error retreiving ID tag.\n");
        return 2;
      }
      
      if (dim_tag && id_tag)
      {
        if (MB_SUCCESS != mb.get_entities_by_type_and_tag( 0, 
                                                        MBENTITYSET, 
                                                        &dim_tag,
                                                        0,
                                                        1,
                                                        entities ))
        {
          fprintf( stderr, "Error retreiving geometry entitities.\n" );
        }
      }
      
      if (entities.empty())
      {
        fprintf( stderr, "No geometry entities defined in file.\n" );
      }
      
      for (Range::iterator i = entities.begin(); i != entities.end(); ++i)
      {
        int id = 0, dim = 0;
        if (MB_SUCCESS != mb.tag_get_data( dim_tag, &*i, 1, &dim ) ||
            MB_SUCCESS != mb.tag_get_data(  id_tag, &*i, 1,  &id ))
        {
          fprintf( stderr, "Error retreiving tag data for geometry entity.\n");
          continue;
        }
        
          // skip vertices
        if (0 == dim)
          continue;
        
        printf( "%s %d:\n", geom_type_names[dim], id );
        if (MB_SUCCESS != gather_set_stats( *i, file_stats ))
          fprintf(stderr, "Error processing mesh from file: %s\n", f->c_str());
        else
          print_stats( file_stats );
        file_stats.clear();
      }
    }


    if (mesh_owners)
    {
      for (int t = 0; t < 3; ++t)
      {
        Range entities;
        Tag tag = 0;
        ErrorCode rval = mb.tag_get_handle( mesh_type_tags[t], tag );
        if (MB_TAG_NOT_FOUND == rval) 
        {
          continue;
        }
        else if (MB_SUCCESS != rval)
        {
          fprintf( stderr, "Error retreiving %s tag.\n", mesh_type_tags[t]);
          return 2;
        }
      
        if (MB_SUCCESS != mb.get_entities_by_type_and_tag( 0, 
                                                        MBENTITYSET, 
                                                        &tag,
                                                        0,
                                                        1,
                                                        entities ))
        {
          fprintf( stderr, "Error retreiving %s entitities.\n", mesh_type_names[t] );
          continue;
        }
      
        for (Range::iterator i = entities.begin(); i != entities.end(); ++i)
        {
          int id = 0;
          if (MB_SUCCESS != mb.tag_get_data( tag, &*i, 1, &id ))
          {
            fprintf( stderr, "Error retreiving tag data for %s entity.\n", mesh_type_names[t]);
            continue;
          }

          printf( "%s %d:\n", mesh_type_names[t], id );
          if (MB_SUCCESS != gather_set_stats( *i, file_stats ))
            fprintf(stderr, "Error processing mesh from file: %s\n", f->c_str());
          else
            print_stats( file_stats );
          file_stats.clear();
        }
      }
    }

    if (just_list) mb.list_entities(0, 1);
   
    if (just_list_basic) mb.list_entities(0, 0);
    
    mb.delete_mesh();
  }
  
  if (file_list.size() > 1)
  {
    printf("Total for all files:\n");
    print_stats( total_stats );
  }
  
  return 0;
}

  
