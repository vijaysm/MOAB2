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
#include <vector>
#include <cstdlib>

#include "parse.hpp"

#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBInterface.hpp"
#define IS_BUILDING_MB
#include "MBInternals.hpp"
#undef IS_BUILDING_MB

#define CALL(A,B) \
  do { MBErrorCode _r = iface->A B ; \
       if (MB_SUCCESS != _r) { \
         std::cerr << #A << #B << " failed at " << __FILE__ << ":" << __LINE__ << std::endl; \
         exit( 5 ); }\
       } \
  while (false)


MBInterface* iface = 0;
const char* exe_name = 0;

void usage( )
{
  std::cerr << "Usage: " << exe_name << " <options> <input_file> <output_file>" << std::endl
            << "Options: " << std::endl
            << "  -t <ident_tag>[=<value>]  " << std::endl
            << "  -d <data_tag>[=<default>] " << std::endl
            << "  -c <data_tag=type:size>[=defult] " << std::endl
            << "  -w <write_tag>            " << std::endl
            << "  -n|-e                     " << std::endl
            << std::endl
            << "At least one ident_tag must be specified.  This tag will be used to "
               "select entity sets to process.  If more than one ident_tag is "
               "specified, entity sets that match A\bAN\bNY\bY the specified tags "
               "will be selected (logical OR).  A tag value to match may optionally "
               "be specified for each ident_tag. " << std::endl << std::endl
            << "The data_tag is the tag on the selected entity sets for which the "
               "value should be propogated to the contained mesn entities.  If no "
               "data_tag is specified, the value of the ident_tag will be used.  If "
               "multiple ident_tags are specified, the data_tag must be specified "
               "explicitly.  If the data tag is specified with the '-d' option, the "
               "tag must exist and the optional value will be set on any entiites for "
               "which the owning set doesn't have the data_tag defined.  If the data_tag "
               "is specified with the '-c' option, the tag will be created." << std::endl << std::endl
            << "If the write_tag is specified, the tag propogated to the mesh "
               "entities contained in the set will be given this name.  If no "
               "write_tag is specified, the data_tag will be used." << std::endl << std::endl
            << "If -\b-n\bn is specified, the tag values will be propogated only "
               "to nodes contained in the identified sets.  If -\b-e\be is specified, "
               "the tags will be propogated to the contained elements.  If neither is "
               "specified, the default is both." << std::endl << std::endl
            << "The syntax for specifying tag values is as follows: " << std::endl;
  tag_syntax(std::cerr);
  std::cerr << std::endl;
  exit(1);
}

void about( )
{
  std::cerr << "A utility to propogate tag values from the entity sets "
               "containing mesh entities to the entities contained in "
               "those sets." << std::endl << std::endl;
  usage();
}
  
void parse_error( const char* msg, const char* val = 0 )
{
  std::cerr << msg;
  if (val)
    std::cerr << ": " << val;
  std::cerr << std::endl;
  std::cerr << "Try '" << exe_name << " -h' for help" << std::endl;
  exit(1);
}

int main( int argc, char* argv[] )
{
  MBCore mb_core;
  exe_name = argv[0];
  iface = &mb_core;

  if (argc == 1)
    about();
    
    // find file names 
    // load input file before processing other options so
    // tags are defined
  const char* input_name = 0;
  const char* output_name = 0;
  for (int i = 1; i < argc; ++i)
  {
    if (argv[i][0] == '-')
    {
      switch (argv[i][1]) { case 't': case 'c': case 'd': case 'w': ++i; }
      continue;
    }
    
    if (!input_name)
      input_name = argv[i];
    else if (!output_name)
      output_name = argv[i];
    else 
      parse_error( "Unexpected argument", argv[i] );
  }
  
  if (!input_name)
    parse_error( "No input file specified." );
  if (!output_name)
    parse_error( "No output file specified." );
  
    // Read the input file
  if (MB_SUCCESS != iface->load_mesh( input_name ))
  {
    std::cerr << "Failed to read file: " << input_name << std::endl;
    std::string message;
    if (MB_SUCCESS == iface->get_last_error(message))
      std::cerr << message << std::endl;
    return 2;
  }

    
    
  bool nodes_spec = false;
  bool elems_spec = false;
  bool have_data_tag = false;
  const char* write_tag_name = 0;
  MBTag write_tag = 0;
  TagSpec data_tag;
  typedef std::vector<TagSpec> TagVect;
  TagVect ident_tags;
  int data_size = 0;
  
  for (int i = 1; i < argc; ++i)
  {
    if (argv[i] == input_name || argv[i] == output_name)
      continue;
    else if (!strcmp(argv[i],"-h"))
      usage();
    else if (!strcmp(argv[i],"-n"))
      nodes_spec = true;
    else if (!strcmp(argv[i], "-e"))
      elems_spec = true;
    else if (!argv[i][0])
      usage();
    else
    {
      char flag = argv[i][1];
      if ((flag != 't' && flag != 'd' && flag != 'w' && flag != 'c') || argv[i][2])
        parse_error( "Invalid argument", argv[i] );
      
      ++i;
      if (i == argc)
        parse_error( "Expected tag spec following option", argv[i-1] );
    
      if (flag == 'w')
      {
        if (write_tag_name)
          parse_error( "Invalid argument", argv[i] );
        write_tag_name = argv[i];
      }
      else if (flag == 'c')
      {
        TagSpec spec;
        if (parse_tag_create( argv[i], spec, iface ))
          parse_error( "Failed to parse tag spec", argv[i] );
        
        if (have_data_tag)
          parse_error( "Invalid argument", argv[i] );
        
        data_tag = spec;
        have_data_tag = true;
      }         
      else
      {
        TagSpec spec;
        if (parse_tag_spec( argv[i], spec, iface))
          parse_error("Failed to parse tag spec", argv[i] );
        
        if (flag == 'd')
        {
          if (have_data_tag)
            parse_error( "Invalid argument", argv[i] );
         
          data_tag = spec;
          have_data_tag = true;
        }
        else
        {
          ident_tags.push_back( spec );
        }
      }
    }
  } // for(args)
  
    // if neither, default to both
  if (!nodes_spec && !elems_spec)
    nodes_spec = elems_spec = true;
  
    // must have at least one identifying tag
  if (ident_tags.empty())
    parse_error ("At least one identifying tag must be specified.");
  
    // If data tag wasn't specified, use identifying tag for data
  if (!have_data_tag)
  {
    if (ident_tags.size() > 1)
      parse_error ("No data tag specified.");
    data_tag.value = 0;
    data_tag.handle = ident_tags[0].handle;
  }
  CALL( tag_get_size, (data_tag.handle, data_size) );
  
    // If write dat wasn't specified, use data tag 
  if (!write_tag_name)
  {
    write_tag = data_tag.handle;
  }
    // If write tag was specified, if it exists its type
    // msut match that of the data tag.  If it doesn't exist,
    // create it.
  else
  {
    MBDataType data_type;
    CALL( tag_get_data_type, (data_tag.handle, data_type) );
    
    MBErrorCode rval = iface->tag_get_handle( write_tag_name, write_tag );
    if (MB_FAILURE == rval)
    {
      std::cerr << "Unknown error retrieving handle for tag: " << write_tag_name << std::endl;
      exit( 5 );
    }
    else if (MB_TAG_NOT_FOUND == rval)
    {
      CALL( tag_create, (write_tag_name, data_size, MB_TAG_SPARSE, data_type, write_tag, 0) );
    }
    else
    {
      MBDataType write_type;
      int write_size;
      CALL( tag_get_data_type, (write_tag, write_type) );
      CALL( tag_get_size, (write_tag, write_size) );
      
      if (data_size != write_size ||
          (write_type != MB_TYPE_OPAQUE &&
           data_type != MB_TYPE_OPAQUE &&
           write_type != data_type)
          )
      {
        std::cerr << "Data and write tags have incompatible types." << std::endl;
        exit( 3 );
      }
    }
  }
  
  
  /**************** Done processing input -- do actual work ****************/
 
    // Get list of sets with identifying tags
  MBRange sets, temp;
  for (TagVect::iterator i = ident_tags.begin(); i != ident_tags.end(); ++i)
  {
    const void* value[] = { i->value };
    CALL( get_entities_by_type_and_tag,
          (0, MBENTITYSET, &i->handle, i->value ? value : 0, 1, temp) );
    sets.merge( temp );
  }
  
    // For each set, set tag on contained entities
  std::vector<unsigned char> tag_data(data_size);
  for (MBRange::iterator i = sets.begin(); i != sets.end(); ++i)
  {
      // Get tag value
    MBErrorCode rval = iface->tag_get_data( data_tag.handle, &*i, 1, &tag_data[0] );
    if (MB_TAG_NOT_FOUND == rval)
    {
      if (!data_tag.value)
      {
        std::cerr << "Data tag not set for entityset " 
                  << iface->id_from_handle(*i) << std::endl;
        continue;
      }
      memcpy( &tag_data[0], data_tag.value, data_size );
    }
    else if (MB_SUCCESS != rval)
    {
      CALL( tag_get_data, (data_tag.handle, &*i, 1, &tag_data[0]) );
    }
    
      // Get entities
    MBRange entities;
    CALL( get_entities_by_handle, (*i, entities, true) );
    int junk;
    MBRange::iterator eb = entities.lower_bound( entities.begin(), 
                                                 entities.end(),
                                                 CREATE_HANDLE(MBEDGE,0,junk) );
    MBRange::iterator j   = nodes_spec ? entities.begin() : eb;
    MBRange::iterator end = elems_spec ? entities.end()   : eb;
    for ( ; j != end; ++j)
      CALL( tag_set_data, (write_tag, &*j, 1, &tag_data[0]) );
  }
  
    // Write the output file
  if (MB_SUCCESS != iface->write_mesh( output_name ))
  {
    std::cerr << "Failed to write file: " << output_name << std::endl;
    std::string message;
    if (MB_SUCCESS == iface->get_last_error(message))
      std::cerr << message << std::endl;
    return 2;
  }
  
  return 0;
}

  
  
    
    
      
  
  
  
