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

#include "Tqdcfr.hpp"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBReadUtilIface.hpp"
#include "GeomTopoTool.hpp"
#include "MBTagConventions.hpp"
#include <assert.h>

const bool debug = false;
const int ACIS_DIMS[] = {-1, 3, -1, 2, -1, -1, 1, 0, -1, -1};
const char default_acis_dump_file[] = "dumped_acis.sat";
const char acis_dump_file_tag_name[] = "__ACISDumpFile";
const char Tqdcfr::geom_categories[][CATEGORY_TAG_NAME_LENGTH] = 
{"Vertex\0", "Curve\0", "Surface\0", "Volume\0"};
  
// acis dimensions for each entity type, to match
// enum {BODY, LUMP, SHELL, FACE, LOOP, COEDGE, EDGE, VERTEX, ATTRIB, UNKNOWN} 

MBReaderIface* Tqdcfr::factory( MBInterface* iface )
  { return new Tqdcfr( iface ); }

Tqdcfr::Tqdcfr(MBInterface *impl) 
    : cubFile(NULL), globalIdTag(0), geomTag(0), uniqueIdTag(0), groupTag(0), 
      blockTag(0), nsTag(0), ssTag(0), attribVectorTag(0), entityNameTag(0),
      categoryTag(0), instance(this)
{
  assert(NULL != impl);
  mdbImpl = impl;
  std::string iface_name = "MBReadUtilIface";
  impl->query_interface(iface_name, reinterpret_cast<void**>(&readUtilIface));
  assert(NULL != readUtilIface);

  currNodeIdOffset = -1;
  for (MBEntityType this_type = MBVERTEX; this_type < MBMAXTYPE; this_type++)
    currElementIdOffset[this_type] = -1;
}

Tqdcfr::~Tqdcfr() 
{
  std::string iface_name = "MBReadUtilIface";
  mdbImpl->release_interface(iface_name, readUtilIface);
}

  
MBErrorCode Tqdcfr::load_file(const char *file_name,
                              const int*, const int) 
{
  Tqdcfr* instance = this;
  
    // open file
  cubFile = fopen(file_name, "rb");
  if (NULL == cubFile) return MB_FAILURE;
  
    // verify magic string
  FREADC(4);
  if (!(char_buf[0] == 'C' && char_buf[1] == 'U' && 
         char_buf[2] == 'B' && char_buf[3] == 'E')) 
    return MB_FAILURE;

    // ***********************
    // read model header type information...
    // ***********************
  read_file_header();

  read_model_entries();
  
    // read model metadata
  read_meta_data(fileTOC.modelMetaDataOffset, modelMetaData);

    // ***********************
    // read mesh...
    // ***********************
  int index = find_model(mesh); assert(-1 != index);
  ModelEntry *mesh_model = &modelEntries[index];
  
    // first the header & metadata info
  mesh_model->read_header_info(this);
  mesh_model->read_metadata_info(this);

    // now read in mesh for each geometry entity
  for (int gindex = 0; 
       gindex < mesh_model->feModelHeader.geomArray.numEntities;
       gindex++) {
    EntityHeader *geom_header = &mesh_model->feGeomH[gindex];

      // read nodes
    read_nodes(mesh_model, geom_header);
    
      // read elements
    read_elements(mesh_model, geom_header);
  }

    // ***********************
    // read acis records...
    // ***********************
  read_acis_records();

    // ***********************
    // read groups...
    // ***********************
  for (int grindex = 0; 
       grindex < mesh_model->feModelHeader.groupArray.numEntities;
       grindex++) {
    EntityHeader *group_header = &mesh_model->feGroupH[grindex];
    read_group(mesh_model, group_header);
  }
  
    // ***********************
    // read blocks...
    // ***********************
  for (int blindex = 0; 
       blindex < mesh_model->feModelHeader.blockArray.numEntities;
       blindex++) {
    EntityHeader *block_header = &mesh_model->feBlockH[blindex];
    read_block(mesh_model, block_header);
  }

    // ***********************
    // read nodesets...
    // ***********************
  for (int nsindex = 0; 
       nsindex < mesh_model->feModelHeader.nodesetArray.numEntities;
       nsindex++) {
    EntityHeader *nodeset_header = &mesh_model->feNodeSetH[nsindex];
    read_nodeset(mesh_model, nodeset_header);
  }

    // ***********************
    // read sidesets...
    // ***********************
  for (int ssindex = 0; 
       ssindex < mesh_model->feModelHeader.sidesetArray.numEntities;
       ssindex++) {
    EntityHeader *sideset_header = &mesh_model->feSideSetH[ssindex];
    read_sideset(mesh_model, sideset_header);
  }

  if (debug) {
    std::cout << "Read the following mesh:" << std::endl;
    std::string dum;
    mdbImpl->list_entities(0, 0);
  }

    // **************************
    // restore geometric topology
    // **************************
  GeomTopoTool gtt(mdbImpl);
  MBErrorCode result = gtt.restore_topology();
  
  return result;
}

void Tqdcfr::read_nodeset(Tqdcfr::ModelEntry *model,
                          Tqdcfr::EntityHeader *nodeseth)  
{
  if (nodeseth->numTypes == 0) return;

  Tqdcfr* instance = this;

    // position file
  FSEEK(model->modelOffset+nodeseth->nodeOffset);
  
    // read ids for each entity type
  int this_type, num_ents;
  for (int i = 0; i < nodeseth->numTypes; i++) {
      // get how many and what type
    FREADI(2);
    this_type = int_buf[0];
    num_ents = int_buf[1];

      // now get the ids
    FREADI(num_ents);
    
      // get the entity handles (sets or entities) corresponding to this type
      // and put into this ns's set
    add_set_entities(this_type, num_ents, &(int_buf[0]), nodeseth->setHandle);
  }
}

void Tqdcfr::read_sideset(Tqdcfr::ModelEntry *model,
                          Tqdcfr::EntityHeader *sideseth)  
{
  if (sideseth->numTypes == 0) return;

    // position file
  FSEEK(model->modelOffset+sideseth->nodeOffset);
  
    // read ids for each entity type
  int this_type, num_ents, sense_size;
  for (int i = 0; i < sideseth->numTypes; i++) {
      // get how many and what type
    FREADI(3);
    this_type = int_buf[0];
    num_ents = int_buf[1];
    sense_size = int_buf[2];

      // now get the ids
    FREADI(num_ents);
    
      // get the entity handles (sets or entities) corresponding to this type
      // and put into this ns's set
    add_set_entities(this_type, num_ents, &(int_buf[0]), sideseth->setHandle);

    if (sense_size == 1) {
        // byte-size sense flags; make sure read ends aligned...
      int read_length = (num_ents / 8) * 8;
      if (read_length < num_ents) read_length += 8;
      FREADC(read_length);
      
        // now do something with them...
    }
    else if (sense_size == 2) {
        // int-size sense flags
      FREADI(num_ents);
      
        // now do something with them...
    }
  }

  if (sideseth->numElements > 0) {
      // have to read dist factors
    FREADD(sideseth->numElements);
    
      // now do something with them...
  }
}

void Tqdcfr::read_block(Tqdcfr::ModelEntry *model,
                        Tqdcfr::EntityHeader *blockh)  
{
  if (blockh->numTypes == 0) return;
  
    // position file
  FSEEK(model->modelOffset+blockh->nodeOffset);
  
    // read ids for each entity type
  int this_type, num_ents;
  for (int i = 0; i < blockh->numTypes; i++) {
      // get how many and what type
    FREADI(2);
    this_type = int_buf[0];
    num_ents = int_buf[1];

      // now get the ids
    FREADI(num_ents);
    
      // get the entity handles (sets or entities) corresponding to this type
      // and put into this block's set
    add_set_entities(this_type, num_ents, &(int_buf[0]), blockh->setHandle);
  }

    // read attribs if there are any
  if (blockh->numElements > 0) {
    MBTag block_attribs;
    
    FREADD(blockh->numElements);
      // now do something with them...
    MBErrorCode result = mdbImpl->tag_create("Block_Attributes", 
                                              blockh->numElements*sizeof(double), MB_TAG_SPARSE, 
                                              block_attribs, NULL);
    assert(MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);
    result = mdbImpl->tag_set_data(block_attribs, &(blockh->setHandle), 1,
                                   &(dbl_buf[0]));
    assert(MB_SUCCESS == result);
  }
}

void Tqdcfr::add_set_entities(const int this_type, const int num_ents,
                              const int *ids, MBEntityHandle &set_handle) 
{
    // parse the entity type
  MBEntityType set_type = MBMAXTYPE;
  if (0 <= this_type && 3 >= this_type) 
    set_type = MBENTITYSET;
  else if (4 == this_type)
    set_type = MBHEX;
  else if (5 == this_type)
    set_type = MBTET;
  else if (6 == this_type)
    set_type = MBPYRAMID;
  else if (7 == this_type)
    set_type = MBQUAD;
  else if (8 == this_type)
    set_type = MBTRI;
  else if (9 == this_type)
    set_type = MBEDGE;
  else if (10 == this_type)
    set_type = MBVERTEX;
  else {
    assert(false);
    return;
  }
  
  MBRange temp_entities, add_entities;
  MBEntityHandle dum_ent;
  int i;
  MBErrorCode result;
  
  if (set_type == MBENTITYSET) {
      // get an entity set with the right dimension & global id tag
    const void *values[2];
    const MBTag tags[2] = {geomTag, globalIdTag};
    int dum_dim = 3 - this_type;
    values[0] = &dum_dim;

    for (i = 0; i < num_ents; i++) {
      values[1] = ids+i;
      temp_entities.clear();
      result = 
        mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, tags, values, 2, temp_entities);
      assert(MB_FAILURE != result);
      if (temp_entities.size() == 1)
        add_entities.insert(*temp_entities.begin());
      else
        std::cout << "Warning: Tqdcfr::add_set_entities supposed to find a set with geom dimension "
                  << dum_dim << " and id " << *(ids+i) << " but didn't!" << std::endl;
    }
  }
  
  else {
      // else offset the id, then look for that entity
    for (i = 0; i < num_ents; i++) {
      result = mdbImpl->handle_from_id(set_type, int_buf[i], dum_ent);
      if (MB_SUCCESS == result) add_entities.insert(dum_ent);
    }
  }

    // now add the entities to the set
  if ((int) add_entities.size() != num_ents)
    std::cout << "Warning: Tqdcfr::add_entities: didn't find the expected number of entities."
              << std::endl;
  
  if (!add_entities.empty()) {
    result = mdbImpl->add_entities(set_handle, add_entities);
    assert(MB_SUCCESS == result);
  }
}

void Tqdcfr::read_group(Tqdcfr::ModelEntry *model,
                        Tqdcfr::EntityHeader *grouph)  
{
    // position file
  FSEEK(model->modelOffset+grouph->nodeOffset);
  
    // read ids for each entity type
  int this_type, num_ents;
  for (int i = 0; i < grouph->numTypes; i++) {
      // get how many and what type
    FREADI(2);
    this_type = int_buf[0];
    num_ents = int_buf[1];

      // now get the ids
    FREADI(num_ents);
    
      // now do something with them...
  }
}

void Tqdcfr::read_nodes(Tqdcfr::ModelEntry *model,
                        Tqdcfr::EntityHeader *entity) 
{
  if (entity->numNodes == 0) return;
  
    // get the ids & coords in separate calls to minimize memory usage
    // position the file
  FSEEK(model->modelOffset+entity->nodeOffset);
    // get node ids in int_buf
  FREADI(entity->numNodes);

      // check to see if ids are contiguous...
    if (0 == check_contiguous(entity->numNodes))
      std::cout << "Node ids are not contiguous!" << std::endl;
    
    // get a space for reading nodal data directly into MB
  MBEntityHandle node_handle = 0;
  std::vector<double*> arrays;
  readUtilIface->get_node_arrays(3, entity->numNodes,
                                 int_buf[0], node_handle, arrays);
  long unsigned int node_offset;
  node_offset = mdbImpl->id_from_handle( node_handle);
  
  node_offset -= int_buf[0];
  if (-1 == currNodeIdOffset)
    currNodeIdOffset = node_offset;
  else
    assert((long unsigned int) currNodeIdOffset == node_offset);

    // get node x's in arrays[0]
  FREADDA(entity->numNodes, arrays[0]);
    // get node y's in arrays[1]
  FREADDA(entity->numNodes, arrays[1]);
    // get node z's in arrays[2]
  FREADDA(entity->numNodes, arrays[2]);

    // add these nodes into the entity's set
  MBRange dum_range(node_handle, 
                     node_handle+entity->numNodes-1);
  MBErrorCode result = mdbImpl->add_entities(entity->setHandle, dum_range);
  assert(MB_SUCCESS == result);

    // set the dimension to at least zero (entity has at least nodes) on the geom tag
  int max_dim = 0;
  result = mdbImpl->tag_set_data(geomTag, &(entity->setHandle), 1, &max_dim);
  assert(MB_SUCCESS == result);
    // set the category tag just in case there're only vertices in this set
  result = mdbImpl->tag_set_data(categoryTag, &entity->setHandle, 1, 
                                 &geom_categories[0]);
  assert(MB_SUCCESS == result);
  
    // don't bother with FixedNode metadata for now...
}

void Tqdcfr::read_elements(Tqdcfr::ModelEntry *model,
                           Tqdcfr::EntityHeader *entity) 
{
  if (entity->numTypes == 0) return;
  
    // get data in separate calls to minimize memory usage
    // position the file
  FSEEK(model->modelOffset+entity->elementOffset);

  int int_type, nodes_per_elem, num_elem;
  int max_dim = -1;
  MBErrorCode result;
  for (int i = 0; i < entity->numTypes; i++) {
      // for this elem type, get the type, nodes per elem, num elems
    FREADI(3);
    int_type = int_buf[0];
    nodes_per_elem = int_buf[1];
    num_elem = int_buf[2];

      // get MB element type from cub file's 
    MBEntityType elem_type = type_from_cub_type(int_type, nodes_per_elem);
    max_dim = (max_dim < MBCN::Dimension(elem_type) ? MBCN::Dimension(elem_type) : max_dim);
    
      // get element ids
    FREADI(num_elem);
    
      // check to see if ids are contiguous...
    if (0 == check_contiguous(num_elem))
      std::cout << "Element ids are not contiguous!" << std::endl;
    
      // get a space for reading connectivity data directly into MB
    MBEntityHandle *conn, start_handle;
    
    readUtilIface->get_element_array(num_elem, nodes_per_elem,
                                     elem_type, int_buf[0], start_handle, conn);
        
    long unsigned int elem_offset;
    elem_offset = mdbImpl->id_from_handle( start_handle) - int_buf[0];
    if (-1 == currElementIdOffset[elem_type])
      currElementIdOffset[elem_type] = elem_offset;
    else
      ;
    
        //assert((long unsigned int) currElementIdOffset[elem_type] == elem_offset);

      // now do something with them...

      // get the connectivity array
    int total_conn = num_elem * nodes_per_elem;
    FREADIA(total_conn, conn);

      // post-process connectivity into handles
    MBEntityHandle new_node_handle, dum_handle;
    int dum_err;
    for (i = 0; i < total_conn; i++) {
        // do it this way to avoid the cost of checking in optimized code
      new_node_handle = CREATE_HANDLE(MBVERTEX, currNodeIdOffset+conn[i], dum_err);
      assert(MB_SUCCESS == mdbImpl->handle_from_id(MBVERTEX, currNodeIdOffset+conn[i],
                                                    dum_handle));
      conn[i] = new_node_handle;
    }

      // add these elements into the entity's set
    MBRange dum_range(start_handle, start_handle+num_elem-1);
    result = mdbImpl->add_entities(entity->setHandle, dum_range);
    assert(MB_SUCCESS == result);
  }

    // set the dimension on the geom tag
  result = mdbImpl->tag_set_data(geomTag, &entity->setHandle, 1, &max_dim);
  assert(MB_SUCCESS == result);
  if (max_dim != -1) {
    result = mdbImpl->tag_set_data(categoryTag, &entity->setHandle, 1, 
                                   &geom_categories[max_dim]);
    assert(MB_SUCCESS == result);
  }
}

int Tqdcfr::check_contiguous(const int num_ents) 
{
  std::vector<int>::iterator id_it;
  int curr_id, i;

    // check in forward-contiguous direction
  id_it = int_buf.begin();
  curr_id = *id_it++ + 1;
  for (i = 1; id_it != int_buf.end() && i < num_ents; id_it++, i++) {
    if (*id_it != curr_id) {
      i = 0;
      break;
    }
    curr_id++;
  }

    // if we got here and we're at the end of the loop, it's forward-contiguous
  if (i == num_ents) return 1;

// check in reverse-contiguous direction
  id_it = int_buf.begin();
  curr_id = *id_it++ - 1;
  for (i = 1; id_it != int_buf.end() && i < num_ents; id_it++, i++) {
    if (*id_it != curr_id) {
      i = 0;
      break;
    }
    curr_id--;
  }


    // if we got here and we're at the end of the loop, it's reverse-contiguous
  if (i == num_ents) return -1;

    // else it's not contiguous at all
  return 0;
}
  
MBEntityType Tqdcfr::type_from_cub_type(const int cub_type, const int nodes_per_elem) 
{
  bool mid_nodes[3];
  
  MBEntityType ret_type = MBMAXTYPE;
  
    // test the cub_type, and verify it's a valid element
  if (0 <= cub_type && 7 >= cub_type)
    ret_type = MBHEX;
  else if (8 <= cub_type && 15 >= cub_type)
    ret_type = MBTET;
  else if (16 <= cub_type && 19 >= cub_type)
    ret_type = MBPYRAMID;
  else if (20 <= cub_type && 23 >= cub_type)
    ret_type = MBQUAD;
  else if (24 <= cub_type && 27 >= cub_type)
    ret_type = MBTRI;
  else if (28 <= cub_type && 29 >= cub_type)
    ret_type = MBEDGE;
  else if (30 <= cub_type && 30 >= cub_type)
    ret_type = MBVERTEX;
  else 
    assert(false);

    // check to see that the number of nodes passed in represents a valid
    // number of nodes wrt mid nodes
  MBCN::HasMidNodes(ret_type, nodes_per_elem, mid_nodes);
  assert(nodes_per_elem == MBCN::HONodeIndex(ret_type, nodes_per_elem, 
                                               MBCN::Dimension(ret_type), -1)+1);

  return ret_type;
}

void Tqdcfr::FEModelHeader::init(const int offset, Tqdcfr* instance ) 
{
  FSEEK(offset);
  FREADI(4);
  feEndian = instance->int_buf[0];
  feSchema = instance->int_buf[1];
  feCompressFlag = instance->int_buf[2];
  feLength = instance->int_buf[3];
  FREADI(3); geomArray.init(instance->int_buf);
  FREADI(2);
  nodeArray.metaDataOffset = instance->int_buf[0];
  elementArray.metaDataOffset = instance->int_buf[1];
  FREADI(3); groupArray.init(instance->int_buf);
  FREADI(3); blockArray.init(instance->int_buf);
  FREADI(3); nodesetArray.init(instance->int_buf);
  FREADI(3); sidesetArray.init(instance->int_buf);
  FREADI(1);
}

void Tqdcfr::read_file_header() 
{
    // read file header
  FSEEK(4);
  FREADI(6);
  fileTOC.fileEndian = int_buf[0];
  fileTOC.fileSchema = int_buf[1];
  fileTOC.numModels = int_buf[2];
  fileTOC.modelTableOffset = int_buf[3];
  fileTOC.modelMetaDataOffset = int_buf[4];
  fileTOC.activeFEModel = int_buf[5];
  if (debug) fileTOC.print();
}

void Tqdcfr::read_model_entries() 
{
  
    // read model entries
  FSEEK(fileTOC.modelTableOffset);
  FREADI(fileTOC.numModels*6);
  modelEntries = new ModelEntry[fileTOC.numModels];
  assert(NULL != modelEntries);
  std::vector<int>::iterator int_it = int_buf.begin();
  for (int i = 0; i < fileTOC.numModels; i++) {
    modelEntries[i].modelHandle = *int_it++;
    modelEntries[i].modelOffset = *int_it++;
    modelEntries[i].modelLength = *int_it++;
    modelEntries[i].modelType = *int_it++;
    modelEntries[i].modelOwner = *int_it++;
    modelEntries[i].modelPad = *int_it++;
    assert(int_it != int_buf.end() || i == fileTOC.numModels-1);
    if (debug) modelEntries[i].print();
  }
}

int Tqdcfr::find_model(const int model_type) 
{
  for (int i = 0; i < fileTOC.numModels; i++) 
    if (modelEntries[i].modelType == model_type) return i;
  
  return -1;
}

void Tqdcfr::read_meta_data(const int metadata_offset, 
                            Tqdcfr::MetaDataContainer &mc) 
{
    // read the metadata header
  FSEEK(metadata_offset);
  FREADI(3);
  mc.mdSchema = int_buf[0];
  mc.compressFlag = int_buf[1];
  mc.numDatums = int_buf[2];

    // allocate space for the entries
  mc.metadataEntries = 
    new Tqdcfr::MetaDataContainer::MetaDataEntry[mc.numDatums];
  
    // now read the metadata values
  for (int i = 0; i < mc.numDatums; i++) {
    FREADI(2);
    mc.metadataEntries[i].mdOwner = int_buf[0];
    mc.metadataEntries[i].mdDataType = int_buf[1];
    
      // read the name string
    read_md_string(mc.metadataEntries[i].mdName);

    if (mc.metadataEntries[i].mdDataType == 0) {
        // integer
      FREADI(1);
      mc.metadataEntries[i].mdIntValue = int_buf[0];
    }
    else if (mc.metadataEntries[i].mdDataType == 1) {
        // string
      read_md_string(mc.metadataEntries[i].mdStringValue);
    }
    else if (mc.metadataEntries[i].mdDataType == 2) {
        // double
      FREADD(1);
      mc.metadataEntries[i].mdDblValue = dbl_buf[0];
    }
    else if (mc.metadataEntries[i].mdDataType == 3) {
        // int array
      FREADI(1);
      mc.metadataEntries[i].mdIntArrayValue.resize(int_buf[0]);
      FREADI(mc.metadataEntries[i].mdIntArrayValue.size());
      std::copy(int_buf.begin(), 
                int_buf.begin() + mc.metadataEntries[i].mdIntArrayValue.size(),
                std::back_inserter(mc.metadataEntries[i].mdIntArrayValue));
    }
    else if (mc.metadataEntries[i].mdDataType == 4) {
        // double array
      FREADI(1);
      mc.metadataEntries[i].mdDblArrayValue.resize(int_buf[0]);
      FREADD(mc.metadataEntries[i].mdDblArrayValue.size());
      std::copy(dbl_buf.begin(), 
                dbl_buf.begin() + mc.metadataEntries[i].mdDblArrayValue.size(),
                std::back_inserter(mc.metadataEntries[i].mdDblArrayValue));
    }
    else
      assert("Bad metadata type" && false);
  }
  if (debug) mc.print();
}

void Tqdcfr::read_md_string(std::string &name) 
{
  FREADI(1);
  int str_size = int_buf[0];
  if (str_size > 0) {
    FREADC(str_size);
    if (char_buf.size() <= (unsigned int) str_size)
      char_buf.resize(str_size+1);
    char_buf[str_size] = '\0';
    name = (char *) &char_buf[0];
      // read pad if any
    int num_word = str_size/4;
    if (4*num_word != str_size) {
        // read extra chars to end of pad
      str_size = (num_word+1)*4 - str_size;
      FREADC(str_size);
    }
  }
}
  
void Tqdcfr::EntityHeader::read_info_header(const int model_offset, 
                                            const Tqdcfr::FEModelHeader::ArrayInfo &info,
                                            const int info_type,
                                            Tqdcfr* instance,
                                            Tqdcfr::EntityHeader *&entity_headers) 
{
  entity_headers = new EntityHeader[info.numEntities];
  FSEEK(model_offset+info.tableOffset);
  int dum_int;
  MBErrorCode result;

  if (0 == instance->categoryTag) {
    static const char val[CATEGORY_TAG_NAME_LENGTH] = "\0";
    result = instance->mdbImpl->tag_create(CATEGORY_TAG_NAME, CATEGORY_TAG_NAME_LENGTH,
                                           MB_TAG_SPARSE, instance->categoryTag, val);
    assert (MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);
  }

  for (int i = 0; i < info.numEntities; i++) {

      // create an entity set for this entity
    result = instance->mdbImpl->create_meshset(MESHSET_SET, entity_headers[i].setHandle);
    assert(MB_SUCCESS == result);
    static const char group_category[CATEGORY_TAG_NAME_LENGTH] = "Group\0";
    static const char material_category[CATEGORY_TAG_NAME_LENGTH] = "Material Set\0";
    static const char dirichlet_category[CATEGORY_TAG_NAME_LENGTH] = "Dirichlet Set\0";
    static const char neumann_category[CATEGORY_TAG_NAME_LENGTH] = "Neumann Set\0";
    
    switch (info_type) {
      case geom:
        FREADI(8);
        entity_headers[i].numNodes = instance->int_buf[0];
        entity_headers[i].nodeOffset = instance->int_buf[1];
        entity_headers[i].numElements = instance->int_buf[2];
        entity_headers[i].elementOffset = instance->int_buf[3];
        entity_headers[i].numTypes = instance->int_buf[4];
        entity_headers[i].entityLength = instance->int_buf[5];
        entity_headers[i].entityID = instance->int_buf[6];
        entity_headers[i].entityType = 
          entity_headers[i].entityColor = 
          entity_headers[i].pyrType = 
          entity_headers[i].matType = 
          entity_headers[i].blockDimension = 
          entity_headers[i].shellsFlag = -1;

          // set the dimension to -1; will have to reset later, after elements are read
        dum_int = -1;
        result = instance->mdbImpl->tag_set_data(instance->geomTag, 
                                                           &(entity_headers[i].setHandle), 1, &dum_int);
        assert(MB_SUCCESS == result);

          // set a unique id tag
        result = instance->mdbImpl->tag_set_data(instance->uniqueIdTag, 
                                                           &(entity_headers[i].setHandle), 1, 
                                                           &(entity_headers[i].entityID));
        assert(MB_SUCCESS == result);
        
        break;
      case group:
        FREADI(6);
        entity_headers[i].entityID = instance->int_buf[0];
        entity_headers[i].entityType = instance->int_buf[1];
        entity_headers[i].numNodes = instance->int_buf[2];
        entity_headers[i].nodeOffset = instance->int_buf[3];
        entity_headers[i].numTypes = instance->int_buf[4];
        entity_headers[i].entityLength = instance->int_buf[5];
        entity_headers[i].numElements = 
          entity_headers[i].elementOffset = 
          entity_headers[i].entityColor = 
          entity_headers[i].pyrType = 
          entity_headers[i].matType = 
          entity_headers[i].blockDimension = 
          entity_headers[i].shellsFlag = -1;

          // set the group tag to 1 to signify this is a group
        dum_int = 1;
        result = instance->mdbImpl->tag_set_data(instance->groupTag, 
                                                           &(entity_headers[i].setHandle), 1, &dum_int);
        assert(MB_SUCCESS == result);
        result = instance->mdbImpl->tag_set_data(instance->categoryTag, 
                                                 &(entity_headers[i].setHandle), 1, 
                                                 group_category);
        assert(MB_SUCCESS == result);

          // set a global id tag
        result = instance->mdbImpl->tag_set_data(instance->globalIdTag, 
                                                           &(entity_headers[i].setHandle), 1, 
                                               &(entity_headers[i].entityID));
        assert(MB_SUCCESS == result);
        
        break;
      case block:
        FREADI(12);
        entity_headers[i].entityID = instance->int_buf[0];
        entity_headers[i].entityType = instance->int_buf[1];
        entity_headers[i].numNodes = instance->int_buf[2];
        entity_headers[i].nodeOffset = instance->int_buf[3];
        entity_headers[i].numTypes = instance->int_buf[4];
        entity_headers[i].numElements = instance->int_buf[5]; // attrib order
        entity_headers[i].entityColor = instance->int_buf[6];
        entity_headers[i].elementOffset = instance->int_buf[7]; // mixed elem type
        entity_headers[i].pyrType = instance->int_buf[8];
        entity_headers[i].matType = instance->int_buf[9];
        entity_headers[i].entityLength = instance->int_buf[10];
        entity_headers[i].blockDimension = instance->int_buf[11];
        entity_headers[i].shellsFlag = -1;

          // set the material set tag and id tag both to id
        result = instance->mdbImpl->tag_set_data(instance->blockTag, &(entity_headers[i].setHandle), 1, 
                                               &(entity_headers[i].entityID));
        assert(MB_SUCCESS == result);
        result = instance->mdbImpl->tag_set_data(instance->globalIdTag, &(entity_headers[i].setHandle), 1, 
                                               &(entity_headers[i].entityID));
        assert(MB_SUCCESS == result);
        result = instance->mdbImpl->tag_set_data(instance->categoryTag, 
                                                 &(entity_headers[i].setHandle), 1, 
                                                 material_category);
        assert(MB_SUCCESS == result);
        
        break;
      case nodeset:
        FREADI(8);
        entity_headers[i].entityID = instance->int_buf[0];
        entity_headers[i].numNodes = instance->int_buf[1];
        entity_headers[i].nodeOffset = instance->int_buf[2];
        entity_headers[i].numTypes = instance->int_buf[3];
        entity_headers[i].pyrType = instance->int_buf[4];  // point sym
        entity_headers[i].entityColor = instance->int_buf[5];
        entity_headers[i].entityLength = instance->int_buf[6];
          // pad

        entity_headers[i].entityType =
          entity_headers[i].numElements =
          entity_headers[i].elementOffset =
          entity_headers[i].matType =
          entity_headers[i].blockDimension =
          entity_headers[i].shellsFlag = -1;

          // set the dirichlet set tag and id tag both to id
        result = instance->mdbImpl->tag_set_data(instance->nsTag, &(entity_headers[i].setHandle), 1, 
                                               &(entity_headers[i].entityID));
        assert(MB_SUCCESS == result);
        result = instance->mdbImpl->tag_set_data(instance->globalIdTag, &(entity_headers[i].setHandle), 1, 
                                               &(entity_headers[i].entityID));
        assert(MB_SUCCESS == result);
        result = instance->mdbImpl->tag_set_data(instance->categoryTag, 
                                                 &(entity_headers[i].setHandle), 1, 
                                                 dirichlet_category);
        assert(MB_SUCCESS == result);
        
        break;
      case sideset:
        FREADI(8);
        entity_headers[i].entityID = instance->int_buf[0];
        entity_headers[i].numNodes = instance->int_buf[1];
        entity_headers[i].nodeOffset = instance->int_buf[2];
        entity_headers[i].numTypes = instance->int_buf[3];
        entity_headers[i].numElements = instance->int_buf[4]; // num dist factors
        entity_headers[i].entityColor = instance->int_buf[5];
        entity_headers[i].shellsFlag = instance->int_buf[6];
        entity_headers[i].entityLength = instance->int_buf[7];

        entity_headers[i].pyrType =
          entity_headers[i].entityType =
          entity_headers[i].elementOffset =
          entity_headers[i].pyrType =
          entity_headers[i].matType =
          entity_headers[i].blockDimension = -1;

          // set the neumann set tag and id tag both to id
        result = instance->mdbImpl->tag_set_data(instance->ssTag, &(entity_headers[i].setHandle), 1, 
                                               &(entity_headers[i].entityID));
        assert(MB_SUCCESS == result);
        result = instance->mdbImpl->tag_set_data(instance->globalIdTag, &(entity_headers[i].setHandle), 1, 
                                               &(entity_headers[i].entityID));
        assert(MB_SUCCESS == result);
        result = instance->mdbImpl->tag_set_data(instance->categoryTag, 
                                                 &(entity_headers[i].setHandle), 1, 
                                                 neumann_category);
        assert(MB_SUCCESS == result);
        
        break;
      default:
        assert(false);
    }
  }
}

void Tqdcfr::ModelEntry::print_header(const char *prefix,
                                      Tqdcfr::FEModelHeader::ArrayInfo &info,
                                      Tqdcfr::EntityHeader *header) 
{
  if (!debug) return;
  std::cout << prefix << std::endl;
  for (int i = 0; i < info.numEntities; i++)
    header[i].print();
}
          
void Tqdcfr::ModelEntry::read_header_info( Tqdcfr* instance )
{
  feModelHeader.init(modelOffset, instance);
  int default_val = -1;
  MBErrorCode result;

  result = instance->mdbImpl->tag_create(GLOBAL_ID_TAG_NAME, 4, MB_TAG_DENSE, 
                                         instance->globalIdTag, &default_val);
  assert(MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);

  if (feModelHeader.geomArray.numEntities > 0) {
    result = instance->mdbImpl->tag_create(GEOM_DIMENSION_TAG_NAME, 4, MB_TAG_SPARSE, 
                                                instance->geomTag, &default_val);
    assert(MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);
    
    result = instance->mdbImpl->tag_create("UNIQUE_ID", 4, MB_TAG_SPARSE, 
                                                            instance->uniqueIdTag, &default_val);
    assert(MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);
    
    Tqdcfr::EntityHeader::read_info_header(modelOffset, 
                                           feModelHeader.geomArray, 
                                           Tqdcfr::EntityHeader::geom, 
                                           instance,
                                           feGeomH);
    print_header("Geom headers:", feModelHeader.geomArray,
                 feGeomH);
  }
  
  if (feModelHeader.groupArray.numEntities > 0) {
    result = instance->mdbImpl->tag_create("GROUP_SET", 4, MB_TAG_SPARSE, 
                                                instance->groupTag, &default_val);
    assert(MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);
    
    Tqdcfr::EntityHeader::read_info_header(modelOffset, 
                                           feModelHeader.groupArray, 
                                           Tqdcfr::EntityHeader::group,
                                           instance,
                                           feGroupH);
    print_header("Group headers:", feModelHeader.groupArray,
                 feGroupH);
  }

  if (feModelHeader.blockArray.numEntities > 0) {
    result = instance->mdbImpl->tag_create(MATERIAL_SET_TAG_NAME, 4, MB_TAG_SPARSE, 
                                                instance->blockTag, &default_val);
    assert(MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);
    
    Tqdcfr::EntityHeader::read_info_header(modelOffset, 
                                           feModelHeader.blockArray, 
                                           Tqdcfr::EntityHeader::block, 
                                           instance,
                                           feBlockH);
    print_header("Block headers:", feModelHeader.blockArray,
                 feBlockH);
  }
  if (feModelHeader.nodesetArray.numEntities > 0) {
    result = instance->mdbImpl->tag_create(DIRICHLET_SET_TAG_NAME, 4, MB_TAG_SPARSE, 
                                                instance->nsTag, &default_val);
    assert(MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);
    
    Tqdcfr::EntityHeader::read_info_header(modelOffset, 
                                           feModelHeader.nodesetArray, 
                                           Tqdcfr::EntityHeader::nodeset, 
                                           instance,
                                           feNodeSetH);
    print_header("Nodeset headers:", feModelHeader.nodesetArray,
                 feNodeSetH);
  }
  if (feModelHeader.sidesetArray.numEntities > 0) {
    result = instance->mdbImpl->tag_create(NEUMANN_SET_TAG_NAME, 4, MB_TAG_SPARSE, 
                                                instance->ssTag, &default_val);
    assert(MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);
    
    Tqdcfr::EntityHeader::read_info_header(modelOffset, 
                                           feModelHeader.sidesetArray, 
                                           Tqdcfr::EntityHeader::sideset, 
                                           instance,
                                           feSideSetH);
    print_header("SideSet headers:", feModelHeader.sidesetArray,
                 feSideSetH);
  }

}

void Tqdcfr::ModelEntry::read_metadata_info(Tqdcfr *tqd) 
{
  if (debug) std::cout << "Geom metadata:" << std::endl;
  tqd->read_meta_data(modelOffset+feModelHeader.geomArray.metaDataOffset,
                      geomMD);
  if (debug) std::cout << "Node metadata:" << std::endl;
  tqd->read_meta_data(modelOffset+feModelHeader.nodeArray.metaDataOffset,
                      nodeMD);
  if (debug) std::cout << "Elem metadata:" << std::endl;
  tqd->read_meta_data(modelOffset+feModelHeader.elementArray.metaDataOffset,
                      elementMD);
  if (debug) std::cout << "Group metadata:" << std::endl;
  tqd->read_meta_data(modelOffset+feModelHeader.groupArray.metaDataOffset,
                      groupMD);
  if (debug) std::cout << "Block metadata:" << std::endl;
  tqd->read_meta_data(modelOffset+feModelHeader.blockArray.metaDataOffset,
                      blockMD);
  if (debug) std::cout << "Nodeset metadata:" << std::endl;
  tqd->read_meta_data(modelOffset+feModelHeader.nodesetArray.metaDataOffset,
                      nodesetMD);
  if (debug) std::cout << "Sideset metadata:" << std::endl;
  tqd->read_meta_data(modelOffset+feModelHeader.sidesetArray.metaDataOffset,
                      sidesetMD);
}

void Tqdcfr::read_acis_records() 
{

    // get the acis model location
  int acis_model_offset = 0, acis_model_length = 0, acis_model_handle = 1,
    acis_sat_type = 1;
  for (int i = 0; i < fileTOC.numModels; i++) {
    if (modelEntries[i].modelHandle == acis_model_handle &&
        modelEntries[i].modelType == acis_sat_type) {
      acis_model_offset = modelEntries[i].modelOffset;
      acis_model_length = modelEntries[i].modelLength;
      break;
    }
  }
  
  if (acis_model_length == 0) return;
  
  std::vector<AcisRecord> records;

    // get name for acis dump file
  MBTag acis_file_tag;
  MBErrorCode rval;
  rval = mdbImpl->tag_get_handle( acis_dump_file_tag_name, acis_file_tag );
  const char* filename = default_acis_dump_file;
  if (MB_SUCCESS == rval)
    mdbImpl->tag_get_data( acis_file_tag, 0, 0, &filename );
  
  acisDumpFile = NULL;
  if (filename && *filename)
  {
    acisDumpFile = fopen( filename, "w+" );
  }

    // position the file at the start of the acis model
  FSEEK(acis_model_offset);

  int bytes_left = acis_model_length;
  
  struct AcisRecord this_record;
  reset_record(this_record);
  char *ret;

    // make the char buffer at least buf_size+1 long, to fit null char
  const int buf_size = 1023;
  
  CHECK_SIZE(char_buf, buf_size+1);
  
  while (0 != bytes_left) {
      // read the next buff characters, or bytes_left if smaller
    int next_buf = (bytes_left > buf_size ? buf_size : bytes_left);
    FREADC(next_buf);

    if (NULL != acisDumpFile)
      fwrite(&char_buf[0], sizeof(char), next_buf, acisDumpFile);
    
      // put null at end of string to stop searches 
    char_buf[next_buf] = '\0';
    int buf_pos = 0;

      // check for first read, and if so, get rid of the header
    if (bytes_left == acis_model_length) {
        // look for 3 newlines
      ret = strchr(&(char_buf[0]), '\n'); ret = strchr(ret+1, '\n'); ret = strchr(ret+1, '\n');
      assert(NULL != ret);
      buf_pos += ret - &(char_buf[0]) + 1;
    }
      
    bytes_left -= next_buf;

      // now start grabbing records
    do {
      
        // get next occurrence of '#' (record terminator)
      ret = strchr(&(char_buf[buf_pos]), '#');
      if (NULL != ret) {
          // grab the string (inclusive of the record terminator and the line feed) and complete the record
        int num_chars = ret-&(char_buf[buf_pos])+2;
        this_record.att_string.append(&(char_buf[buf_pos]), num_chars);
        buf_pos += num_chars;
        process_record(this_record);

          // put the record in the list...
        records.push_back(this_record);

          // and reset the record
        reset_record(this_record);
      }
      else {
          // reached end of buffer; cache string then go get another; discard last character,
          // which will be the null character
        this_record.att_string.append(&(char_buf[buf_pos]), next_buf-buf_pos);
        buf_pos = next_buf;
      }
      
    }
    while (buf_pos < next_buf);
  }

  if (NULL != acisDumpFile)
    fwrite("\n======================\nSorted acis records:\n======================\n", 1, 68, acisDumpFile);
    
    // now interpret the records
  interpret_acis_records(records);
  
  if (NULL != acisDumpFile)
    fclose(acisDumpFile);
}

void Tqdcfr::interpret_acis_records(std::vector<AcisRecord> &records) 
{
    // make a tag for the vector holding unrecognized attributes
  void *default_val = NULL;
  MBErrorCode result = 
    mdbImpl->tag_create("ATTRIB_VECTOR", sizeof(void*), MB_TAG_SPARSE, 
                         attribVectorTag, &default_val);
  assert(MB_SUCCESS == result || MB_ALREADY_ALLOCATED == result);

  int current_record = 0;

#define REC records[current_record]

  while (current_record != (int) records.size()) {

      // if this record's been processed, or if it's an attribute, continue
    if (REC.processed || REC.rec_type == Tqdcfr::ATTRIB) {
      current_record++;
      continue;
    }

    if (REC.rec_type == Tqdcfr::UNKNOWN) {
      REC.processed = true;
      current_record++;
      continue;
    }
    
      // it's a known, non-attrib rec type; parse for any attribs
    parse_acis_attribs(current_record, records);

    REC.processed = true;
    
    current_record++;
  }
}

void Tqdcfr::parse_acis_attribs(const int entity_rec_num,
                                std::vector<AcisRecord> &records) 
{
  int num_read;
  std::vector<std::string> *attrib_vec = NULL;
  char temp_name[80], *name_tag = NULL;
  int id = -1;
  int uid = -1;
  int next_attrib = -1;
  MBErrorCode result;
  
  int current_attrib = records[entity_rec_num].first_attrib;
  if (-1 == current_attrib) return;

  if (NULL != acisDumpFile) {
    fwrite("-----------------------------------------------------------------------\n", 1, 72, acisDumpFile);
    fwrite(records[entity_rec_num].att_string.c_str(), sizeof(char), 
           records[entity_rec_num].att_string.length(), acisDumpFile);
  }

  while (-1 != current_attrib) {
    assert(records[current_attrib].rec_type == Tqdcfr::UNKNOWN ||
           (records[current_attrib].att_next == next_attrib &&
            records[current_attrib].att_ent_num == entity_rec_num));
    
    if (NULL != acisDumpFile)
      fwrite(records[current_attrib].att_string.c_str(), sizeof(char), 
             records[current_attrib].att_string.length(), acisDumpFile);

      // is the attrib one we already recognize?
    if (strncmp(records[current_attrib].att_string.c_str(), "ENTITY_NAME", 11) == 0) {
        // parse name
      int num_chars;
      num_read = sscanf(records[current_attrib].att_string.c_str(), "ENTITY_NAME @%d %s", &num_chars, temp_name);
      if (num_read != 2)
        num_read = sscanf(records[current_attrib].att_string.c_str(), "ENTITY_NAME %d %s", &num_chars, temp_name);
      assert(num_read == 2);

        // put the name on the entity
      name_tag = new char[num_chars];
      strcpy(name_tag, temp_name);
    }
    else if (strncmp(records[current_attrib].att_string.c_str(), "ENTITY_ID", 9) == 0) {
        // parse id
      int bounding_uid, bounding_sense;
      num_read = sscanf(records[current_attrib].att_string.c_str(), "ENTITY_ID 0 3 %d %d %d", 
                        &id, &bounding_uid, &bounding_sense);
      if (3 != num_read)
        std::cout << "Warning: bad ENTITY_ID attribute in .sat file, record number " << entity_rec_num
                  << ", record follows:" << std::endl
                  << records[current_attrib].att_string.c_str() << std::endl;
;
    }
    else if (strncmp(records[current_attrib].att_string.c_str(), "UNIQUE_ID", 9) == 0) {
        // parse uid
      num_read = sscanf(records[current_attrib].att_string.c_str(), "UNIQUE_ID 1 0 1 %d", &uid);
      assert(1 == num_read);
    }
    else {
      if (attrib_vec == NULL) attrib_vec = new std::vector<std::string>;
      attrib_vec->push_back(records[current_attrib].att_string);
    }

    records[current_attrib].processed = true;
    next_attrib = current_attrib;
    current_attrib = records[current_attrib].att_prev;
  }

    // at this point, there aren't entity sets for entity types which don't contain mesh
    // in this case, just return
  if (records[entity_rec_num].rec_type == BODY)
    return;

  else if (records[entity_rec_num].entity == 0 && uid == -1) {
    std::cout << "Warning: couldn't resolve entity of type " << records[entity_rec_num].rec_type
              << " because no uid was found." << std::endl;
    return;
  }
  
    // parsed the data; now put on mdb entities; first we need to find the entity
  if (records[entity_rec_num].entity == 0) {
      // get the choices of entity
    MBRange entities;
    assert(uid != -1);
    const void *dum_ptr = &uid;
    result = mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &uniqueIdTag, 
                                                   &dum_ptr, 1, entities);
    assert(MB_SUCCESS == result);
    if (entities.size() != 1) return;
    else records[entity_rec_num].entity = *entities.begin();
  }
  
    // set the id
  if (id != -1) {
    result = mdbImpl->tag_set_data(globalIdTag, &(records[entity_rec_num].entity), 1, &id);
    assert(MB_SUCCESS == result);
  }
  
    // set the name
  if (NULL != name_tag) {
    if (0 == entityNameTag) {
      result = mdbImpl->tag_get_handle("NAME", entityNameTag);
      if (MB_SUCCESS != result || 0 == entityNameTag) {
        char *dum_val = NULL;
        result = mdbImpl->tag_create("NAME", sizeof(char*), MB_TAG_SPARSE, 
                                     entityNameTag, &dum_val);
      }
    }
    assert(0 != entityNameTag);

    result = mdbImpl->tag_set_data(entityNameTag, &(records[entity_rec_num].entity), 1, &name_tag);
    assert(MB_SUCCESS == result);
  }

  if (NULL != attrib_vec) {
      // put the attrib vector in a tag on the entity
    std::vector<std::string> *dum_vec;
    result = mdbImpl->tag_get_data(attribVectorTag, &(records[entity_rec_num].entity), 1, &dum_vec);
    assert(MB_SUCCESS == result || MB_TAG_NOT_FOUND == result);
    if (MB_TAG_NOT_FOUND == result || dum_vec == NULL) {
        // put this list directly on the entity
      result = mdbImpl->tag_set_data(attribVectorTag, &(records[entity_rec_num].entity), 1, &attrib_vec);
      assert(MB_SUCCESS == result);
    }
    else {
        // copy this list over, and delete this list
      std::copy(attrib_vec->begin(), attrib_vec->end(), 
                std::back_inserter(*dum_vec));
      delete attrib_vec;
    }
  }
  
}

void Tqdcfr::reset_record(AcisRecord &this_record) 
{
  this_record.rec_type = Tqdcfr::UNKNOWN;
  static std::string blank;
  this_record.att_string = blank;
  this_record.first_attrib = this_record.att_prev = 
    this_record.att_next = this_record.att_ent_num = -1;
  this_record.processed = false;
  this_record.entity = 0;
}
  
void Tqdcfr::process_record(AcisRecord &this_record)
{
    // get the entity type
  const char *type_substr;

    // try attribs first, since the others have some common processing between them
  if ((type_substr = strstr(this_record.att_string.c_str(), "attrib")) != NULL && 
      type_substr-this_record.att_string.c_str() < 20) {
    this_record.rec_type = Tqdcfr::ATTRIB;
    bool simple_attrib = false;
    if ((type_substr = strstr(this_record.att_string.c_str(), "simple-snl-attrib")) != NULL)
      simple_attrib = true;
    else {
      this_record.rec_type = Tqdcfr::UNKNOWN;
      return;
    }

      // find next space
    type_substr = strchr(type_substr, ' ');
    assert(NULL != type_substr);
    
      // read the numbers from there
    int num_converted = sscanf(type_substr, " $-1 -1 $%d $%d $%d -1", &(this_record.att_prev), 
                               &(this_record.att_next), &(this_record.att_ent_num));
    assert(num_converted == 3);
    
      // trim the string to the attribute, if it's a simple attrib
    if (simple_attrib) {
      type_substr = strstr(this_record.att_string.c_str(), "NEW_SIMPLE_ATTRIB");
      assert(NULL != type_substr);
      type_substr = strstr(type_substr, "@");
      assert(NULL != type_substr);
      type_substr = strstr(type_substr, " ") + 1;
      assert(NULL != type_substr);
        // copy the rest of the string to a dummy string
      std::string dum_str(type_substr);
      this_record.att_string = dum_str;
    }
  }
  else {
      // else it's a topological entity, I think
    if ((type_substr = strstr(this_record.att_string.c_str(), "body")) != NULL 
        && type_substr-this_record.att_string.c_str() < 20) {
      this_record.rec_type = Tqdcfr::BODY;
    }
    else if ((type_substr = strstr(this_record.att_string.c_str(), "lump")) != NULL  && 
             type_substr-this_record.att_string.c_str() < 20) {
      this_record.rec_type = Tqdcfr::LUMP;
    }
    else if ((type_substr = strstr(this_record.att_string.c_str(), "shell")) != NULL && 
             type_substr-this_record.att_string.c_str() < 20) {
        // don't care about shells
      this_record.rec_type = Tqdcfr::UNKNOWN;
    }
    else if ((type_substr = strstr(this_record.att_string.c_str(), "surface")) != NULL && 
             type_substr-this_record.att_string.c_str() < 20) {
        // don't care about surfaces
      this_record.rec_type = Tqdcfr::UNKNOWN;
    }
    else if ((type_substr = strstr(this_record.att_string.c_str(), "face")) != NULL && 
             type_substr-this_record.att_string.c_str() < 20) {
      this_record.rec_type = Tqdcfr::FACE;
    }
    else if ((type_substr = strstr(this_record.att_string.c_str(), "loop")) != NULL && 
             type_substr-this_record.att_string.c_str() < 20) {
        // don't care about loops
      this_record.rec_type = Tqdcfr::UNKNOWN;
    }
    else if ((type_substr = strstr(this_record.att_string.c_str(), "coedge")) != NULL && 
             type_substr-this_record.att_string.c_str() < 20) {
        // don't care about coedges
      this_record.rec_type = Tqdcfr::UNKNOWN;
    }
    else if ((type_substr = strstr(this_record.att_string.c_str(), "edge")) != NULL && 
             type_substr-this_record.att_string.c_str() < 20) {
      this_record.rec_type = Tqdcfr::EDGE;
    }
    else if ((type_substr = strstr(this_record.att_string.c_str(), "vertex")) != NULL && 
             type_substr-this_record.att_string.c_str() < 20) {
      this_record.rec_type = Tqdcfr::VERTEX;
    }
    else 
      this_record.rec_type = Tqdcfr::UNKNOWN;
    
    if (this_record.rec_type != Tqdcfr::UNKNOWN) {

        // print a warning if it looks like there are sequence numbers
      if (type_substr != this_record.att_string.c_str())
        std::cout << "Warning: acis file has sequence numbers!" << std::endl;

        // scan ahead to the next white space
      type_substr = strchr(type_substr, ' ');
      assert(NULL != type_substr);
      
        // get the id of the first attrib
      int num_converted = sscanf(type_substr, " $%d", &(this_record.first_attrib));
      assert(num_converted == 1);
    }
  }
}

#ifdef TEST_TQDCFR
#include "MBCore.hpp"
int main(int argc, char* argv[])
{

    // Check command line arg
  if (argc < 2)
  {
    std::cout << "Usage: tqdcfr <cub_file_name>" << std::endl;
    exit(1);
  }

  MBCore my_impl;
  mdbImpl = &my_impl;
  Tqdcfr my_tqd(&my_impl);

  my_tqd.load_file(std::string(argv[1]));
  
}
#endif
