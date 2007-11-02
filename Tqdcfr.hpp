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

/*
 * Tim's Quick 'N Dirty Cub File Reader (Tqdcfr)
 *
 */

#ifndef TQDCFR
#define TQDCFR

#include "MBForward.hpp"
#include "MBReaderIface.hpp"
#include "MBTagConventions.hpp"

#include <stdio.h>
#include <string>
#include <vector>
#include <map>

class MBReadUtilIface;
class FEModelHeader;
class GeomHeader;
class GroupHeader;
class BlockHeader;
class NodesetHeader;
class SidesetHeader;


class Tqdcfr : public MBReaderIface
{
public:  

  void FSEEK( long offset );        // set cubFile offset to specified value
  void FREADI( unsigned num_ents ); // read integers into int_buf
  void FREADD( unsigned num_ents ); // read doubles into dbl_buf
  void FREADC( unsigned num_ents ); // read characters into char_buf
  void FREADIA( unsigned num_ents,    int* array ); // read integers
  void FREADDA( unsigned num_ents, double* array ); // read doubles
  void FREADCA( unsigned num_ents,   char* arrat ); // read bytes

    // class for holding the file table of contents
  class FileTOC
  {
  public:
    int fileEndian, fileSchema, numModels, modelTableOffset, 
      modelMetaDataOffset, activeFEModel;

    FileTOC();
    void print();
  };

    // 
  class FEModelHeader 
  {
  public:
    int feEndian, feSchema, feCompressFlag, feLength;

    class ArrayInfo 
    {
    public:
      int numEntities, tableOffset, metaDataOffset;

      ArrayInfo();
      
      void print();
      void init(const std::vector<int>& int_buf);
    };
    
    ArrayInfo geomArray, nodeArray, elementArray, groupArray, 
      blockArray, nodesetArray, sidesetArray;

    void init(const int offset, Tqdcfr* instance );
        
    void print();
  };

  class MetaDataContainer
  {
  public:
    int mdSchema, compressFlag, numDatums;

    class MetaDataEntry
    {
    public:
      int mdOwner, mdDataType, mdIntValue;
      std::string mdName, mdStringValue;
      std::vector<int> mdIntArrayValue;
      double mdDblValue;
      std::vector<double> mdDblArrayValue;
      
      MetaDataEntry();

      void print();
    };

    void print();

    int get_md_entry(const int owner, const std::string &name);
    
    MetaDataEntry *metadataEntries;
    MetaDataContainer();
    ~MetaDataContainer();
  };

  class GeomHeader 
  {
  public:
    int geomID, nodeCt, nodeOffset, elemCt, elemOffset, 
      elemTypeCt, elemLength;

    MBEntityHandle setHandle;

    void print();

    static MBErrorCode read_info_header(const int model_offset, 
                                        const FEModelHeader::ArrayInfo &info,
                                        Tqdcfr* instance,
                                        GeomHeader *&entity_headers);

    GeomHeader();
  };
  
  class GroupHeader 
  {
  public:
    int grpID, grpType, memCt, memOffset, memTypeCt, grpLength;

    MBEntityHandle setHandle;

    void print();

    static MBErrorCode read_info_header(const int model_offset, 
                                 const FEModelHeader::ArrayInfo &info,
                                 Tqdcfr* instance,
                                 GroupHeader *&entity_headers);

    GroupHeader();
  };
  
  class BlockHeader 
  {
  public:
    int blockID, blockElemType, memCt, memOffset, memTypeCt, attribOrder, blockCol,
      blockMixElemType, blockPyrType, blockMat, blockLength, blockDim;

    MBEntityHandle setHandle;

    MBEntityType blockEntityType;

    void print();

    static MBErrorCode read_info_header(const double data_version,
                                        const int model_offset, 
                                        const FEModelHeader::ArrayInfo &info,
                                        Tqdcfr* instance,
                                        BlockHeader *&block_headers);

    BlockHeader();
  };
  
  class NodesetHeader 
  {
  public:
    int nsID, memCt, memOffset, memTypeCt, pointSym, nsCol, nsLength;

    MBEntityHandle setHandle;

    void print();

    static MBErrorCode read_info_header(const int model_offset, 
                                 const FEModelHeader::ArrayInfo &info,
                                 Tqdcfr* instance,
                                 NodesetHeader *&entity_headers);

    NodesetHeader();
  };
  
  class SidesetHeader 
  {
  public:
    int ssID, memCt, memOffset, memTypeCt, numDF, ssCol, useShell, ssLength;

    MBEntityHandle setHandle;

    void print();

    static MBErrorCode read_info_header(const int model_offset, 
                                 const FEModelHeader::ArrayInfo &info,
                                 Tqdcfr* instance,
                                 SidesetHeader *&entity_headers);

    SidesetHeader();
  };
  
    // class to hold model entry data for various kinds of models 
    // (acis, free mesh, etc.)
  class ModelEntry
  {
  public:
    ModelEntry();

    ~ModelEntry();
    
    int modelHandle, modelOffset, modelLength, modelType, modelOwner, modelPad;

    FEModelHeader feModelHeader;
    GeomHeader *feGeomH;
    GroupHeader *feGroupH;
    BlockHeader *feBlockH;
    NodesetHeader *feNodeSetH;
    SidesetHeader *feSideSetH;
    
    MetaDataContainer geomMD, nodeMD, elementMD, groupMD, blockMD, nodesetMD, sidesetMD;
    
    void print();

    void print_geom_headers(const char *prefix,
                            GeomHeader *header,
                            int num_headers);

    void print_group_headers(const char *prefix,
                             GroupHeader *header,
                             const int num_headers);

    void print_block_headers(const char *prefix,
                             BlockHeader *header,
                             const int num_headers);

    void print_nodeset_headers(const char *prefix,
                               NodesetHeader *header,
                               const int num_headers);

    void print_sideset_headers(const char *prefix,
                               SidesetHeader *header,
                               const int num_headers);
    
    MBErrorCode read_header_info( Tqdcfr* instance, const double data_version);
    MBErrorCode read_metadata_info(Tqdcfr *tqd);
  };

  enum {aBODY, LUMP, SHELL, FACE, LOOP, COEDGE, aEDGE, aVERTEX, ATTRIB, UNKNOWN};
  
  const int *ACIS_DIMS;
  
  struct AcisRecord 
  {
    int rec_type;
    std::string att_string;
    bool processed;
    int first_attrib;
    int att_prev, att_next, att_ent_num;
    MBEntityHandle entity;
  };

  ~Tqdcfr();

  MBReadUtilIface *readUtilIface;
  MBInterface *mdbImpl;
  FILE *cubFile;
  FileTOC fileTOC;
  ModelEntry *modelEntries;
  MetaDataContainer modelMetaData;
  int currNodeIdOffset;
  int currElementIdOffset[MBMAXTYPE];
  MBTag globalIdTag, cubIdTag, geomTag, uniqueIdTag, blockTag, nsTag, ssTag,
    attribVectorTag, entityNameTag, categoryTag;
  std::map<int, MBEntityHandle> uidSetMap;
  std::map<int, MBEntityHandle> gidSetMap[6];

  std::vector<int> int_buf;
  std::vector<double> dbl_buf;
  std::vector<char> char_buf;

  static MBReaderIface* factory( MBInterface* );
  
    // read cub file
  MBErrorCode load_file(const char *file_name,
                        MBEntityHandle& file_set,
                        const FileOptions& opts,
                        const int* block_list,
                        int num_blocks );
  MBErrorCode read_nodeset(ModelEntry *model,
                    NodesetHeader *nodeseth);
  MBErrorCode read_sideset(const double data_version,
                    ModelEntry *model,
                    SidesetHeader *sideseth);
  MBErrorCode read_block(const double data_version,
                  ModelEntry *model,
                  BlockHeader *blockh);
  MBErrorCode read_group(const int gr_index,
                         ModelEntry *model,
                         GroupHeader *grouph);
  MBErrorCode read_nodes(const int gindex,
                  ModelEntry *model,
                  GeomHeader *entity);
  MBErrorCode read_elements(ModelEntry *model,
                     GeomHeader *entity);
  MBErrorCode read_file_header();
  MBErrorCode read_model_entries();
  int find_model(const int model_type);
  MBErrorCode read_meta_data(const int metadata_offset, 
                      MetaDataContainer &mc);
  MBErrorCode read_md_string(std::string &name);
  
  enum {mesh, acist, acisb, facet, exodusmesh};
  MBEntityType type_from_cub_type(const int cub_type, const int nodes_per_elem);
  void check_contiguous(const int num_ents, int &contig, int &max_id);

  Tqdcfr(MBInterface *impl);

  static std::string BLOCK_NODESET_OFFSET_TAG_NAME;
  static std::string BLOCK_SIDESET_OFFSET_TAG_NAME;
  
  MBErrorCode create_set( MBEntityHandle& h, int flags = MESHSET_SET );

private:

  MBEntityHandle mFileSet; // set containing read entities.

  MBErrorCode convert_nodesets_sidesets();

  MBErrorCode read_acis_records( const char* sat_file_name = 0 );
  
  MBErrorCode parse_acis_attribs(const int entity_rec_num,
                          std::vector<AcisRecord> &records);
  MBErrorCode interpret_acis_records(std::vector<AcisRecord> &records);

  MBErrorCode reset_record(AcisRecord &this_record);
  
  MBErrorCode process_record(AcisRecord &this_record);
  
  static const char geom_categories[][CATEGORY_TAG_SIZE];
  
  FILE* acisDumpFile;

    // map between cub ids and MOAB handles
  std::vector<MBEntityHandle> *cubMOABVertexMap;

    // enum used to identify element/entity type in groups
  enum {GROUP = 0, BODY, VOLUME, SURFACE, CURVE, VERTEX, HEX, TET, PYRAMID, QUAD, TRI, EDGE, NODE};
  static const MBEntityType group_type_to_mb_type[];

  enum {SPHERE_EXO=0,
        BAR, BAR2, BAR3,
        BEAM, BEAM2, BEAM3,
        TRUSS, TRUSS2, TRUSS3,
        SPRING,
        TRIthree, TRI3, TRI6, TRI7,
        TRISHELL, TRISHELL3, TRISHELL6, TRISHELL7,
        SHEL, SHELL4, SHELL8, SHELL9,
        QUADfour, QUAD4, QUAD5, QUAD8, QUAD9,
        TETRAfour, TETRA4, TETRA8, TETRA10, TETRA14,
        PYRAMIDfive, PYRAMID5, PYRAMID8, PYRAMID13, PYRAMID18,
        HEXeight, HEX8, HEX9, HEX20, HEX27, HEXSHELL, 
        INVALID_ELEMENT_TYPE};
  static const MBEntityType block_type_to_mb_type[];
  static const int cub_elem_num_verts[];

    //! mapping from mesh packet type to moab type
  static const MBEntityType mp_type_to_mb_type[];
  
    //! get entities with individually-specified types; if is_group is false, 
    //! increment each mem_type by 2 since they're CSOEntityType's and not group types
  MBErrorCode get_entities(const int *mem_types,
                           int *id_buf, const int id_buf_size,
                           const bool is_group,
                           std::vector<MBEntityHandle> &entities);
  
    //! get entities specified by type and ids, append to entities
  MBErrorCode get_entities(const int this_type, 
                           int *id_buf, const int id_buf_size,
                           std::vector<MBEntityHandle> &entities,
                           std::vector<MBEntityHandle> &excl_entities);
  
    //! get ref entity sets with specified type and ids
  MBErrorCode get_ref_entities(const int this_type, 
                               int *id_buf, const int id_buf_size,
                               std::vector<MBEntityHandle> &entities);
  
    //! get mesh entities with specified type and ids
  MBErrorCode get_mesh_entities(const int this_type, 
                                int *id_buf, const int id_buf_size,
                                std::vector<MBEntityHandle> &entities,
                                std::vector<MBEntityHandle> &excl_entities);
  
    //! process entities in a sideset according to sense flags stored in int_buf
    //! or char_buf (depending on sense_size)
  MBErrorCode process_sideset_10(const int this_type, const int num_ents,
                                 const int sense_size,
                                 std::vector<MBEntityHandle> &ss_entities,
                                 Tqdcfr::SidesetHeader *sideseth);

  MBErrorCode process_sideset_11(std::vector<MBEntityHandle> &ss_entities,
                                 std::vector<int> &wrt_ents,
                                 Tqdcfr::SidesetHeader *sideseth);
  
    // put entities into the specfied set, and excluded entities into a 
    // std::vector pointed to by the "Exclude_Entities" tag on that set
  MBErrorCode put_into_set(MBEntityHandle set_handle,
                           std::vector<MBEntityHandle> &entities,
                           std::vector<MBEntityHandle> &excl_entities);
  
};

#endif
