/*
 * Tim's Quick 'N Dirty Cub File Reader (Tqdcfr)
 *
 */

#ifndef TQDCFR
#define TQDCFR

#include "MBInternals.hpp"
#include "MBReaderIface.hpp"

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
/*
#define CHECK_SIZE(buf_vec, new_size)  \
    if (buf_vec.size() < (unsigned int) new_size) buf_vec.resize(new_size)
#define FSEEK(offset) assert(0 == fseek(Tqdcfr::instance()->cubFile, offset, SEEK_SET))
#define FREADI(num_ents) \
   {CHECK_SIZE(Tqdcfr::instance()->int_buf, num_ents); \
   assert(fread(&Tqdcfr::instance()->int_buf[0], 4, num_ents, Tqdcfr::instance()->cubFile) == (unsigned int) num_ents);}
#define FREADD(num_ents) \
   {CHECK_SIZE(dbl_buf, num_ents); \
    assert(fread(&dbl_buf[0], 8, num_ents, Tqdcfr::instance()->cubFile) == (unsigned int) num_ents);}
#define FREADC(num_ents) \
   {CHECK_SIZE(char_buf, num_ents); \
    assert(fread(&char_buf[0], 1, num_ents, Tqdcfr::instance()->cubFile) == (unsigned int) num_ents);}
#define FREADIA(num_ents,array) \
   {assert(fread(array, 4, num_ents, Tqdcfr::instance()->cubFile) == (unsigned int) num_ents);}
#define FREADDA(num_ents,array) \
   {assert(fread(array, 8, num_ents, Tqdcfr::instance()->cubFile) == (unsigned int) num_ents);}
#define FREADCA(num_ents,array) \
   {assert(fread(array, 1, num_ents, Tqdcfr::instance()->cubFile) == (unsigned int) num_ents);}
*/
#define CHECK_SIZE(buf_vec, new_size)  \
    if (buf_vec.size() < (unsigned int) new_size) buf_vec.resize(new_size)
#define FSEEK(offset) assert(0 == fseek(instance->cubFile, offset, SEEK_SET))
#define FREADI(num_ents) \
   {CHECK_SIZE(instance->int_buf, num_ents); \
   assert(fread(&instance->int_buf[0], 4, num_ents, instance->cubFile) == (unsigned int) num_ents);}
#define FREADD(num_ents) \
   {CHECK_SIZE(dbl_buf, num_ents); \
    assert(fread(&dbl_buf[0], 8, num_ents, instance->cubFile) == (unsigned int) num_ents);}
#define FREADC(num_ents) \
   {CHECK_SIZE(char_buf, num_ents); \
    assert(fread(&char_buf[0], 1, num_ents, instance->cubFile) == (unsigned int) num_ents);}
#define FREADIA(num_ents,array) \
   {assert(fread(array, 4, num_ents, instance->cubFile) == (unsigned int) num_ents);}
#define FREADDA(num_ents,array) \
   {assert(fread(array, 8, num_ents, instance->cubFile) == (unsigned int) num_ents);}
#define FREADCA(num_ents,array) \
   {assert(fread(array, 1, num_ents, instance->cubFile) == (unsigned int) num_ents);}


class MBReadUtilIface;

class Tqdcfr : public MBReaderIface
{
public:  
  class FileTOC
  {
  public:
    int fileEndian, fileSchema, numModels, modelTableOffset, 
      modelMetaDataOffset, activeFEModel;

    FileTOC() 
        : fileEndian(0), fileSchema(0), numModels(0), modelTableOffset(0), 
          modelMetaDataOffset(0), activeFEModel(0) {}

    void print() 
      {
        std::cout << "FileTOC:End, Sch, #Mdl, TabOff, "
                  << "MdlMDOff, actFEMdl = ";
        std::cout << fileEndian << ", " << fileSchema << ", " << numModels 
                  << ", " << modelTableOffset << ", " 
                  << modelMetaDataOffset << ", " << activeFEModel << std::endl;
      }
  };

  class FEModelHeader 
  {
  public:
    int feEndian, feSchema, feCompressFlag, feLength;

    class ArrayInfo 
    {
    public:
      int numEntities, tableOffset, metaDataOffset;

      ArrayInfo() : numEntities(0), tableOffset(0), metaDataOffset(0) 
        {}
      
      void print() 
        {
          std::cout << "ArrayInfo:numEntities, tableOffset, metaDataOffset = "
                    << numEntities << ", " << tableOffset << ", " << metaDataOffset << std::endl;
        }
      void init( const std::vector<int>& int_buf )
        {
          numEntities = int_buf[0]; tableOffset = int_buf[1]; metaDataOffset = int_buf[2];
        }
    };
    
    ArrayInfo geomArray, nodeArray, elementArray, groupArray, 
      blockArray, nodesetArray, sidesetArray;

    void init(const int offset, Tqdcfr* instance );
        
    void print() 
      {
        std::cout << "FEModelHeader:feEndian, feSchema, feCompressFlag, feLength = "
                  << feEndian << ", " << feSchema << ", " << feCompressFlag << ", " << feLength << std::endl;
        
        std::cout << "geomArray: "; geomArray.print();
        std::cout << "nodeArray: "; nodeArray.print();
        std::cout << "elementArray: "; elementArray.print();
        std::cout << "groupArray: "; groupArray.print();
        std::cout << "blockArray: "; blockArray.print();
        std::cout << "nodesetArray: "; nodesetArray.print();
        std::cout << "sidesetArray: "; sidesetArray.print();
      }
  };

  class EntityHeader 
  {
  public:
    int entityID, entityType, numNodes, nodeOffset, numElements, elementOffset,
      numTypes, entityColor, entityLength,
      pyrType, matType, blockDimension, shellsFlag;

    MBEntityHandle setHandle;

    void print() 
      {
        std::cout << "EH:uid, Typ, #N, NO, #E, EO, #Typ, Col, EntL, "
                  << "pyrTyp, matTyp, Dim, shlFlag = " << std::endl;
        std::cout << entityID << ", " << entityType << ", " << numNodes << ", " 
                  << nodeOffset << ", " << numElements 
                  << ", " << elementOffset << ", " << numTypes << ", " << entityColor 
                  << ", " << entityLength << ", " << pyrType << ", " << matType << ", "
                  << blockDimension << ", " << shellsFlag << std::endl;
      }

    enum {geom, group, block, nodeset, sideset};
    
    static void read_info_header(const int model_offset, 
                                 const FEModelHeader::ArrayInfo &info,
                                 const int info_type,
                                 Tqdcfr* instance,
                                 EntityHeader *&entity_headers);

    EntityHeader() 
        : entityID(0), entityType(0), numNodes(0), nodeOffset(0), numElements(0), 
          elementOffset(0), numTypes(0), entityColor(0), entityLength(0),
          pyrType(0), matType(0), blockDimension(0), shellsFlag(0),
          setHandle(0)
      {}
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
      
      MetaDataEntry() : mdOwner(0), mdDataType(0), mdIntValue(0), 
                        mdName("(uninit)"), mdStringValue("(uninit)"), mdDblValue(0) 
        {}

      void print() 
        {
          std::cout << "MetaDataEntry:own, typ, name, I, D, S = "
                    << mdOwner << ", " << mdName << ", " << mdIntValue << ", " 
                    << mdDblValue << ", " << mdStringValue;
          unsigned int i;
          if (mdIntArrayValue.size()) {
            std::cout << std::endl << "IArray = " << mdIntArrayValue[0];
            for (i = 1; i < mdIntArrayValue.size(); i++)
              std::cout << ", " << mdIntArrayValue[i];
          }
          if (mdDblArrayValue.size()) {
            std::cout << std::endl << "DArray = " << mdDblArrayValue[0];
            for (i = 1; i < mdDblArrayValue.size(); i++)
              std::cout << ", " << mdDblArrayValue[i];
          }
          std::cout << std::endl;
        }
    };

    void print() 
      {
        std::cout << "MetaDataContainer:mdSchema, compressFlag, numDatums = "
                  << mdSchema << ", " << compressFlag << ", " << numDatums << std::endl;

        for (int i = 0; i < numDatums; i++)
          metadataEntries[i].print();
      }

    MetaDataEntry *metadataEntries;
    MetaDataContainer() : mdSchema(0), compressFlag(0), numDatums(0), metadataEntries(NULL) {}
    ~MetaDataContainer() {if (NULL != metadataEntries) delete [] metadataEntries;}
  };

  class ModelEntry
  {
  public:
    ModelEntry() 
        : modelHandle(0), modelOffset(0), 
          modelLength(0), modelType(0), modelOwner(0), modelPad(0),
          feGeomH(NULL), feGroupH(NULL), feBlockH(NULL), 
          feNodeSetH(NULL), feSideSetH(NULL)
      {}

    ~ModelEntry() 
      {
        delete [] feGeomH; delete [] feGroupH; delete [] feBlockH;
        delete [] feNodeSetH; delete [] feSideSetH;
      }
    
    int modelHandle, modelOffset, modelLength, modelType, modelOwner, modelPad;

    FEModelHeader feModelHeader;
    EntityHeader *feGeomH, *feGroupH, *feBlockH, *feNodeSetH, *feSideSetH;
    MetaDataContainer geomMD, nodeMD, elementMD, groupMD, blockMD, nodesetMD, sidesetMD;
    
    void print() 
      {
        std::cout << "ModelEntry: Han, Of, Len, Tp, Own, Pd = "
                  << modelHandle << ", " << modelOffset << ", " << modelLength 
                  << ", " << modelType << ", " << modelOwner << ", " << modelPad
                  << std::endl;
      }

    void print_header(const char *prefix,
                      FEModelHeader::ArrayInfo &info,
                      EntityHeader *header);
    
    void read_header_info( Tqdcfr* instance );
    void read_metadata_info(Tqdcfr *tqd);
  };

  enum {BODY, LUMP, SHELL, FACE, LOOP, COEDGE, EDGE, VERTEX, ATTRIB, UNKNOWN};
  
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
  MBTag globalIdTag, geomTag, uniqueIdTag, groupTag, blockTag, nsTag, ssTag,
    attribVectorTag, entityNameTag;

  std::vector<int> int_buf;
  std::vector<double> dbl_buf;
  std::vector<char> char_buf;

  static MBReaderIface* factory( MBInterface* );
  
    // read cub file
  MBErrorCode load_file(const char *file_name,
                        const int* block_list,
                        int num_blocks );
  void read_nodeset(ModelEntry *model,
                    EntityHeader *nodeseth);
  void read_sideset(ModelEntry *model,
                    EntityHeader *sideseth);
  void read_block(ModelEntry *model,
                  EntityHeader *blockh);
  void read_group(ModelEntry *model,
                  EntityHeader *grouph);
  void read_nodes(ModelEntry *model,
                  EntityHeader *entity);
  void read_elements(ModelEntry *model,
                     EntityHeader *entity);
  void read_file_header();
  void read_model_entries();
  int find_model(const int model_type);
  void read_meta_data(const int metadata_offset, 
                      MetaDataContainer &mc);
  void read_md_string(std::string &name);
  
  enum {mesh, acist, acisb, facet, exodusmesh};
  MBEntityType type_from_cub_type(const int cub_type, const int nodes_per_elem);
  void add_set_entities(const int this_type, const int num_ents,
                        const int *ints, MBEntityHandle &set_handle);
  int check_contiguous(const int num_ents);

  Tqdcfr(MBInterface *impl);
  
private:

  void read_acis_records();
  
  void parse_acis_attribs(const int entity_rec_num,
                          std::vector<AcisRecord> &records);
  void interpret_acis_records(std::vector<AcisRecord> &records);

  void reset_record(AcisRecord &this_record);
  
  void process_record(AcisRecord &this_record);
  
  Tqdcfr *const instance;
};

#endif
