//-------------------------------------------------------------------------
// Filename      : ReadNCDF.hpp
//
// Purpose       : ExodusII reader
//
// Special Notes : Lots of code taken from verde implementation
//
// Creator       : Tim Tautges & Corey Ernst
//
// Date          : 3/02
//
// Owner         : Tim Tautges & Corey Ernst
//-------------------------------------------------------------------------

#ifndef READNCDF_HPP
#define READNCDF_HPP

#ifndef IS_BUILDING_MB
#error "ReadNCDF.hpp isn't supposed to be included into an application"
#endif

#include <set>
#include <vector>
#include <deque>
#include <functional>
#include <string>

#include "MBReadUtilIface.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"
#include "ExoIIUtil.hpp"

struct ReadBlockData
{
  int blockId;
  int startExoId; 
  MBEntityHandle startMBId; 
  int numElements;
  bool reading_in;
  ExoIIElementType elemType;
};

class NcFile;

//! Output Exodus File for VERDE
class ReadNCDF 
{
   
public:
    //! load an ExoII file
  MBErrorCode load_file(const char *exodus_file_name,
                         const int* blocks_to_load,
                         const int num_blocks);
  
   //! Constructor
   ReadNCDF(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadNCDF();

private:

  MBReadUtilIface* readMeshIface;

  bool dimension_exists(const char *attrib_name);
  
  void reset();

    //! read the header from the ExoII file
  MBErrorCode read_exodus_header(const char *exodus_file_name);
  
    //! read the nodes
  MBErrorCode read_nodes();
  
    //! read block headers, containing info about element type, number, etc.
  MBErrorCode read_block_headers(const int *blocks_to_load,
                                  const int num_blocks);
  
    //! read the element blocks
  MBErrorCode read_elements();
  
    //! read in the global element ids
  MBErrorCode read_global_ids();

    //! read the nodesets into meshsets
  MBErrorCode read_nodesets();
  
    //! read the sidesets (does nothing for now)
  MBErrorCode read_sidesets();

    //! exodus file bound to this object
  int exodus_file();

    //! number of dimensions in this exo file
  int number_dimensions();

  //! map a character exodusII element type to a TSTT type & topology
  MBErrorCode get_type(char *exo_element_type,
                        MBEntityType &elem_type);
 
  MBErrorCode get_type(MBEntityType &elem_type,
                        std::string &exo_element_type);

  /* 
  int get_int_tag(const MB_MeshSet *this_ms,
                  const TagHandle tag_id);
 */

  //qa record stuff 
  MBErrorCode read_qa_records();
  MBErrorCode read_qa_information( std::vector<char*> &qa_record_list);

  MBErrorCode read_qa_string(char *string,
                              int record_number,
                              int record_position); 

  MBErrorCode create_ss_elements( int *element_ids, int *side_list,
                                   int num_sides, int num_dist_factors,
                                   std::vector<MBEntityHandle> &entities_to_add,
                                   std::vector<MBEntityHandle> &reverse_entities,
                                   std::vector<double> &dist_factor_vector,
                                   int ss_seq_id);

  MBErrorCode find_side_element_type( const int element_id, ExoIIElementType &type, 
                                       ReadBlockData &block_data, int &df_index, int side_id );

  MBErrorCode remove_previously_loaded_blocks(const int *blocks_to_load,
                                               const int num_blocks,
                                               std::vector<int> &new_blocks);
  
 /* MBErrorCode assign_block_ids_to_ssets(MBEntityHandle ss_handle,
                                         MB_MeshSet *ss_mesh_set);
                                         */

  //! creates an element with the given connectivity
  MBErrorCode create_sideset_element( std::vector<MBEntityHandle>, MBEntityType, MBEntityHandle&);

  //! I think this ought to be moved to MBCore. KGM
  MBErrorCode check_file_status(const char *&file_name,
                                 bool& previously_read);

  int get_number_nodes( MBEntityHandle handle );



  //------------member variables ------------//

    //! interface instance
  MBInterface* mdbImpl;
  
  NcFile *ncFile;        // netcdf/exodus file

  int CPU_WORD_SIZE;
  int IO_WORD_SIZE;

    //! int to offset vertex ids with
  int vertexOffset;

    //! file name
  std::string exodusFile;

    //! number of nodes in the current exoII file
  int numberNodes_loading;

    //! number of elements in the current exoII file
  int numberElements_loading;

    //! number of blocks in the current exoII file
  int numberElementBlocks_loading; 

    //! number of nodesets in the current exoII file
  int numberNodeSets_loading; 

    //! number of sidesets in the current exoII file
  int numberSideSets_loading; 

  int numberDimensions_loading;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;

  //keeps track of the exodus ids of the elements and nodes just loaded
  std::vector<char> nodesInLoadedBlocks;
  //note- vector<bool> has limited capabilities

  //vector of blocks that are loading 
  std::vector< ReadBlockData > blocksLoading;

  //! Cached tags for reading.  Note that all these tags are defined when the
  //! core is initialized.
  MBTag mMaterialSetTag;
  MBTag mDirichletSetTag;
  MBTag mNeumannSetTag;
  MBTag mHasMidNodesTag;
  MBTag mDistFactorTag;
  MBTag mGlobalIdTag;
  MBTag mQaRecordTag;

  int max_line_length, max_str_length;
};

// inline functions
inline int ReadNCDF::number_dimensions() 
{
   return numberDimensions_loading;
}


#endif




