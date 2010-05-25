#ifndef READCCMIO_HPP
#define READCCMIO_HPP

#ifndef IS_BUILDING_MB
#error "ReadCCMIO.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <set>
#include <map>
#include <string>

#include "moab/Forward.hpp"
#include "moab/Range.hpp"
#include "moab/ExoIIInterface.hpp"
#include "moab/ReaderIface.hpp"
#include "moab/TupleList.hpp"
#include "ccmio.h"

namespace moab {

#undef READCCMIO_USE_TUPLE_LIST

class ReadUtilIface;

class ReadCCMIO : public ReaderIface
{
 
public:

  typedef std::map<int, std::vector<EntityHandle> > TupleList;
  typedef std::map<int, std::vector<int> > SenseList;

    //! Constructor
  ReadCCMIO(Interface *impl);

   //! Destructor
  virtual ~ReadCCMIO();
  
  static ReaderIface* factory( Interface* );

  ErrorCode load_file(const char *file_name,
                        const EntityHandle* file_set,
                        const FileOptions& opts,
                        const ReaderIface::IDTag* subset_list,
                        int subset_list_length,
                        const Tag* file_id_tag);

private:
  
  ErrorCode read_processor(CCMIOID rootID, CCMIOID problemID,
                           CCMIOID processorID, CCMIOID verticesID,
                           CCMIOID topologyID, CCMIOSize_t proc,
                             Range *new_ents);

  ErrorCode read_cells(CCMIOSize_t proc, CCMIOID processorID,
                         CCMIOID verticesID, CCMIOID topologyID,
                         TupleList &vert_map, Range *new_cells);


  ErrorCode construct_cells(TupleList &face_map, 
#ifndef READCCMIO_USE_TUPLE_LIST
                            SenseList &sense_map,
#endif
                            TupleList &vert_map, 
                            std::map<int,int> &cell_topo_types,
                            std::vector<EntityHandle> &new_cells);

  ErrorCode ccmio_to_moab_type(int ccm_type, EntityType &moab_type, bool &has_mid_nodes);

  ErrorCode create_cell_from_faces(std::vector<EntityHandle> &facehs,
                                     std::vector<int> &senses,
                                   EntityType this_type,
                                   bool has_mid_nodes,
                                   EntityHandle &cell);

  ErrorCode read_gids_and_types(CCMIOID problemID,
                                  CCMIOID topologyID,
                                  std::vector<EntityHandle> &cells);

  ErrorCode read_all_faces(CCMIOID topologyID, TupleList &vert_map, 
                             TupleList &face_map
#ifndef READCCMIO_USE_TUPLE_LIST
                             ,SenseList &sense_map
#endif
                             , Range *new_faces);


  ErrorCode read_faces(CCMIOID faceID, 
                       CCMIOEntity bdy_or_int,
                         TupleList &vert_map,
                         TupleList &face_map
#ifndef READCCMIO_USE_TUPLE_LIST
                         ,SenseList &sense_map
#endif
                         , Range *new_faces);

  ErrorCode make_faces(int *farray, 
                       TupleList &vert_map,
                       std::vector<EntityHandle> &new_faces,
                       int num_faces);

  ErrorCode read_vertices(CCMIOSize_t proc, CCMIOID processorID, CCMIOID verticesID,
                            CCMIOID topologyID, 
                            Range *verts, TupleList &vert_map);


  ErrorCode get_processors(CCMIOID stateID, CCMIOID &processorID,
                           CCMIOID &verticesID, CCMIOID &topologyID,
                           CCMIOID &solutionID, 
                           std::vector<CCMIOSize_t> &procs, 
                           bool &has_solution);

  ErrorCode get_state(CCMIOID rootID, CCMIOID &problemID, CCMIOID &stateID);


  virtual ErrorCode read_tag_values( const char* file_name,
                                       const char* tag_name,
                                       const FileOptions& opts,
                                       std::vector<int>& tag_values_out,
                                       const IDTag* subset_list = 0,
                                       int subset_list_length = 0 );
  
  ErrorCode load_matset_data(CCMIOID problemID);
  
  ErrorCode load_neuset_data(CCMIOID problemID);
  
  ErrorCode load_metadata(CCMIOID rootID, CCMIOID problemID, 
                          CCMIOID stateID, CCMIOID processorID,
                            const EntityHandle *file_set);
  
  ErrorCode create_matset_tags(Tag &matNameTag, Tag &matPorosityTag, 
                                 Tag &matSpinTag, Tag &matGroupTag);

  ErrorCode read_topology_types(CCMIOID &topologyID, 
                                std::map<int,int> &cell_topo_types);

  ErrorCode get_opt_string(const char *opt_name, CCMIOID node, std::vector<char> &opt_string);
  
    //! Cached tags for reading.  Note that all these tags are defined when the
    //! core is initialized.
  Tag mMaterialSetTag;
  Tag mDirichletSetTag;
  Tag mNeumannSetTag;
  Tag mHasMidNodesTag;
  Tag mGlobalIdTag;
  Tag mNameTag;
  
  Interface *mbImpl;

  ReadUtilIface* readMeshIface;

  std::map<int,EntityHandle> newMatsets, newNeusets;
  
  bool hasSolution;
};

} // namespace moab

#endif
