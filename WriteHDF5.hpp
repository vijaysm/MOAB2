/**
 * \class  WriteHDF5
 * \brief  Write mesh database to TSTT HDF5 file.
 * \author Jason Kraftcheck
 * \date   01 April 2004
 */

#ifndef WRITE_HDF5_HPP
#define WRITE_HDF5_HPP

#include <list>
#include "mhdf.h"
#include "MBInterface.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBRange.hpp"
#include "MBWriterIface.hpp"

class MB_DLL_EXPORT WriteHDF5 : public MBWriterIface
{
public:

  static MBWriterIface* factory( MBInterface* );

  /** The type of the global ID tag data 
   * 
   * NOTE:  If this is changed, the value of id_type 
   *        MUST be changed accordingly.
   */
  typedef int id_t;
  
  /** HDF5 type corresponding to type of id_t */
  static const hid_t id_type;

  WriteHDF5( MBInterface* iface );
  
  virtual ~WriteHDF5();
  
  /** Export specified meshsets to file
   * \param filename     The filename to export.  Must end in <em>.mhdf</em>
   * \param export_sets  Array of handles to sets to export, or NULL to export all.
   * \param export_set_count Length of <code>export_sets</code> array.
   */
  MBErrorCode write_file( const char* filename,
                          bool overwrite,
                          const MBEntityHandle* export_sets,
                          const int export_set_count,
                          std::vector<std::string>& qa_records,
                          int user_dimension = 3 );

  /** Create attributes holding the HDF5 type handle for the 
   *  type of a bunch of the default tags.
   */
  static MBErrorCode register_known_tag_types( MBInterface* );
  

private:
  MBErrorCode init();
  
  //! Range of entities, grouped by type, to export 
  struct ExportSet 
  {
    //! The range of entities.
    MBRange range;
    //! The type of the entities in the range
    MBEntityType type;
    //! The number of nodes per entity - not used for nodes and sets
    int num_nodes;
    //! The type handle for the mhdf library.
    mhdf_ElemHandle type2;
    //! The first Id allocated by the mhdf library.  Entities in range have sequential IDs.
    id_t first_id;
  };
  
  //! The size of the data buffer (<code>dataBuffer</code>).
  const int bufferSize;
  //! A memory buffer to use for all I/O operations.
  char* dataBuffer;

  //! MBInterface pointer passed to constructor
  MBInterface* iFace;
  
  //! Cached pointer to writeUtil interface.
  MBWriteUtilIface* writeUtil;
  
  //! The file handle from the mhdf library
  mhdf_FileHandle filePtr;
  
  //! True if created the ID tag in init()
  bool createdIdTag;
  
  //! Handle for the ID tag.
  MBTag idTag;
  
  //! The list elements to export.
  std::list<ExportSet> exportList;
  
  //! The list of nodes to export
  ExportSet nodeSet;
  
  //! The list of sets to export
  ExportSet setSet;

  //! Zero the ID tag on all entities in the mesh.
  MBErrorCode clear_all_id_tags();
  
  /** Get the subset of an entity range that is a specified element type
   * \param input_range  The range to subset
   * \param type         The base element type
   * \param number_nodes The number of nodes per element
   * \param output_range The result subset of the input range.
   */
  MBErrorCode subrange_by_type_and_conn( const MBRange& input_range,
                                         MBEntityType type,
                                         int number_nodes,
                                         MBRange& output_range );
  
  /** Get the subrange of a given range containing the specified elem type
   * \param input_range  The range to subset
   * \param type         The element type
   * \param output_range The subset of input_range of type <code>type</code>
   */
  MBErrorCode subrange_by_type( const MBRange& input_range,
                                MBEntityType type,
                                MBRange& output_range );
  
  /** Get higher-order types
   *
   * For each higher-order type of the element, append the number of
   * of nodes beyond those of the base/simplest type that that higher-
   * order type has.  Also appends zero to the list for the base type.
   */
  MBErrorCode midnode_combinations( MBEntityType type,
                                    std::vector<int>& combinations );
  
  /** Get all mesh to export from given list of sets.
   *
   * Populate exportSets, nodeSet and setSet with lists of
   * entities to write.
   *
   * \param export_sets  The list of meshsets to export
   */
  MBErrorCode gather_mesh_info( const std::vector<MBEntityHandle>& export_sets );
  
  //! Same as gather_mesh_info, except for entire mesh
  MBErrorCode gather_all_mesh( );
 
  /** Write out the nodes.
   *
   * Note: Assigns IDs to nodes.
   */
  MBErrorCode write_nodes( int dimension );
  
  /** Write out element connectivity.
   *
   * Write connectivity for passed set of elements.
   *
   * Note: Assigns element IDs.
   * Note: Must do write_nodes first so node IDs get assigned.
   */
  MBErrorCode write_elems( ExportSet& elemset );
  
  /** Write out poly(hedr/g)on connectivity
   *
   * Write connectivity for passed set of poly
   *
   * Note: Assigns element IDs.
   * Note: Must do write of lower-dimension entities first so
   *       IDs get assigned to them.
   */
  MBErrorCode write_poly( ExportSet& elemset );
  
  /** Write out meshsets
   * 
   * Write passed set of meshsets, including parent/child relations.
   *
   * Note: Must have written nodes and element connectivity
   *       so entities have assigned IDs.
   */
  MBErrorCode write_sets( );
  
  /** Write adjacency info for passed set of elements
   *
   * Note: Must have written element connectivity so elements
   *       have IDs assigned.
   */
  MBErrorCode write_adjacencies( const ExportSet& export_set );
  
  /** Write tag information and data.
   * 
   * Note: Must have already written nodes, elem connectivity and
   *       sets so that entities have IDs assigned.
   */
  MBErrorCode write_tag( MBTag tag_handle );
  
  /** Get possibly compacted list of IDs for passed entities
   *
   * For the passed range of entities, determine if IDs
   * can be compacted and optionally write IDs to passed list.
   *
   * If the IDs are not compacted, the output list will contain
   * a simple ordered list of IDs.
   *
   * If IDs are compacted, the output list will contain 
   * {start,count} pairs.
   *
   * If the ID list is compacted, output_length will be less than
   * range.size().
   */
  MBErrorCode range_to_id_list( const MBRange& input_range,
                                std::vector<id_t>* output_id_list,
                                id_t& output_length );
  
  //! Get IDs for entities 
  MBErrorCode vector_to_id_list( const std::vector<MBEntityHandle>& input,
                                 std::vector<id_t>& output );

  /** Get IDs of adjacent entities.
   * 
   * For all entities adjacent to the passed entity, if the
   * adjacent entity is to be exported (ID is not zero), append
   * the ID to the passed list.
   */
  MBErrorCode get_adjacencies( MBEntityHandle entity, std::vector<id_t>& adj );
  
  //! Write dense tag for all entities 
  MBErrorCode write_dense_tag( MBTag tag_handle,
                               hid_t hdf_write_type );

  //! Write dense tag for specified entity set
  MBErrorCode write_dense_tag( ExportSet& set,
                               MBTag tag_handle,
                               hid_t hdf_write_type );
  
  //! Write sparse tag for all entities.
  MBErrorCode write_sparse_tag( MBTag tag_handle,
                                hid_t hdf_write_type );

  //! Get information about a meshset
  MBErrorCode get_set_info( MBEntityHandle set,
                            long& num_entities,
                            long& num_children,
                            unsigned long& flags );
                            
  //! Get element connectivity
  MBErrorCode get_connectivity( MBRange::const_iterator begin,
                                MBRange::const_iterator end,
                                int nodes_per_element,
                                id_t* id_data_out );

  //! Write exodus-type QA info
  MBErrorCode write_qa( std::vector<std::string>& list );
};

#endif
