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


#ifndef MB_WRITE_UTIL_HPP
#define MB_WRITE_UTIL_HPP

#ifndef IS_BUILDING_MB
#error "MBWriteUtil.hpp isn't supposed to be included into an application"
#endif

#include "MBWriteUtilIface.hpp"

class MBCore;
class MBError;

class MBWriteUtil : public MBWriteUtilIface
{
private:
  //! pointer to the MBCore
  MBCore* mMB;
  MBError* mError;
public:

  //! constructor takes MBCore pointer
  MBWriteUtil(MBCore* mdb, MBError* error_handler);

  //! destructor
  ~MBWriteUtil(){}
  
    //! Check if the specified file already exists.
    //! Returns MB_SUCCESS if file does not exist, MB_ALREADY_ALLOCATED
    //! if file does exist, or MB_FAILURE for some other error condition.
  virtual MBErrorCode check_doesnt_exist( const char* file_name );

  //! gets arrays for coordinate data from the MB
  MBErrorCode get_node_arrays(
      const int num_arrays,
      const int num_nodes, 
      const MBRange& entities,
      MBTag node_id_tag,
      const int start_node_id,
      std::vector<double*>& arrays
      );
      
  /** Get an array of coordinate values for nodes
   *
   * Given a range of node handles, retreive a single coordinate
   * value for each. 
   *
   * Failure conditions:
   *  - invalid entity handles (not vertices, non-existant entity, etc.)
   *  - range is empty (<code>iter == end</code>)
   *  - <code>output_array</code> is null
   *  - insufficient space in <code>output_array</code>
   *
   *\param which_array  The coordinate to retreive (0-&gt;X, 1-&gt;Y, 2-&gt;Z)
   *\param begin        The first node handle.
   *\param end          One past the last node handle.
   *\param output_size  The size of <code>output_array</code>.
   *\param output_array The memory in which to write the node coordinates.
   *\author Jason Kraftcheck
   */
  MBErrorCode get_node_array(
      const int which_array, 
      MBRange::const_iterator begin,
      const MBRange::const_iterator end,
      const size_t output_size,
      double* const output_array
      );

  //! get array for connectivity data from the MB
  MBErrorCode get_element_array(
      const int num_elements, 
      const int verts_per_element,
      MBTag node_id_tag,
      const MBRange& entities, 
      MBTag element_id_tag,
      int start_element_id,
      int* array
      );

  /** Get connectivity for elements 
   *
   * Get the connectivity list for a range of elements.
   *
   * Failure cases:
   *  - Passed range is empty (<code>begin == end</code>).
   *  - <code>vertices_per_elem</code> is less than one
   *  - <code>element_array</code> is null.
   *  - The range contains invalid handles (non-existant entities,
   *      not an element, etc.)
   *  - Retreiving ID tag for an entity failed.
   *  - Insufficient space in passed array.
   *
   *\param begin        The first element handle
   *\param end          One past the last element handle
   *\param vertices_per_elem Number of vertices to retreive for each
   *                    element.  If the element has more vertices, the
   *                    element connectivity will be truncated.  If 
   *                    <code>vertices_per_elem</code> is greater than the
   *                    number of nodes for an eleement, the data will be
   *                    padded with zeros.
   *\param node_id_tag  A tag with integer values.  
   *\param array_size   The length of <code>element_array</code>
   *\param element_array The memory location at which to store the 
   *                    connectivity list.
   *\author Jason Kraftcheck
   */
  MBErrorCode get_element_array(
      MBRange::const_iterator begin,
      const MBRange::const_iterator end,
      const int vertices_per_elem,
      MBTag node_id_tag,
      const size_t array_size, 
      int *const element_array
      );

  /** Get connectivity for elements 
   *
   * Get the connectivity list for a range of elements.
   *
   * Failure cases:
   *  - Passed range is empty (<code>begin == end</code>).
   *  - <code>vertices_per_elem</code> is less than one
   *  - <code>element_array</code> is null.
   *  - The range contains invalid handles (non-existant entities,
   *      not an element, etc.)
   *  - Insufficient space in passed array.
   *
   *\param begin        The first element handle
   *\param end          One past the last element handle
   *\param vertices_per_elem Number of vertices to retreive for each
   *                    element.  If the element has more vertices, the
   *                    element connectivity will be truncated.  If 
   *                    <code>vertices_per_elem</code> is greater than the
   *                    number of nodes for an eleement, the data will be
   *                    padded with zeros.
   *\param array_size   The length of <code>element_array</code>
   *\param element_array The memory location at which to store the 
   *                    connectivity list.
   *\author Jason Kraftcheck
   */
  virtual MBErrorCode get_element_array(
      MBRange::const_iterator begin,
      const MBRange::const_iterator end,
      const int vertices_per_elem,
      const size_t array_size, 
      MBEntityHandle *const element_array
      );

  /** Get poly (polygon or polyhedron) connectivity size
   *\param begin  First iterator in range of poly
   *\param end    One past last in range of poly.
   *\param connectivity_size  The lenght of the connectivity list
   *              For the specified range of polyhedra.
   *\author Jason Kraftcheck
   */
  virtual MBErrorCode get_poly_array_size(
      MBRange::const_iterator begin,
      const MBRange::const_iterator end,
      int& connectivity_size 
      );
   
   
  /** Get poly (polygon or polyhedron) connectivity.
   *
   * This function will add as many polys as possible to the
   * passed arrays given the sizes of those arrays.  It will
   * then pass back position at which it stoped and the sizes
   * of the data written to the arrays.
   *
   *\param iter               As input, the first element handle.
   *                          As output, one past the last element handle
   *                          for which data was written to the arrays.
   *\param end                The iterator at which to stop.
   *\param node_id_tag        A tag with integer values.  
   *\param element_array_len  As input, length of <code>element_array</code>.
   *                          As output, the number of entries written in that
   *                          array.
   *\param element_array      The memory location at which to store the 
   *                          connectivity list.
   *\param index_array_len    As input, the length of <code>index_array</code>.
   *                          As output, the number of entries written in that
   *                          array.
   *\param index_array        The memory location at which to store offsets.
   *\param index_offset       Value to offset (add to) index values.  As output
   *                          the input value plus the amount of data 
   *                          written to the element array.  (The value you
   *                          presumably want to pass to the next call.)
   *\author Jason Kraftcheck
   */
  virtual MBErrorCode get_poly_arrays(
      MBRange::const_iterator& iter,
      const MBRange::const_iterator end,
      const MBTag node_id_tag,
      size_t& handle_array_len,
      int *const handle_array,
      size_t& index_array_len,
      int *const index_array,
      int& index_offset 
      );
 
  //! get a set of nodes that represent a set of elements
  MBErrorCode gather_nodes_from_elements(
      const MBRange& elements,
      const MBTag node_bit_mark_tag,
      MBRange& nodes
      );
  
    //! assign ids to input elements starting with start_id, written to id_tag
    //! if zero, assigns to GLOBAL_ID_TAG_NAME
  MBErrorCode assign_ids(MBRange &elements,
                         MBTag id_tag,
                         const int start_id);

  
  /** Get explict adjacencies 
   *
   * Get explicit adjacences stored in database.
   * Does not create any explicit adjacencies or search for
   * implicit ones.
   *
   *\param entity  The entity to retreive adjacencies for.
   *\param id_tag  The global ID tag
   *\param adj     The output list of global IDs of adjacent entities.
   */
  MBErrorCode get_adjacencies(
      MBEntityHandle entity,
      MBTag id_tag,
      std::vector<int>& adj 
  );
  
  MBErrorCode get_adjacencies( MBEntityHandle entity,
                               const MBEntityHandle*& adj_array,
                               int& num_adj );


  /**\brief Get list of tags to write.
   *
   * Get the list of tags to write to the file, possibly using
   * an optional user-specifed tag list.
   *
   *\author Jason Kraftcheck
   */
  virtual MBErrorCode 
  get_tag_list( std::vector<MBTag>& result_list,
                const MBTag* user_tag_list = 0, 
                int user_tag_list_length = 0,
                bool include_variable_length_tags = true );

  //! tell MB there was an error when writing the mesh
  //! it makes sense to have this as long as MBInterface has a write_mesh function
  MBErrorCode report_error( const std::string& error );

  MBErrorCode report_error( const char* error, ... )
#ifdef __GNUC__
 __attribute__((format(printf,2,3)))
#endif
  ;
};

#endif


