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

/**
 * \class ReadGmsh
 * \brief Gmsh (http://www.geuz.org/gmsh) file reader
 * \author Jason Kraftcheck
 */


#ifndef READ_GMSH_HPP
#define READ_GMSH_HPP

#include "MBInterface.hpp"
#include "MBReaderIface.hpp"
#include "MBRange.hpp"

class MBReadUtilIface;

// Base class for binary and ASCII readers
class ReadGmsh : public MBReaderIface
{
   
public:

    //! factory method 
  static MBReaderIface* factory( MBInterface* );

  MBErrorCode load_file(const char *file_name,
                        const int* material_set_list,
                        const int num_material_sets );
  
    //! Constructor
  ReadGmsh(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadGmsh();

  struct ElementType {
    MBEntityType mbtype;
    int nodes;
    const int* node_order;
  };

private:

  MBErrorCode load_file_impl( const char *file_name,
                              const int* material_set_list,
                              const int num_material_sets );
  
  MBErrorCode create_elements( const ElementType& type,
                               const std::vector<int>& elem_ids,
                               const std::vector<int>& matl_ids,
                               const std::vector<int>& geom_ids,
                               const std::vector<int>& prtn_ids,
                               const std::vector<MBEntityHandle>& connectivity );

  MBErrorCode create_sets( MBEntityType element_type,
                           const MBRange& elements,
                           const std::vector<int>& set_ids,
                           int set_type );
  
  MBErrorCode create_geometric_topology();
  

  MBReadUtilIface* readMeshIface;

    //! interface instance
  MBInterface* mdbImpl;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;
  
  MBTag globalId;
  MBRange geomSets;
};


#endif
