/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2008 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#ifndef SWEPT_ELEMENT_SEQUENCE
#define SWEPT_ELEMENT_SEQUENCE

//
// Class: SweptElementSequence
//
// Purpose: represent a swept element of mesh
//

#include "ElementSequence.hpp"
#include "ScdElementData.hpp"

class SweptElementSeq : public ElementSequence
{
public:

    //! constructor
  SweptElementSeq(
		  MBEntityHandle start_handle,
		  const int imin, const int jmin, const int kmin,
		  const int imax, const int jmax, const int kmax,
		  const int* Cq );
  
  virtual ~SweptElementSeq();

    //! given a handle, get the corresponding parameters
  MBErrorCode get_params(const MBEntityHandle ehandle,
                          int &i, int &j, int &k) const
    { }
  
    //! get connectivity of an entity given entity's parameters
  MBErrorCode get_params_connectivity(const int i, const int j, const int k,
                                std::vector<MBEntityHandle>& connectivity) const
  { }
  
  
    /***************** Methods from ElementSeq *****************/

  virtual MBErrorCode get_connectivity( MBEntityHandle handle,
                                        std::vector<MBEntityHandle>& connect,
                                        bool topological = false ) const;
  
  virtual MBErrorCode get_connectivity( MBEntityHandle handle,
                                        MBEntityHandle const*& connect,
                                        int &connect_length,
                                        bool topological = false,
                                        std::vector<MBEntityHandle>* storage = 0
                                       ) const;

  virtual MBErrorCode set_connectivity( MBEntityHandle handle,
                                        MBEntityHandle const* connect,
                                        int connect_length );
  
  virtual MBEntityHandle* get_connectivity_array();
 
   /***************** Methods from EntitySequence *****************/

    /* Replace the ElementSequence implementation of this method with
     * one that always returns zero, because we cannot re-use handles
     * that are within a ScdElementData
     */
  virtual int values_per_entity() const;

  virtual EntitySequence* split( MBEntityHandle here );

  virtual SequenceData* create_data_subset( MBEntityHandle start_handle,
                                            MBEntityHandle end_handle ) const;

  virtual void get_const_memory_use( unsigned long& bytes_per_entity,
                                     unsigned long& size_of_sequence ) const;

protected:
  SweptElementSeq( SweptElementSeq& split_from, MBEntityHandle here )
    : ElementSequence( split_from, here )
    {}
};

#endif
