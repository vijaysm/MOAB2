//-------------------------------------------------------------------------
// Filename      : ReadVtk.hpp
//
// Purpose       : Vtk reader
//
// Special Notes : Lots of code taken from verde implementation
//
//-------------------------------------------------------------------------

#ifndef READVTK_HPP
#define READVTK_HPP

#ifndef IS_BUILDING_MB
#error "ReadVtk.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"
#include "MBReaderIface.hpp"

class MBReadUtilIface;

class ReadVtk : public MBReaderIface
{
   
public:

  static MBReaderIface* factory( MBInterface* );

    //! load a file
  MBErrorCode load_file(const char *file_name,
                        const int* material_set_list,
                        const int num_material_sets );
  
    //! Constructor
  ReadVtk(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadVtk() {}

private:

  MBReadUtilIface* readMeshIface;

  //------------member variables ------------//

    //! interface instance
  MBInterface* mdbImpl;
  
    //! file name
  std::string fileName;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;
};

#endif




