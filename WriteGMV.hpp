//-------------------------------------------------------------------------
// Filename      : WriteGMV.hpp
//
// Purpose       : Writer template
//
// Special Notes : 
//
// Creator       : Tim Tautges
//
// Date          : 2/04
//
//-------------------------------------------------------------------------

#ifndef WRITEGMV_HPP
#define WRITEGMV_HPP

#include "MBInterface.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBWriterIface.hpp"

//! Output Exodus File for VERDE
class MB_DLL_EXPORT WriteGMV : public MBWriterIface
{
 
public:

   //! Constructor
   WriteGMV(MBInterface *impl);

   //! Destructor
  virtual ~WriteGMV();

  static MBWriterIface* factory( MBInterface* );

  MBErrorCode write_file( const char* filename,
                          const bool overwite,
                          const MBEntityHandle* output_sets,
                          const int num_output_sets,
                          std::vector<std::string>& qa_list,
                          int requested_dimension );

    //! writes out a mesh file
  MBErrorCode write_file(const char *file_name,
                         const MBEntityHandle output_set,
                         const int user_dimension = 3,
                         const bool mesh = true,
                         const bool poly_mesh = true);
  
protected:

private:

    //! interface instance
  MBInterface *mbImpl;
  MBWriteUtilIface* mWriteIface;
  
    //! Meshset Handle for the mesh that is currently being written
  MBEntityHandle mCurrentMeshHandle;

  //! Cached tags for reading.  Note that all these tags are defined when the
  //! core is initialized.
  MBTag mMaterialSetTag;
  MBTag mDirichletSetTag;
  MBTag mNeumannSetTag;
  MBTag mHasMidNodesTag;
  MBTag mGeomDimensionTag;
  MBTag mGlobalIdTag;

  static const char *gmvTypeNames[MBMAXTYPE];
  
  MBErrorCode local_write_mesh(const char *file_name,
                               const MBEntityHandle output_set,
                               const int user_dimension,
                               const bool mesh,
                               const bool poly_mesh);
};

#endif
