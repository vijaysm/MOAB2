#ifndef READ_SMS_HPP
#define READ_SMS_HPP

#include "MBForward.hpp"
#include "MBReaderIface.hpp"
#include "MBRange.hpp"
#include <vector>

class MBReadUtilIface;

// Base class for binary and ASCII readers
class ReadSms : public MBReaderIface
{
   
public:

    //! factory method 
  static MBReaderIface* factory( MBInterface* );

  MBErrorCode load_file(const char *file_name,
                        MBEntityHandle& file_set,
                        const FileOptions& opts,
                        const int* material_set_list,
                        int num_material_sets );
  
    //! Constructor
  ReadSms(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadSms();

private:

  MBErrorCode load_file_impl( const char *file_name,
                              const int* material_set_list,
                              const int num_material_sets);
  
  MBErrorCode get_set(std::vector<MBEntityHandle> *sets,
                      int set_type, int set_id,
                      MBTag set_tag,
                      MBEntityHandle &this_set);

  MBErrorCode read_parallel_info(FILE *file_ptr);

  MBReadUtilIface* readMeshIface;

    //! interface instance
  MBInterface* mdbImpl;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;
  
  MBTag globalId, paramCoords, geomDimension;
};


#endif
