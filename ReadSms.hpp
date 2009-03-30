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
                        const char* set_tag_name,
                        const int* set_tag_values,
                        int num_set_tag_values );
  
    //! Constructor
  ReadSms(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadSms();

private:

  MBErrorCode load_file_impl( const char *file_name );
  
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
