#ifndef READ_PARALLEL_HPP
#define READ_PARALLEL_HPP

#include "MBForward.hpp"
#include "MBReaderIface.hpp"

class MBReadUtilIface;

class ReadParallel : public MBReaderIface
{
   
public:

  static MBReaderIface* factory( MBInterface* );

    //! load a file
  MBErrorCode load_file(const char *file_name,
                        MBEntityHandle& file_set,
                        const FileOptions &opts,
                        const int* material_set_list,
                        const int num_material_sets );
  
    //! Constructor
  ReadParallel(MBInterface* impl = NULL) {mbImpl = impl;};

   //! Destructor
  virtual ~ReadParallel() {}

protected:

private:
  MBInterface *mbImpl;
  
  MBErrorCode delete_nonlocal_entities(std::string &partition_name,
                                       MBEntityHandle file_set);
};

#endif
