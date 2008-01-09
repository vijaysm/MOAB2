#ifndef READ_PARALLEL_HPP
#define READ_PARALLEL_HPP

#include "MBForward.hpp"
#include "MBReaderIface.hpp"

#include <string>

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
  
  MBErrorCode load_file(const char *file_name,
                        MBEntityHandle& file_set,
                        int parallel_mode, 
                        std::string &partition_tag_name, 
                        std::vector<int> &partition_tag_vals, 
                        bool distrib,
                        std::vector<int> &pa_vec,
                        const int* material_set_list,
                        const int num_material_sets,
                        const FileOptions &opts,
                        int reader_rank,
                        bool cputime);

  MBErrorCode delete_nonlocal_entities(std::string &ptag_name,
                                       std::vector<int> &ptag_vals,
                                       bool distribute,
                                       MBEntityHandle file_set);
  
  MBErrorCode delete_nonlocal_entities(MBRange &partition_sets,
                                       MBEntityHandle file_set);
};

#endif
