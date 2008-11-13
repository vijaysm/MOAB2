#ifndef READ_PARALLEL_HPP
#define READ_PARALLEL_HPP

#include "MBForward.hpp"
#include "MBReaderIface.hpp"

#include <string>

class MBReadUtilIface;
class MBParallelComm;

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
  
    //! load multiple files
  MBErrorCode load_file(const char **file_names,
                        const int num_files,
                        MBEntityHandle& file_set,
                        const FileOptions &opts,
                        const int* material_set_list,
                        const int num_material_sets );
  
    //! Constructor
  ReadParallel(MBInterface* impl = NULL, MBParallelComm *pc = NULL);

   //! Destructor
  virtual ~ReadParallel() {}

  static const char *parallelOptsNames[];
  
  enum ParallelOpts {POPT_NONE=0, POPT_BCAST, POPT_BCAST_DELETE, 
                     POPT_READ_DELETE, POPT_READ_PARALLEL,
                     POPT_FORMAT, POPT_DEFAULT};
  
protected:

private:
  MBErrorCode load_file(const char **file_names,
                        const int num_files,
                        MBEntityHandle& file_set,
                        int parallel_mode, 
                        std::string &partition_tag_name, 
                        std::vector<int> &partition_tag_vals, 
                        bool distrib,
                        std::vector<int> &pa_vec,
                        const int* material_set_list,
                        const int num_material_sets,
                        const FileOptions &opts,
                        const int reader_rank,
                        const bool cputime,
                        const int resolve_dim,
                        const int shared_dim,
                        const int ghost_dim,
                        const int bridge_dim,
                        const int num_layers);

  MBErrorCode delete_nonlocal_entities(std::string &ptag_name,
                                       std::vector<int> &ptag_vals,
                                       bool distribute,
                                       MBEntityHandle file_set);
  
  MBErrorCode delete_nonlocal_entities(MBEntityHandle file_set);

  MBInterface *mbImpl;

    // each reader can keep track of its own pcomm
  MBParallelComm *myPcomm;
};

inline MBErrorCode ReadParallel::load_file(const char *file_name,
                                           MBEntityHandle& file_set,
                                           const FileOptions &opts,
                                           const int* material_set_list,
                                           const int num_material_sets ) 
{
  return load_file(&file_name, 1, file_set, opts, 
                   material_set_list, num_material_sets);
}
  
#endif
