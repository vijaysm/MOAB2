#ifndef READ_PARALLEL_HPP
#define READ_PARALLEL_HPP

#include "moab/Forward.hpp"
#include "moab/ReaderIface.hpp"
#include "DebugOutput.hpp"

#include <string>

namespace moab {

class ReadUtilIface;
class ParallelComm;

class ReadParallel
{
   
public:

  static ReaderIface* factory( Interface* );

    //! load a file
  ErrorCode load_file(const char *file_name,
                        const EntityHandle* file_set,
                        const FileOptions &opts,
                        const ReaderIface::IDTag* subset_list = 0,
                        int subset_list_length = 0,
                        const Tag* file_id_tag = 0 );
  
    //! load multiple files
  ErrorCode load_file(const char **file_names,
                        const int num_files,
                        const EntityHandle* file_set,
                        const FileOptions &opts,
                        const ReaderIface::IDTag* subset_list = 0,
                        int subset_list_length = 0,
                        const Tag* file_id_tag = 0 );
  
  ErrorCode load_file(const char **file_names,
                        const int num_files,
                        const EntityHandle* file_set,
                        int parallel_mode, 
                        std::string &partition_tag_name, 
                        std::vector<int> &partition_tag_vals, 
                        bool distrib,
                        bool partition_by_rank,
                        std::vector<int> &pa_vec,
                        const FileOptions &opts,
                        const ReaderIface::IDTag* subset_list,
                        int subset_list_length,
                        const Tag* file_id_tag,
                        const int reader_rank,
                        const bool cputime,
                        const int resolve_dim,
                        const int shared_dim,
                        const int ghost_dim,
                        const int bridge_dim,
                        const int num_layers);
    //! Constructor
  ReadParallel(Interface* impl = NULL, ParallelComm *pc = NULL);

   //! Destructor
  virtual ~ReadParallel() {}

  static const char *parallelOptsNames[];
  
  enum ParallelActions {PA_READ=0, 
                        PA_READ_PART, 
                        PA_BROADCAST, 
                        PA_DELETE_NONLOCAL,
                        PA_CHECK_GIDS_SERIAL, 
                        PA_GET_FILESET_ENTS, 
                        PA_RESOLVE_SHARED_ENTS,
                        PA_EXCHANGE_GHOSTS, 
                        PA_PRINT_PARALLEL};

  static const char *ParallelActionsNames[];
  
  enum ParallelOpts { POPT_NONE=0, 
                      POPT_BCAST, 
                      POPT_BCAST_DELETE, 
                      POPT_READ_DELETE, 
                      POPT_READ_PART, 
                      POPT_DEFAULT};

    //! PUBLIC TO ALLOW TESTING
  ErrorCode delete_nonlocal_entities(std::string &ptag_name,
                                       std::vector<int> &ptag_vals,
                                       bool distribute,
                                       EntityHandle file_set);
  
  ErrorCode delete_nonlocal_entities(EntityHandle file_set);

protected:
  ErrorCode create_partition_sets( std::string &ptag_name,
                                     EntityHandle file_set );

private:

  Interface *mbImpl;

    // each reader can keep track of its own pcomm
  ParallelComm *myPcomm;
  
  DebugOutput myDebug;
};

inline ErrorCode ReadParallel::load_file(const char *file_name,
                                           const EntityHandle* file_set,
                                           const FileOptions &opts,
                                           const ReaderIface::IDTag* subset_list,
                                           int subset_list_length,
                                           const Tag* file_id_tag )
{
  return load_file(&file_name, 1, file_set, opts, 
                   subset_list, subset_list_length, file_id_tag);
}

} // namespace moab

#endif
