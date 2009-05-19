#include "MBInterface.hpp"
#include "MBReaderIface.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

class MBReadUtilIface;

class ReadMCNP5 : public MBReaderIface
{

public:
  // factory method
  static MBReaderIface* factory( MBInterface* );
  
  MBErrorCode load_file(const char*       fname,
                        MBEntityHandle    &input_meshset,
                        const FileOptions &options,
                        const char*       set_tag_name, /* not used */
                        const int*        material_set_list,
                        const int         num_material_sets );

  // Constructor
  ReadMCNP5(MBInterface* impl = NULL);

  // Destructor
  virtual ~ReadMCNP5();
  
protected:
  
private:
  // Constants
  static const double PI;
  static const double C2PI;
  static const double CPI;

  enum coordinate_system { NO_SYSTEM,
                           CARTESIAN,
                           CYLINDRICAL,
                           SPHERICAL };
  enum particle { NEUTRON,
                  PHOTON,
                  ELECTRON };
  
  // Read mesh interface
  MBReadUtilIface* readMeshIface;
  
  // MOAB Interface
  MBInterface* MBI;
  
  MBErrorCode create_tags( MBTag &date_and_time_tag,  
                           MBTag &title_tag,
                           MBTag &nps_tag, 
                           MBTag &tally_number_tag,   
                           MBTag &tally_comment_tag,
                           MBTag &tally_particle_tag, 
                           MBTag &tally_coord_sys_tag,
                           MBTag &tally_tag,          
                           MBTag &error_tag );

  MBErrorCode read_file_header( std::fstream      &file,
                                bool              debug,
                                char              date_and_time[100], 
                                char              title[100], 
                                unsigned long int &nps );

  MBErrorCode set_header_tags( MBEntityHandle             output_meshset, 
                                        char              date_and_time[100],
                                        char              title[100],
                                        unsigned long int nps,
                                        MBTag             data_and_time_tag,
                                        MBTag             title_tag,
                                        MBTag             nps_tag );

  MBErrorCode read_tally_header( std::fstream &file,
                                 bool         debug,
                                 unsigned int &tally_number,
                                 char         tally_comment[100],
                                 particle     &tally_particle );

  MBErrorCode get_tally_particle( std::string a,
                                  bool        debug,
                                  particle    &tally_particle );

  MBErrorCode read_mesh_planes( std::fstream         &file, 
                                bool                 debug, 
                                std::vector<double>  planes[3], 
                                coordinate_system    &coord_sys);

  MBErrorCode get_mesh_plane( std::istringstream &ss, std::vector<double> &plane);

  MBErrorCode read_element_values_and_errors( std::fstream        &file,
                                              bool                debug,
                                              std::vector<double> planes[3],
                                              int                 n_chopped_x2_planes,
                                              double              values[],
                                              double              errors[] );

  MBErrorCode set_tally_tags( MBEntityHandle    tally_meshset,
                              unsigned int      tally_number,
                              char              tally_comment[100],
                              particle          tally_particle,
                              coordinate_system tally_coord_sys,
                              MBTag             tally_number_tag, 
                              MBTag             tally_comment_tag,
                              MBTag             tally_particle_tag,
                              MBTag             tally_coord_sys_tag );

  MBErrorCode create_vertices( std::vector<double> planes[3],
                               //MBRange             &vert_handles,
                               MBEntityHandle      &start_vert,
                               coordinate_system   coord_sys,
                               MBEntityHandle      tally_meshset );
 
  MBErrorCode create_elements( bool                debug, 
                               std::vector<double> planes[3],
                               int                 n_chopped_x2_planes,
                               //MBRange             vert_handles,
                               MBEntityHandle      start_vert,
                               double              values[],
                               double              errors[],
                               MBTag               tally_tag,
                               MBTag               error_tag,
                               MBEntityHandle      tally_meshset );

  MBErrorCode average_with_existing_tally( unsigned long int &new_nps,
                                           unsigned long int nps,
                                           unsigned int      tally_number,
                                           MBTag             tally_number_tag,
                                           MBTag             nps_tag,
                                           MBTag             tally_tag,
                                           MBTag             error_tag,
                                           double            values[],
                                           double            errors[],
                                           unsigned int      n_elements );
  
  MBErrorCode transform_point_to_cartesian(double *in, 
                                           double *out, 
                                           coordinate_system coord_sys);

  MBErrorCode average_tally_values(const unsigned long int nps0, 
                                   const unsigned long int nps1,
                                   double                  *values0,
                                   const double            *values1,
                                   double                  *errors0,
                                   const double            *errors1,
                                   const unsigned long int n_values);
  
};
