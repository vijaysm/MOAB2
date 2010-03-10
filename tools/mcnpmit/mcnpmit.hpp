#include "MBCore.hpp"
#include "MBRange.hpp"
#include <iostream>
#define MCNP mc_instance()
#define BOXMIN_TAG "BOXMIN_TAG"
#define BOXMAX_TAG "BOXMAX_TAG"
#define TALLY_TAG  "TALLY_TAG"
#define ERROR_TAG  "ERROR_TAG"

#define MBI mb_instance()

enum MCNPError { MCNP_SUCCESS, MCNP_FAILURE, DONE };
enum { NOSYS, CARTESIAN, CYLINDRICAL, SPHERICAL };

class McnpData {
      
      public:
            // Constructor and Destructor
            McnpData ();
            ~McnpData ();

            // Coordinate system and rotation matrix
            int coord_system;
            double rotation_matrix[16];

            // Vertices and elements
            std::vector<MBEntityHandle> MCNP_vertices;
            std::vector<MBEntityHandle> MCNP_elems;
            MBRange                     vert_handles;
            MBRange                     elem_handles;

            // Tally data
            MBTag box_min_tag, box_max_tag;
            MBTag tally_tag;
            MBTag relerr_tag;

            // MCNP Meshtal file name
            std::string MCNP_filename;

            // Setting and retrieving coordinate sysem
            MCNPError set_coord_system(int);
            int get_coord_system();

            // Setting and retrieving roation matrix
            MCNPError set_rotation_matrix(double[16]); 
            double* get_rotation_matrix();

            // Set the filename
            MCNPError set_filename(std::string);
            std::string get_filename();

            // MCNP reading routines
            MCNPError read_mcnpfile(bool);
            MCNPError read_coord_system(std::string);
            MCNPError read_rotation_matrix(std::string, int);
            MCNPError make_elements(std::vector<double> [3], int*);
            MCNPError make_adjacencies(int*);
            MCNPError initialize_tags();
            MCNPError extract_tally_data(std::string, MBEntityHandle);

            // Transformation routine
            MCNPError transform_point(double*, double*, int, double*);

            // Parameters
            static const double pi   = 3.141592653589793;
            static const double c2pi = 0.1591549430918954;
            static const double cpi  = 0.3183098861837907;

};
