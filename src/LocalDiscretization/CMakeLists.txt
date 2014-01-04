set( MOAB_LOC_DSCR_SRCS ElemEvaluator.cpp
                        LinearHex.cpp
                        LinearQuad.cpp
                        LinearTet.cpp
                        QuadraticHex.cpp )

set( MOAB_LOCALDSRC_Install_Headers moab/ElemEvaluator.hpp
                                    moab/LinearHex.hpp
                                    moab/LinearQuad.hpp
	                                  moab/LinearTet.hpp
                                    moab/QuadraticHex.hpp )

add_library( moab_loc_discr STATIC
                ${MOAB_LOC_DSCR_SRCS}
                ${MOAB_LOCALDSRC_Install_Headers} )

target_link_libraries( moab_loc_discr )

set_target_properties( moab_loc_discr 
                       PROPERTIES COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}")

include(AutoconfHeader)
moab_install_headers(${MOAB_LOCALDSRC_Install_Headers})
