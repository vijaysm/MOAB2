add_executable( test_oldinc test_oldinc.cpp)
set_target_properties( test_oldinc PROPERTIES COMPILE_FLAGS "${MOAB_DEFINES}" )
target_link_libraries( test_oldinc MOAB )
