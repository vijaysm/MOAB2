MOAB_DIR=../..
include $(MOAB_DIR)/MB.common

h5file_test: h5file_test.cpp
	g++ $(CXXFLAGS) $(DEBUG_FLAG) -Wl,-rpath,$(MOAB_DIR) \
	-o $@ -I $(MOAB_DIR) $< $(PCOM_APP_LIBS)

