MOAB_DIR=../..
include $(MOAB_DIR)/MB.common

dual_test: dual_test.cpp
	g++ $(CXXFLAGS) $(DEBUG_FLAG) -Wl,-rpath,$(MOAB_DIR) \
	-o $@ -I $(MOAB_DIR) $< $(PCOM_APP_LIBS)

