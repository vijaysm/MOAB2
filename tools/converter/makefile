MOAB_DIR = ../..
include $(MOAB_DIR)/MB.common

convert: convert.cpp 
	$(CXX) $(DEBUG_FLAG) $(MACH_CXXFLAGS) \
	-I$(MOAB_DIR) -L$(MOAB_DIR) -lMOAB \
	convert.cpp -o convert $(MOAB_APP_LIBS)

../../libMOAB.so:
	cd ../.. && make libMOAB.so

clean:
	rm -f convert
