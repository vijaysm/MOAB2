DEBUG_FLAG = -g
#DEBUG_FLAG = -pg -O2

CXXFLAGS =

dual_test: dual_test.cpp
	g++ $(CXXFLAGS) $(DEBUG_FLAG) -Wl,-rpath,../.. -o dual_test -DLINUX -I ../.. dual_test.cpp -L ../.. -lMOAB -L/cubit/netcdf/lib -lnetcdf_c++ -lnetcdf

