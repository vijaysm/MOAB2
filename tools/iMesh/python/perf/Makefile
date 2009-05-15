all: perf.c
	gcc -o perf perf.c -I/home/porterj/moab/include \
	-L/home/porterj/moab/lib -L/usr/local/hdf5serial/lib \
	-liMesh -lm -lgcc_s -lc -lMOAB -lnetcdf_c++ -lnetcdf -lhdf5 -lstdc++ \
	-O2 -DNDEBUG -pthread -g -Wall -Wstrict-prototypes -fPIC \
	-fno-strict-aliasing -fwrapv