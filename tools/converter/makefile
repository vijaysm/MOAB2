include ../../MB.$(MACHINE_TYPE)
include ../../MB.local


convert: convert.cpp ../../libMOAB.so
	$(CXX) $(DEBUG_FLAG) $(MACH_CXXFLAGS) \
        convert.cpp -o convert -I../.. -L../.. -lMOAB

../../libMOAB.so:
	cd ../.. && make libMOAB.so

clean:
	rm -f convert
