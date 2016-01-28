#
# The following few lines get the library install location for MOAB
# and puts it in MOAB_LIB_DIR; if you just want
# a simple application makefile, just comment out the next 3 lines and
# substitute your own value for $(MOAB_LIB_DIR), or set an environment
# variable pointing to the right place.
MOAB_LIB_DIR = @libdir@

include $(MOAB_LIB_DIR)/moab.make

CPPFLAGS = ${MOAB_INCLUDES}
CXXFLAGS = -g
CXX      = $(MOAB_CXX)

default: all

all: moabuse1 moabuse2 moabuse3 moabuse4

moabuse1: moabuse1.o
	$(CXX) -o $@ $< $(MOAB_LIBS_LINK)

moabuse2: moabuse2.o
	$(CXX) -o $@ $< $(MOAB_LIBS_LINK)

moabuse3: moabuse3.o
	$(CXX) -o $@ $< $(MOAB_LIBS_LINK)

moabuse4: moabuse4.o
	$(CXX) -o $@ $< $(MOAB_LIBS_LINK)

%.o : %.cpp
	$(CXX) -c $(MOAB_INCLUDES) $<

clean:
	rm -rf *.o *.vtk *.h5m moabuse1 moabuse2 moabuse3 moabuse4

