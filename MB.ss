#
# location of exodus libs and makedepend
#
CUBITROOT = $(CUBIT)

NETCDF_INCLUDE = -I$(CUBITROOT)/include
NETCDF_LIB_DIR = $(CUBITROOT)/lib
NETCDF_LIBS = -lnetcdf_c++ -lnetcdf

# Machine-dependent include paths
CC_INCLUDE = -I/net/sasn232/opt/SUNWspro.old/SC5.0/include/CCios \
        -I/opt/SUNWspro/WS6U1/include/CCios \
        -I/opt/G111G1SUNWspro/WS6U1/include/CC \
        -I/opt/SUNWspro/WS6U1/include/CC/Cstd

MACH_CXXFLAGS = -DSOLARIS
#MACH_LINK_LIBS = -liostream -lC

# Some controls on echoing of compile commands; for verbose output, comment
# out the following two definitions
PREFIX = @
ECHO_COMMAND = @ echo "Compiling(${DEBUG_FLAG}): $<"


# Some compiler definitions - check your environment to make sure the correct
# compilers are being referenced
#
# CC:
CC = cc
CXX = CC
LINKER = ${CXX}
MAKE = make 
ARCHIVER = ${CXX} -xar -o
PIC = -Kpic
SO_LINKER = $(LINKER) $(PIC) -G -lCstd -h $(notdir $@)
RUNTIME_LINKER_LIB = -ldl

# MAKEDEPEND:
MAKEDEPEND = $(CUBITROOT)/bin/makedepend -f -

# location of templates for template classes, helps make clean_all
TEMPLATE_DIR = SunWS_cache

