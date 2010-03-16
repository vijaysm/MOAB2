include ${MOAB_LIB_DIR}/moab.make

SurfArea : SurfArea.o
	${CXX} $< ${MOAB_LIBS_LINK} -o $@

.cpp.o : 
	${CXX} ${MOAB_INCLUDES} -c $<

