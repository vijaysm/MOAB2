DEBUG_FLAG = -O2
#DEBUG_FLAG = -pg -O2

CXXFLAGS = -DREALMAP

CUBIT = ${HOME}/cubit

perf: perf.cpp
	g++ $(CXXFLAGS) $(DEBUG_FLAG) -o perf -DLINUX -DIS_BUILDING_MB -I ../.. perf.cpp -L ../.. -lMOAB -L/cubit/netcdf/lib -lnetcdf_c++ -lnetcdf

cubit_perf: cubit_perf.o
	g++ ${DEBUG_FLAG} -I${CUBIT}/mdb -I${CUBIT}/util cubit_perf.o -L$(CUBIT) -lcubiti -o cubit_perf -L$(CUBIT)/showviz/lib -lshowviz_geom_mesh -lshowviz_mesh -lshowviz_geom -lshowviz \
-L/cubit/VTK-4.2.2/lib/vtk -lvtkHybrid -lvtkIO -lvtkRendering -lvtkftgl -lvtkfreetype -lvtkImaging -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtkexpat -lvtkjpeg -lvtkpng -lvtktiff -lvtkzlib -L$(CUBIT)/df -lcbt_df  -L$(CUBIT)/mdb -lcubit_mdb -L$(CUBIT)/geom/virtual -lcubit_virtual -L$(CUBIT)/geom/facet -lcubit_facet -L$(CUBIT)/geom/facetbool -lcubit_facetboolstub -L$(CUBIT)/geom -lcubit_geom -L$(CUBIT)/geom/Cholla -lcholla -L$(CUBIT)/util -lcubit_util -L$(CUBIT)/camal/lib/Linux -lCMLTri102 \
        -L$(CUBIT)/geom/ACIS -lcubit_ACIS -L/cubit/acis/acis11.5/bin/linux_so -lxstep -lxacis2k -lxcore2k -lxutil -lacisstep -lxiges -lacisiges -lSpaAWarp -lSpaAVis -lSpaASurf -lSpaALops -lSpaABlend -lSpaACIS -lSpaBase -L/cubit/lib -llpk -L/cubit/lib -lexoIIv2c326 -L/cubit/netcdf/lib -lnetcdf  -L/cubit/verdict1.0.1/lib -lverdict101-3.2  -L/usr/X11R6/lib -lXt -lXext -lX11 

cubit_perf.o: cubit_perf.cpp
	g++ -c -DLINUX ${DEBUG_FLAG} -I${CUBIT}/mdb -I${CUBIT}/util cubit_perf.cpp
