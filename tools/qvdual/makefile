QT_DIR = /usr/local/qt-x11-free-3.3.2

VTK_LIB_DIR = /usr/local/ParaViewComplete/bin
VTK_DIR = /usr/local/ParaViewComplete/VTK

DEBUG_FLAG = -g

GRAPHVIZ_LIB_DIR = /usr/local/lib/graphviz
GRAPHVIZ_INCLUDE_DIR = /usr/local/include/graphviz

NETCDF_DIR = /cubit/netcdf

MOAB_DIR = ../..

X11_LIB_DIR = /usr/X11R6/lib 

INCLUDE_FLAGS = -I../.. -I. -I${GRAPHVIZ_INCLUDE_DIR} -I${VTK_DIR} -I${VTK_DIR}/Utilities -I${VTK_DIR}/Parallel -I${VTK_DIR}/Hybrid -I${VTK_DIR}/Rendering -I${VTK_DIR}/IO -I${VTK_DIR}/Imaging -I${VTK_DIR}/Graphics -I${VTK_DIR}/GenericFiltering -I${VTK_DIR}/Filtering -I${VTK_DIR}/Common -I${VTK_DIR}/Common/Testing/Cxx -I${VTK_DIR}/Utilities/DICOMParser -I${VTK_DIR}/Utilities/freetype/include -I/usr/share/doc/NVIDIA_GLX-1.0/include -I/usr/X11R6/include -I${QT_DIR}/include 

CXXFLAGS = -Wno-deprecated ${DEBUG_FLAG} $(INCLUDE_FLAGS)

QVTK_SRC_OBJS = \
	QVTKWidget.o \
	moc_QVTKWidget.o

QVTKWidgetPlugin_SRC_OBJS = \
	QVTKWidget.o \
	QVTKWidgetPlugin.o \
	moc_QVTKWidget.o

QVDual_SRC_OBJS = \
	main.o \
	CropTool.o \
	DrawDual.o \
	uiQVDual.o \
	moc_uiQVDual.o \
	CropToolpopup.o \
	moc_CropToolpopup.o \
	vtkMOABUtils.o \
	vtkMOABReader.o

GENERATED_QVDUAL_QT_FILES =  \
	uiQVDual.h \
	uiQVDual.cxx \
	moc_uiQVDual.cxx \
	CropToolpopup.h \
	CropToolpopup.cxx \
	moc_CropToolpopup.cxx

GENERATED_QVTK_QT_FILES = moc_QVTKWidget.cxx 

default: QVDual

libQVTK.a: $(QVTK_SRC_OBJS) $(QVTK_EXTERNAL_OBJS)  $(QVTK_DEPEND_LIBS)
	ar cr libQVTK.a  $(QVTK_SRC_OBJS)
	ranlib libQVTK.a 

libQVTKWidgetPlugin.so: $(QVTKWidgetPlugin_SRC_OBJS)
	${CXX} -fPIC -shared -Wl,-soname,libQVTKWidgetPlugin.so -o libQVTKWidgetPlugin.so $(QVTKWidgetPlugin_SRC_OBJS) -L${VTK_LIB_DIR} -L${QT_DIR}/lib -L/usr/X11R6/lib -lvtkRendering -lvtkImaging -lvtkGraphics -lvtkFiltering -lvtkCommon -lqt -lGLU -lGL -lSM -lICE -lSM -lICE -lSM -lICE -lSM -lICE -lSM -lICE -lX11 -lXext -lX11 -lXext -lX11 -lXext -lX11 -lXext -lX11 -lXext -lvtkIO -lvtkFiltering -lvtkDICOMParser -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lvtkexpat -lvtkftgl -lGL -lvtkfreetype -lSM -lICE -lX11 -lXext -lm -lXt -lpthread -ldl -lm -lgcc

QVDual: $(QVDual_SRC_OBJS) libQVTK.a
	${CXX} ${DEBUG_FLAG} -Wno-deprecated $(QVDual_SRC_OBJS) -o QVDual -rdynamic -L. -L${VTK_LIB_DIR} -L${MOAB_DIR} -L${NETCDF_DIR}/lib -L${GRAPHVIZ_LIB_DIR} -L${QT_DIR}/lib -L${X11_LIB_DIR} -lMOAB -lnetcdf_c++ -lnetcdf -ldotneato -lgraph -lQVTK -lqt -lvtkRendering -lvtkGraphics -lvtkIO -lvtkCommon -lvtkFiltering -lGLU -lGL -lSM -lICE -lX11 -lXext -lvtkParallel -lvtkHybrid -lvtkRendering -lvtkIO -lvtkImaging -lvtkFiltering -lvtkftgl -lGL -lvtkfreetype -lvtkDICOMParser -lvtkpng -lvtktiff -lvtkjpeg -lvtkexpat -lSM -lICE -lX11 -lXext -lvtkCommon -lm -lXt -lpthread -ldl -lm

clean:
	rm -rf QVDual libQVTK.a libQVTKWidgetPlugin.so ${GENERATED_QVDUAL_QT_FILES} ${GENERATED_QVTK_QT_FILES} 

.SUFFIXES: .o .cpp .cxx
uiQVDual.h: uiQVDual.ui
	${QT_DIR}/bin/uic -o uiQVDual.h uiQVDual.ui 

uiQVDual.o: uiQVDual.cxx uiQVDual.h CropToolpopup.h
	${CXX} ${CXXFLAGS} -c -o uiQVDual.o uiQVDual.cxx

uiQVDual.cxx: uiQVDual.h
	${QT_DIR}/bin/uic -impl uiQVDual.h -o uiQVDual.cxx uiQVDual.ui 

moc_uiQVDual.o: moc_uiQVDual.cxx
	${CXX} ${CXXFLAGS} -c -o moc_uiQVDual.o moc_uiQVDual.cxx

moc_uiQVDual.cxx: uiQVDual.h
	${QT_DIR}/bin/moc -o moc_uiQVDual.cxx uiQVDual.h 

CropToolpopup.h: CropToolpopup.ui
	${QT_DIR}/bin/uic -o CropToolpopup.h CropToolpopup.ui 

CropToolpopup.o: CropToolpopup.cxx
	${CXX} ${CXXFLAGS} -c CropToolpopup.cxx

CropToolpopup.cxx: CropToolpopup.h
	${QT_DIR}/bin/uic -impl CropToolpopup.h -o CropToolpopup.cxx CropToolpopup.ui 

moc_CropToolpopup.cxx: CropToolpopup.h
	${QT_DIR}/bin/moc -o moc_CropToolpopup.cxx CropToolpopup.h 

moc_CropToolpopup.o: moc_CropToolpopup.cxx
	${CXX} ${CXXFLAGS} -c moc_CropToolpopup.cxx

moc_QVTKWidget.cxx: QVTKWidget.h
	${QT_DIR}/bin/moc -o moc_QVTKWidget.cxx QVTKWidget.h 

QVTKWidget.o: QVTKWidget.cxx
	${CXX} ${CXXFLAGS} -c -o QVTKWidget.o    QVTKWidget.cxx

moc_QVTKWidget.o: moc_QVTKWidget.cxx
	${CXX} ${CXXFLAGS} -o moc_QVTKWidget.o    -Wno-deprecated -g $(INCLUDE_FLAGS)  -c moc_QVTKWidget.cxx

QVTKWidgetPlugin.o: QVTKWidgetPlugin.cxx
	${CXX} -o QVTKWidgetPlugin.o -DQVTKWidgetPlugin_EXPORTS   -Wno-deprecated -g -fPIC $(INCLUDE_FLAGS)  -c QVTKWidgetPlugin.cxx

.cxx.o:
	${CXX} -c ${CXXFLAGS} $<

.cpp.o:
	${CXX} -c ${CXXFLAGS} $<
