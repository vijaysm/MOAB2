AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_builddir)
bin_PROGRAMS = qvdual
qvdual_SOURCES = \
  CropTool.cpp \
  CropTool.hpp \
  CropToolpopup.ui \
  CropToolpopup.ui.h \
  DrawDual.cpp \
  DrawDual.hpp \
  QVTKWidget.cxx \
  QVTKWidget.h \
  QVTKWidgetPlugin.cxx \
  QVTKWidgetPlugin.h \
  main.cpp \
  uiQVDual.ui \
  uiQVDual.ui.h \
  vtkMOABReader.cxx \
  vtkMOABReader.h \
  vtkMOABUtils.cxx \
  vtkMOABUtils.h

extra_depends = \
  uiQVDual.o \
  moc_uiQVDual.o \
  CropToolpopup.o \
  CropToolpopup.moc.o \
  QVTKWidget.moc.o 

qvdual_LDADD = $(extra_depends) $(LIBS)

.ui.h: 
	$(UIC) -o $@ $<

.cxx.moc.cxx: 
	$(MOC) -o $@ $<

