AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) $(QVDUAL_INCLUDES)
bin_PROGRAMS = qvdual
qvdual_SOURCES = \
  CropTool.cpp \
  CropTool.hpp \
  CropToolpopup.ui \
  CropToolpopup.ui.h \
  DrawDual.cpp \
  DrawDual.hpp \
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
  moc_CropToolpopup.o

LDADD = $(extra_depends) $(QVDUAL_LIBS) $(top_builddir)/libMOAB.la 
AM_LDFLAGS = -R $(VTK_LIB_DIR)
qvdual_DEPENDENCIES = $(extra_depends) $(top_builddir)/libMOAB.la

uiQVDual.cpp: uiQVDual.h uiQVDual.ui 
	$(UIC) -impl uiQVDual.h -o $@ uiQVDual.ui 

CropToolpopup.cpp: CropToolpopup.h CropToolpopup.ui 
	$(UIC) -impl CropToolpopup.h -o $@ CropToolpopup.ui 

moc_uiQVDual.cpp: uiQVDual.h
	$(MOC) -o $@ uiQVDual.h 

moc_CropToolpopup.cpp: CropToolpopup.h 
	$(MOC) -o $@ CropToolpopup.h 

uiQVDual.ui : uiQVDual.ui.h

uiQVDual.o: uiQVDual.cpp

.ui.h: 
	$(UIC) -o $@ $<

main.o: uiQVDual.h

