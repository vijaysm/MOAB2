AM_CPPFLAGS += -I$(top_srcdir)/src -I$(top_builddir)/src $(QVDUAL_INCLUDES)
bin_PROGRAMS = qvdual
qvdual_SOURCES = \
  CropTool.cpp \
  CropTool.hpp \
  CropToolpopup.ui \
  CropToolpopup.ui.h \
  DrawDual.cpp \
  DrawDual.hpp \
  main.cpp \
  SheetDiagramPopup.ui.h \
  SheetDiagramPopup.ui \
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
  moc_CropToolpopup.o \
  SheetDiagramPopup.o \
  moc_SheetDiagramPopup.o


LDADD = $(extra_depends) $(QVDUAL_LIBS) $(top_builddir)/src/libMOAB.la 
AM_LDFLAGS = -R $(VTK_LIB_DIR)
qvdual_DEPENDENCIES = $(extra_depends) $(top_builddir)/src/libMOAB.la

uiQVDual.cpp: uiQVDual.h uiQVDual.ui 
	$(UIC) -impl uiQVDual.h -o $@ uiQVDual.ui 

CropToolpopup.cpp: CropToolpopup.h CropToolpopup.ui 
	$(UIC) -impl CropToolpopup.h -o $@ CropToolpopup.ui 

SheetDiagramPopup.cpp: SheetDiagramPopup.h SheetDiagramPopup.ui 
	$(UIC) -impl SheetDiagramPopup.h -o $@ SheetDiagramPopup.ui 

moc_uiQVDual.cpp: uiQVDual.h
	$(MOC) -o $@ uiQVDual.h 

moc_CropToolpopup.cpp: CropToolpopup.h 
	$(MOC) -o $@ CropToolpopup.h 

moc_SheetDiagramPopup.cpp: SheetDiagramPopup.h 
	$(MOC) -o $@ SheetDiagramPopup.h 

uiQVDual.ui : uiQVDual.ui.h

uiQVDual.o: uiQVDual.cpp

CropToolpopup.ui : CropToolpopup.ui.h

CropToolpopup.o: CropToolpopup.cpp

SheetDiagramPopup.ui : SheetDiagramPopup.ui.h

SheetDiagramPopup.o: SheetDiagramPopup.cpp

.ui.h: 
	$(UIC) -o $@ $<

main.o: uiQVDual.h

