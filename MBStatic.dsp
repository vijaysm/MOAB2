# Microsoft Developer Studio Project File - Name="MDBStatic" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=MDBStatic - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "MDBStatic.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "MDBStatic.mak" CFG="MDBStatic - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "MDBStatic - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "MDBStatic - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "MDBStatic - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "MDBStatic___Win32_Release"
# PROP BASE Intermediate_Dir "MDBStatic___Win32_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "MDBStatic___Win32_Release"
# PROP Intermediate_Dir "MDBStatic___Win32_Release"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "$(CUBITROOT)\exodus\exodus3.07\include" /I "$(CUBITROOT)\netcdf\netcdf-3.4\include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /D "IS_BUILDING_MDB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "MDBStatic - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "MDBStatic___Win32_Debug"
# PROP BASE Intermediate_Dir "MDBStatic___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "MDBStatic___Win32_Debug"
# PROP Intermediate_Dir "MDBStatic___Win32_Debug"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "$(CUBITROOT)\exodus\exodus3.07\include" /I "$(CUBITROOT)\netcdf\netcdf-3.4\include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "IS_BUILDING_MDB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "MDBStatic - Win32 Release"
# Name "MDBStatic - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\AEntityFactory.cpp
# End Source File
# Begin Source File

SOURCE=.\DenseTagCollections.cpp
# End Source File
# Begin Source File

SOURCE=.\EntitySequence.cpp
# End Source File
# Begin Source File

SOURCE=.\EntitySequenceManager.cpp
# End Source File
# Begin Source File

SOURCE=.\ExoIIUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\HigherOrderFactory.cpp
# End Source File
# Begin Source File

SOURCE=.\HomXform.cpp
# End Source File
# Begin Source File

SOURCE=.\MDBBits.cpp
# End Source File
# Begin Source File

SOURCE=.\MDBCore.cpp
# End Source File
# Begin Source File

SOURCE=.\MDBDefines.cpp
# End Source File
# Begin Source File

SOURCE=.\MDBFactory.cpp
# End Source File
# Begin Source File

SOURCE=.\MDBMeshSet.cpp
# End Source File
# Begin Source File

SOURCE=.\MDBRange.cpp
# End Source File
# Begin Source File

SOURCE=.\MDBReadUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\MDBUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\MDBWriteUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\MOABCN.cpp
# End Source File
# Begin Source File

SOURCE=.\ReadExoII.cpp
# End Source File
# Begin Source File

SOURCE=.\ScdElementSeq.cpp
# End Source File
# Begin Source File

SOURCE=.\ScdVertexSeq.cpp
# End Source File
# Begin Source File

SOURCE=.\SparseTagCollections.cpp
# End Source File
# Begin Source File

SOURCE=.\TagServer.cpp
# End Source File
# Begin Source File

SOURCE=.\WriteExoII.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\AEntityFactory.hpp
# End Source File
# Begin Source File

SOURCE=.\DenseTagCollections.hpp
# End Source File
# Begin Source File

SOURCE=.\EntitySequence.hpp
# End Source File
# Begin Source File

SOURCE=.\EntitySequenceManager.hpp
# End Source File
# Begin Source File

SOURCE=.\ExoIIDefines.hpp
# End Source File
# Begin Source File

SOURCE=.\ExoIIInterface.hpp
# End Source File
# Begin Source File

SOURCE=.\ExoIIUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\HigherOrderFactory.hpp
# End Source File
# Begin Source File

SOURCE=.\HomXform.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBBits.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBCore.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBError.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBInterface.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBInternals.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBMeshSet.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBRange.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBReadUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBReadUtilIface.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBUnknownInterface.h
# End Source File
# Begin Source File

SOURCE=.\MDBWriteUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBWriteUtilIface.hpp
# End Source File
# Begin Source File

SOURCE=.\ReadExoII.hpp
# End Source File
# Begin Source File

SOURCE=.\ReadWriteDefines.h
# End Source File
# Begin Source File

SOURCE=.\ScdElementSeq.hpp
# End Source File
# Begin Source File

SOURCE=.\ScdVertexSeq.hpp
# End Source File
# Begin Source File

SOURCE=.\SparseTagCollections.hpp
# End Source File
# Begin Source File

SOURCE=.\TagServer.hpp
# End Source File
# Begin Source File

SOURCE=.\WriteExoII.hpp
# End Source File
# End Group
# End Target
# End Project
