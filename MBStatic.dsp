# Microsoft Developer Studio Project File - Name="MBStatic" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=MBStatic - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "MBStatic.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "MBStatic.mak" CFG="MBStatic - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "MBStatic - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "MBStatic - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "MBStatic - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Static_Release"
# PROP BASE Intermediate_Dir "Static_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Static_Release"
# PROP Intermediate_Dir "Static_Release"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "$(CUBITROOT)\exodus\exodus3.07\include" /I "$(CUBITROOT)\netcdf\netcdf-3.4.snl\include" /D "NDEBUG" /D "IS_BUILDING_MDB" /D "WIN32" /D "_MBCS" /D "_LIB" /D "IS_BUILDING_MB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "MBStatic - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "MBStatic___Win32_Debug"
# PROP BASE Intermediate_Dir "MBStatic___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Static_Debug"
# PROP Intermediate_Dir "Static_Debug"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "$(CUBITROOT)\exodus\exodus3.07\include" /I "$(CUBITROOT)\netcdf\netcdf-3.4.snl\include" /D "_DEBUG" /D "IS_BUILDING_MB" /D "WIN32" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
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

# Name "MBStatic - Win32 Release"
# Name "MBStatic - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\AEntityFactory.cpp
# End Source File
# Begin Source File

SOURCE=.\DenseTagCollections.cpp
# End Source File
# Begin Source File

SOURCE=.\DualTool.cpp
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

SOURCE=.\GeomTopoTool.cpp
# End Source File
# Begin Source File

SOURCE=.\HigherOrderFactory.cpp
# End Source File
# Begin Source File

SOURCE=.\HomXform.cpp
# End Source File
# Begin Source File

SOURCE=.\MBBits.cpp
# End Source File
# Begin Source File

SOURCE=.\MBCN.cpp
# End Source File
# Begin Source File

SOURCE=.\MBCore.cpp
# End Source File
# Begin Source File

SOURCE=.\MBFactory.cpp
# End Source File
# Begin Source File

SOURCE=.\MBMem.cpp
# End Source File
# Begin Source File

SOURCE=.\MBMeshSet.cpp
# End Source File
# Begin Source File

SOURCE=.\MBRange.cpp
# End Source File
# Begin Source File

SOURCE=.\MBReadUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\MBSkinner.cpp
# End Source File
# Begin Source File

SOURCE=.\MBUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\MBWriteUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\MeshTopoUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\PolyEntitySequence.cpp
# End Source File
# Begin Source File

SOURCE=.\ReadNCDF.cpp
# End Source File
# Begin Source File

SOURCE=.\ReadVtk.cpp
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

SOURCE=.\Tqdcfr.cpp
# End Source File
# Begin Source File

SOURCE=.\WriteGMV.cpp
# End Source File
# Begin Source File

SOURCE=.\WriteNCDF.cpp
# End Source File
# Begin Source File

SOURCE=.\WriteSLAC.cpp
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

SOURCE=.\DualTool.hpp
# End Source File
# Begin Source File

SOURCE=.\EntitySequence.hpp
# End Source File
# Begin Source File

SOURCE=.\EntitySequenceManager.hpp
# End Source File
# Begin Source File

SOURCE=.\ExoIIUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\GeomTopoTool.hpp
# End Source File
# Begin Source File

SOURCE=.\HigherOrderFactory.hpp
# End Source File
# Begin Source File

SOURCE=.\HomXform.hpp
# End Source File
# Begin Source File

SOURCE=.\MBBits.hpp
# End Source File
# Begin Source File

SOURCE=.\MBCN.hpp
# End Source File
# Begin Source File

SOURCE=.\MBCore.hpp
# End Source File
# Begin Source File

SOURCE=.\MBError.hpp
# End Source File
# Begin Source File

SOURCE=.\MBInterface.hpp
# End Source File
# Begin Source File

SOURCE=.\MBInternals.hpp
# End Source File
# Begin Source File

SOURCE=.\MBMem.hpp
# End Source File
# Begin Source File

SOURCE=.\MBMeshSet.hpp
# End Source File
# Begin Source File

SOURCE=.\MBRange.hpp
# End Source File
# Begin Source File

SOURCE=.\MBReadUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\MBReadUtilIface.hpp
# End Source File
# Begin Source File

SOURCE=.\MBSkinner.hpp
# End Source File
# Begin Source File

SOURCE=.\MBUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\MBWriteUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\MBWriteUtilIface.hpp
# End Source File
# Begin Source File

SOURCE=.\MeshTopoUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\PolyEntitySequence.hpp
# End Source File
# Begin Source File

SOURCE=.\ReadNCDF.hpp
# End Source File
# Begin Source File

SOURCE=.\ReadVtk.hpp
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

SOURCE=.\Tqdcfr.hpp
# End Source File
# Begin Source File

SOURCE=.\WriteGMV.hpp
# End Source File
# Begin Source File

SOURCE=.\WriteNCDF.hpp
# End Source File
# Begin Source File

SOURCE=.\WriteSLAC.hpp
# End Source File
# End Group
# End Target
# End Project
