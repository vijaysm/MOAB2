# Microsoft Developer Studio Project File - Name="MDB" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=MDB - Win32 Release
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "MDB.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "MDB.mak" CFG="MDB - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "MDB - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "MDB - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "MDB - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "MDB___Win32_Release"
# PROP BASE Intermediate_Dir "MDB___Win32_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "MDB_EXPORTS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "%CUBITROOT%\exodus\exodus3.07\include" /I "%CUBITROOT%\netcdf\netcdf-3.4\include" /I "$(CUBITROOT)\exodus\exodus3.07\include" /I "$(CUBITROOT)\netcdf\netcdf-3.4\include" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "MDB_EXPORTS" /D "IS_BUILDING_MDB" /YX /FD /c
# SUBTRACT CPP /Fr
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 libexoIIv2c.lib netcdfs.lib /nologo /dll /machine:I386 /out:"components/MDB.dll" /libpath:"$(CUBITROOT)\exodus\exodus3.07\lib\NT" /libpath:"$(CUBITROOT)\netcdf\netcdf-3.4\lib"

!ELSEIF  "$(CFG)" == "MDB - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "MDB___Win32_Debug"
# PROP BASE Intermediate_Dir "MDB___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "MDB_EXPORTS" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "$(CUBITROOT)\exodus\exodus3.07\include" /I "$(CUBITROOT)\netcdf\netcdf-3.4\include" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "MDB_EXPORTS" /D "IS_BUILDING_MDB" /YX /FD /GZ /c
# SUBTRACT CPP /Fr
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 libexoIIv2c_db.lib netcdfsd.lib /nologo /dll /debug /machine:I386 /out:"components/MDBd.dll" /pdbtype:sept /libpath:"$(CUBITROOT)\exodus\exodus3.07\lib\NT" /libpath:"$(CUBITROOT)\netcdf\netcdf-3.4\lib"
# SUBTRACT LINK32 /profile

!ENDIF 

# Begin Target

# Name "MDB - Win32 Release"
# Name "MDB - Win32 Debug"
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

SOURCE=.\EntitySequence.hpp
# End Source File
# Begin Source File

SOURCE=.\EntitySequenceManager.hpp
# End Source File
# Begin Source File

SOURCE=.\ExoIIInterface.hpp
# End Source File
# Begin Source File

SOURCE=.\HigherOrderFactory.hpp
# End Source File
# Begin Source File

SOURCE=.\HomXform.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBCore.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBDefines.h
# End Source File
# Begin Source File

SOURCE=.\MDBEntitySequence.hpp
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

SOURCE=.\MDBSkinner.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBWriteUtil.hpp
# End Source File
# Begin Source File

SOURCE=.\MDBWriteUtilIface.hpp
# End Source File
# Begin Source File

SOURCE=.\MOABCN.hpp
# End Source File
# Begin Source File

SOURCE=.\ReadExoII.hpp
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

SOURCE=.\TSTT_Base.h
# End Source File
# Begin Source File

SOURCE=.\WriteExoII.hpp
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
