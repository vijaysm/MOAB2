include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_BINARY_DIR}
                    parallel
                    LocalDiscretization
                    io io/mhdf
                    ${MOAB_BINARY_DIR}/src/parallel )

if(MOAB_USE_MPI AND MPI_FOUND)
  add_subdirectory(parallel)
  set( MOAB_MPI moabpar )
  set( MPI_LIBRARIES ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY} )
endif()

#create the FCDefs.h and src/MOAB_FCDefs.h files
configure_file(FCDefs.h.in "${CMAKE_CURRENT_BINARY_DIR}/FCDefs.h" @ONLY)
#sed -e "s/FC_FUNC/MOAB_FC_FUNC/" src/FCDefs.h >src/MOAB_FCDefs.h
file(READ ${CMAKE_CURRENT_BINARY_DIR}/FCDefs.h infile)
string( REGEX REPLACE "FC_FUNC" "MOAB_FC_FUNC" outstring "${infile}" )
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/MOAB_FCDefs.h" "${outstring}")

#Create the moab/EntityHandle.hpp file
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/moab/EntityHandle.hpp.in infile2)
string( REGEX REPLACE "#undef" "#cmakedefine" outstring2 "${infile2}" )
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/moab/EntityHandle.hpp.in" "${outstring2}")
configure_file(${CMAKE_CURRENT_BINARY_DIR}/moab/EntityHandle.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/moab/EntityHandle.hpp)

#generate the Version.h file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/moab/Version.h.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/moab/Version.h" @ONLY)

set(moab_link_lib)
add_subdirectory(io)
add_subdirectory(LocalDiscretization)
add_subdirectory(oldinc)

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/MBCN_protos.h
  COMMAND ${CMAKE_COMMAND}
          -DDEF:STRING=MOAB
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/MBCN.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/MBCN_protos.h
          -Dprefix:STRING=MBCN
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/MBCN.h )

if( MOAB_USE_HDF )
  include_directories(${HDF5_INCLUDES})
  set( MOAB_HDF_LIBS ${HDF5_LIBRARIES} )
endif()
if( MOAB_USE_NETCDF )
  include_directories(${NetCDF_INCLUDES})
  if ( MOAB_USE_MPI )
    set( MOAB_NETCDF_LIBS ${PNetCDF_LIBRARIES} ${NetCDF_LIBRARIES} )
  else ( MOAB_USE_MPI )
    set( MOAB_NETCDF_LIBS ${NetCDF_LIBRARIES} )
  endif ( MOAB_USE_MPI )
endif()

set(MOAB_LIB_SOURCES
        AdaptiveKDTree.cpp
        AEntityFactory.hpp AEntityFactory.cpp
        AffineXform.hpp    AffineXform.cpp
        AxisBox.hpp        AxisBox.cpp
        BitPage.hpp        BitPage.cpp
        BitTag.hpp         BitTag.cpp
        BoundBox.cpp
        BSPTree.cpp
        BSPTreePoly.cpp
        BVHTree.cpp
        CN.cpp
        CartVect.cpp
        Core.cpp
        DebugOutput.hpp DebugOutput.cpp
        DenseTag.hpp    DenseTag.cpp
        DualTool.cpp
        ElementSequence.hpp
        EntitySequence.hpp  EntitySequence.cpp
        Factory.cpp
        FBEngine.cpp
        FileOptions.cpp
        GeomUtil.cpp
        GeomTopoTool.cpp
        HigherOrderFactory.cpp
        HomXform.cpp
        Internals.hpp
        MBCNArrays.hpp
        MergeMesh.cpp
        MeshSet.hpp         MeshSet.cpp
        MeshSetSequence.hpp MeshSetSequence.cpp
        MeshTag.hpp         MeshTag.cpp
        MeshTopoUtil.cpp
        OrientedBox.hpp     OrientedBox.cpp
        OrientedBoxTreeTool.cpp
        lotte/poly.c
        lotte/findpt.c
        lotte/errmem.c
        lotte/tensor.c
        PolyElementSeq.hpp  PolyElementSeq.cpp
        ProgOptions.cpp
        Range.cpp
        RangeSeqIntersectIter.hpp RangeSeqIntersectIter.cpp
        ReadUtil.hpp              ReadUtil.cpp
        ReaderWriterSet.cpp
        ReorderTool.cpp
        ScdElementData.hpp        ScdElementData.cpp
        ScdInterface.cpp
        ScdVertexData.hpp         ScdVertexData.cpp
        SequenceData.hpp          SequenceData.cpp
        SequenceManager.hpp       SequenceManager.cpp
        SetIterator.cpp
        Skinner.cpp
        SmoothCurve.hpp           SmoothCurve.cpp
        SmoothFace.hpp            SmoothFace.cpp
        SparseTag.hpp             SparseTag.cpp
        SpatialLocator.cpp
        SpectralMeshTool.cpp
        StructuredElementSeq.hpp  StructuredElementSeq.cpp
        SweptElementData.hpp      SweptElementData.cpp
        SweptElementSeq.hpp       SweptElementSeq.cpp
        SweptVertexData.hpp       SweptVertexData.cpp
        SysUtil.hpp               SysUtil.cpp
        TagCompare.hpp
        TagInfo.hpp               TagInfo.cpp
        TupleList.cpp
        Tree.cpp
        Types.cpp
        TypeSequenceManager.hpp   TypeSequenceManager.cpp
        UnstructuredElemSeq.hpp   UnstructuredElemSeq.cpp
        Util.cpp
        VarLenDenseTag.hpp        VarLenDenseTag.cpp
        VarLenSparseTag.hpp       VarLenSparseTag.cpp
        VarLenTag.hpp      
        VertexSequence.hpp        VertexSequence.cpp
        WriteUtil.hpp             WriteUtil.cpp
        moab_mpe.h )

set( MOAB_INSTALL_HEADERS
        moab/AdaptiveKDTree.hpp
        moab/BoundBox.hpp
        moab/BSPTree.hpp
        moab/BSPTreePoly.hpp
        moab/BVHTree.hpp
        moab/CN.hpp
        moab/CartVect.hpp
        moab/Compiler.hpp
        moab/Core.hpp
        moab/CpuTimer.hpp
        moab/DualTool.hpp
        moab/Error.hpp
        moab/GeomTopoTool.hpp
        moab/HigherOrderFactory.hpp
        moab/HomXform.hpp
        moab/EntityType.hpp
        moab/FBEngine.hpp
        moab/FileOptions.hpp
        moab/FindPtFuncs.h
        moab/Forward.hpp
        moab/GeomUtil.hpp
        moab/Interface.hpp
        moab/point_locater/tree/common_tree.hpp
        moab/point_locater/tree/element_tree.hpp
        moab/point_locater/tree/bvh_tree.hpp
        moab/point_locater/io.hpp 
        moab/point_locater/element_maps/linear_hex_map.hpp 
        moab/point_locater/element_maps/linear_tet_map.hpp 
        moab/point_locater/element_maps/spectral_hex_map.hpp 
        moab/point_locater/element_maps/quadratic_hex_map.hpp 
        moab/point_locater/point_locater.hpp
        moab/point_locater/parametrizer.hpp
        moab/Matrix3.hpp
        moab/MergeMesh.hpp
        moab/MeshTopoUtil.hpp
        moab/OrientedBoxTreeTool.hpp
        moab/ProgOptions.hpp
        moab/Range.hpp
        moab/RangeMap.hpp
        moab/ReadUtilIface.hpp
        moab/ReaderIface.hpp
        moab/ReaderWriterSet.hpp
        moab/ReorderTool.hpp
        moab/ScdInterface.hpp
        moab/SetIterator.hpp
        moab/Skinner.hpp
        moab/SpatialLocator.hpp
        moab/SpectralMeshTool.hpp
        moab/Tree.hpp
        moab/TreeStats.hpp
        moab/Types.hpp
        moab/UnknownInterface.hpp
        moab/Util.hpp
        moab/WriteUtilIface.hpp
        moab/WriterIface.hpp
        MBEntityType.h
        MBCN.h
        MBTagConventions.hpp )

set( MOAB_GENERATED_INSTALL_HEADERS
     ${CMAKE_CURRENT_BINARY_DIR}/MOAB_FCDefs.h
     ${CMAKE_CURRENT_BINARY_DIR}/MBCN_protos.h
     ${CMAKE_CURRENT_BINARY_DIR}/moab/EntityHandle.hpp )

if( CGM_FOUND )
  set( MOAB_CGM_LIBS ${CGM_LIBRARIES} )
endif()

if (MOAB_USE_ZOLTAN)
  set( MOAB_LIB_SOURCES ${MOAB_LIB_SOURCES}
          ZoltanPartitioner.cpp )
  set( MOAB_INSTALL_HEADERS ${MOAB_INSTALL_HEADERS}
          moab/ZoltanPartitioner.hpp )
  include_directories( ${Zoltan_INCLUDES} )
  list(APPEND MOAB_DEP_LIBRARIES ${Zoltan_LIBRARIES})
endif()
if(MOAB_USE_MPI AND MPI_FOUND)
  add_library( MOAB ${MOAB_LIB_SOURCES} ${MOAB_INSTALL_HEADERS} ${MOAB_GENERATED_INSTALL_HEADERS}
                    $<TARGET_OBJECTS:moabio>
                    $<TARGET_OBJECTS:moab_loc_discr>
                    $<TARGET_OBJECTS:moabpar>
                    ${MDHF_OBJ} )
else(MOAB_USE_MPI AND MPI_FOUND)
  add_library( MOAB ${MOAB_LIB_SOURCES} ${MOAB_INSTALL_HEADERS} ${MOAB_GENERATED_INSTALL_HEADERS}
                    $<TARGET_OBJECTS:moabio>
                    $<TARGET_OBJECTS:moab_loc_discr>
                    ${MDHF_OBJ} )
endif(MOAB_USE_MPI AND MPI_FOUND)

target_link_libraries( MOAB ${MOAB_DEP_LIBRARIES} ${MOAB_LIBS} ${MOAB_CGM_LIBS} ${MPI_LIBRARIES} )
set_target_properties( MOAB PROPERTIES COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}")

install(TARGETS MOAB EXPORT MOABTargets
                     RUNTIME DESTINATION bin
                     LIBRARY DESTINATION lib
                     ARCHIVE DESTINATION lib )
set( MOAB_INSTALL_TARGETS ${MOAB_INSTALL_TARGETS} MOAB)
include(AutoconfHeader)
moab_install_headers(${MOAB_INSTALL_HEADERS})
foreach (header ${MOAB_GENERATED_INSTALL_HEADERS})
    STRING(REGEX REPLACE "^${CMAKE_CURRENT_BINARY_DIR}/" "" relative ${header})
    string(REGEX MATCH "(.*)/" DIR ${relative})
    install(FILES ${header} DESTINATION include/${DIR})
endforeach()
