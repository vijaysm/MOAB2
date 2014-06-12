set(moab_old_include_headers MBAdaptiveKDTree.hpp
                             MBBSPTree.hpp
                             MBBSPTreePoly.hpp
                             MBCartVect.hpp
                             MBCN.hpp
                             MBCore.hpp
                             MBEntityHandle.h
                             MBForward.hpp
                             MBGeomUtil.hpp
                             MBInterface.hpp
                             MBmpi.h
                             MBOrientedBoxTreeTool.hpp
                             MBParallelComm.hpp
                             MBParallelData.hpp
                             MBProcConfig.hpp
                             MBRange.hpp
                             MBReaderIface.hpp
                             MBReaderWriterSet.hpp
                             MBReadUtilIface.hpp
                             MBSkinner.hpp
                             MBTypes.h
                             MBUnknownInterface.hpp
                             MBUtil.hpp
                             MBVersion.h
                             MBWriterIface.hpp
                             MBWriteUtilIface.hpp )

include(AutoconfHeader)
moab_install_headers(${moab_old_include_headers})
