SRC_MESH_DIR = Mesh
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_MESH_DIR)
INCLUDED_MAKEFILES += $(SRC_MESH_DIR)/Makefile.am

if HAVE_MPI
  MPI_SRCS = $(SRC_MESH_DIR)/ParallelHelper.cpp
  MPI_INCS = $(SRC_MESH_DIR)/ParallelHelper.hpp
else
  MPI_SRCS = $(SRC_MESH_DIR)/ParallelHelperNoMPIStub.cpp
  MPI_INCS = $(SRC_MESH_DIR)/ParallelHelperNoMPIStub.hpp
endif


MSQ_SRCS += \
	$(SRC_MESH_DIR)/ArrayMesh.cpp \
	$(SRC_MESH_DIR)/ElementPatches.cpp \
	$(SRC_MESH_DIR)/ExtraData.cpp \
	$(SRC_MESH_DIR)/GlobalPatch.cpp \
	$(SRC_MESH_DIR)/IdealElements.cpp \
	$(SRC_MESH_DIR)/MeshDecorator.cpp \
	$(SRC_MESH_DIR)/MeshImpl.cpp \
	$(SRC_MESH_DIR)/MeshImplData.cpp \
	$(SRC_MESH_DIR)/MeshImplTags.cpp \
	$(SRC_MESH_DIR)/MeshUtil.cpp \
	$(SRC_MESH_DIR)/MeshWriter.cpp \
	$(SRC_MESH_DIR)/MsqMeshEntity.cpp \
	$(SRC_MESH_DIR)/MsqVertex.cpp \
	$(SRC_MESH_DIR)/ParallelMeshImpl.cpp \
	$(SRC_MESH_DIR)/PatchData.cpp \
	$(SRC_MESH_DIR)/PatchIterator.cpp \
	$(SRC_MESH_DIR)/PatchSet.cpp \
	$(SRC_MESH_DIR)/TagVertexMesh.cpp \
	$(SRC_MESH_DIR)/TopologyInfo.cpp \
	$(SRC_MESH_DIR)/VertexPatches.cpp \
	$(MPI_SRCS) 

MSQ_HDRS += \
	$(SRC_MESH_DIR)/ArrayMesh.hpp \
	$(SRC_MESH_DIR)/ElementPatches.hpp \
	$(SRC_MESH_DIR)/ExtraData.hpp \
	$(SRC_MESH_DIR)/ExtraDataUser.hpp \
	$(SRC_MESH_DIR)/GlobalPatch.hpp \
	$(SRC_MESH_DIR)/IdealElements.hpp \
	$(SRC_MESH_DIR)/MeshDecorator.hpp \
	$(SRC_MESH_DIR)/MeshImpl.hpp \
	$(SRC_MESH_DIR)/MeshImplData.hpp \
	$(SRC_MESH_DIR)/MeshImplTags.hpp \
	$(SRC_MESH_DIR)/MeshUtil.hpp \
	$(SRC_MESH_DIR)/MeshWriter.hpp \
	$(SRC_MESH_DIR)/MsqFreeVertexIndexIterator.hpp \
	$(SRC_MESH_DIR)/MsqMeshEntity.hpp \
	$(SRC_MESH_DIR)/MsqVertex.hpp \
	$(SRC_MESH_DIR)/ParallelMeshImpl.hpp \
	$(SRC_MESH_DIR)/PatchData.hpp \
	$(SRC_MESH_DIR)/PatchIterator.hpp \
	$(SRC_MESH_DIR)/PatchSet.hpp \
	$(SRC_MESH_DIR)/TagVertexMesh.hpp \
	$(SRC_MESH_DIR)/TopologyInfo.hpp \
	$(SRC_MESH_DIR)/VertexPatches.hpp \
	$(MPI_INCS) 

SRC_CONTROL_DIR = Control
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_CONTROL_DIR)
INCLUDED_MAKEFILES += $(SRC_CONTROL_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_CONTROL_DIR)/Instruction.cpp \
	$(SRC_CONTROL_DIR)/InstructionQueue.cpp \
	$(SRC_CONTROL_DIR)/IQInterface.cpp \
	$(SRC_CONTROL_DIR)/TerminationCriterion.cpp \
        $(SRC_CONTROL_DIR)/Settings.cpp \
        $(SRC_CONTROL_DIR)/SlaveBoundaryVertices.cpp 

MSQ_HDRS += \
	$(SRC_CONTROL_DIR)/Instruction.hpp \
	$(SRC_CONTROL_DIR)/InstructionQueue.hpp \
	$(SRC_CONTROL_DIR)/IQInterface.hpp \
	$(SRC_CONTROL_DIR)/TerminationCriterion.hpp \
        $(SRC_CONTROL_DIR)/Settings.hpp  \
        $(SRC_CONTROL_DIR)/SlaveBoundaryVertices.hpp \
        $(SRC_CONTROL_DIR)/VertexSlaver.hpp 

SRC_WRAPPERS_DIR = Wrappers
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_WRAPPERS_DIR)
INCLUDED_MAKEFILES += $(SRC_WRAPPERS_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_WRAPPERS_DIR)/DeformingDomainWrapper.cpp \
	$(SRC_WRAPPERS_DIR)/LaplaceWrapper.cpp \
	$(SRC_WRAPPERS_DIR)/ShapeImprovementWrapper.cpp \
	$(SRC_WRAPPERS_DIR)/ShapeImprover.cpp \
	$(SRC_WRAPPERS_DIR)/SizeAdaptShapeWrapper.cpp \
	$(SRC_WRAPPERS_DIR)/PaverMinEdgeLengthWrapper.cpp \
	$(SRC_WRAPPERS_DIR)/UntangleWrapper.cpp \
	$(SRC_WRAPPERS_DIR)/ViscousCFDTetShapeWrapper.cpp \
	$(SRC_WRAPPERS_DIR)/Wrapper.cpp

MSQ_HDRS += \
	$(SRC_WRAPPERS_DIR)/DeformingDomainWrapper.hpp \
	$(SRC_WRAPPERS_DIR)/LaplaceWrapper.hpp \
	$(SRC_WRAPPERS_DIR)/ShapeImprovementWrapper.hpp \
	$(SRC_WRAPPERS_DIR)/ShapeImprover.hpp \
	$(SRC_WRAPPERS_DIR)/SizeAdaptShapeWrapper.hpp \
	$(SRC_WRAPPERS_DIR)/PaverMinEdgeLengthWrapper.hpp \
	$(SRC_WRAPPERS_DIR)/UntangleWrapper.hpp \
	$(SRC_WRAPPERS_DIR)/ViscousCFDTetShapeWrapper.hpp \
	$(SRC_WRAPPERS_DIR)/Wrapper.hpp
SRC_MAPPINGFUNCTION_DIR = MappingFunction
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_MAPPINGFUNCTION_DIR)
INCLUDED_MAKEFILES += $(SRC_MAPPINGFUNCTION_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_MAPPINGFUNCTION_DIR)/MappingFunction.cpp \
	$(SRC_MAPPINGFUNCTION_DIR)/NodeSet.cpp \
        $(SRC_MAPPINGFUNCTION_DIR)/Sample.cpp

MSQ_HDRS += \
	$(SRC_MAPPINGFUNCTION_DIR)/MappingFunction.hpp \
	$(SRC_MAPPINGFUNCTION_DIR)/NodeSet.hpp  \
        $(SRC_MAPPINGFUNCTION_DIR)/Sample.hpp

SRC_MAPPINGFUNCTION_LAGRANGE_DIR = MappingFunction/Lagrange
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)
INCLUDED_MAKEFILES += $(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)/HexLagrangeShape.cpp \
	$(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)/QuadLagrangeShape.cpp \
	$(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)/TetLagrangeShape.cpp \
	$(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)/TriLagrangeShape.cpp 

MSQ_HDRS += \
	$(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)/HexLagrangeShape.hpp \
	$(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)/QuadLagrangeShape.hpp \
	$(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)/TetLagrangeShape.hpp \
	$(SRC_MAPPINGFUNCTION_LAGRANGE_DIR)/TriLagrangeShape.hpp 

SRC_MAPPINGFUNCTION_LINEAR_DIR = MappingFunction/Linear
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_MAPPINGFUNCTION_LINEAR_DIR)
INCLUDED_MAKEFILES += $(SRC_MAPPINGFUNCTION_LINEAR_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearHexahedron.cpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearPrism.cpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearPyramid.cpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearQuadrilateral.cpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearTetrahedron.cpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearTriangle.cpp 

MSQ_HDRS += \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearHexahedron.hpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearPrism.hpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearPyramid.hpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearQuadrilateral.hpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearTetrahedron.hpp \
	$(SRC_MAPPINGFUNCTION_LINEAR_DIR)/LinearTriangle.hpp 

SRC_MISC_DIR = Misc
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_MISC_DIR)
INCLUDED_MAKEFILES += $(SRC_MISC_DIR)/Makefile.am

if ENABLE_imesh
  MSQ_ITAPS_COMMON_SRC = $(SRC_MISC_DIR)/MsqIBase.cpp \
                         $(SRC_MISC_DIR)/MsqIBase.hpp
else
  MSQ_ITAPS_COMMON_SRC =
endif

MSQ_SRCS += \
	$(SRC_MISC_DIR)/BoundedCylinderDomain.cpp \
	$(SRC_MISC_DIR)/ConicDomain.cpp \
	$(SRC_MISC_DIR)/CylinderDomain.cpp \
	$(SRC_MISC_DIR)/DomainClassifier.cpp \
	$(SRC_MISC_DIR)/DomainUtil.cpp \
	$(SRC_MISC_DIR)/DomainUtil.hpp \
	$(SRC_MISC_DIR)/EdgeIterator.cpp \
	$(SRC_MISC_DIR)/Exponent.cpp \
	$(SRC_MISC_DIR)/FileTokenizer.cpp \
	$(SRC_MISC_DIR)/MeshDomain1D.cpp \
	$(SRC_MISC_DIR)/MeshTransform.cpp \
	$(SRC_MISC_DIR)/MesquiteVersion.cpp \
	$(SRC_MISC_DIR)/MsqDebug.cpp \
	$(SRC_MISC_DIR)/MsqError.cpp \
	$(SRC_MISC_DIR)/MsqFPE.cpp \
	$(SRC_MISC_DIR)/MsqGeomPrim.cpp \
	$(SRC_MISC_DIR)/MsqHessian.cpp \
	$(SRC_MISC_DIR)/MsqInterrupt.cpp \
	$(SRC_MISC_DIR)/MsqTimer.cpp \
	$(SRC_MISC_DIR)/PlanarDomain.cpp \
	$(SRC_MISC_DIR)/SimpleStats.cpp \
	$(SRC_MISC_DIR)/SphericalDomain.cpp \
	$(SRC_MISC_DIR)/Vector3D.cpp \
	$(SRC_MISC_DIR)/VtkTypeInfo.cpp \
	$(SRC_MISC_DIR)/XYPlanarDomain.cpp \
	$(SRC_MISC_DIR)/XYRectangle.cpp \
        $(MSQ_ITAPS_COMMON_SRC)


MSQ_HDRS += \
	$(SRC_MISC_DIR)/BoundedCylinderDomain.hpp \
	$(SRC_MISC_DIR)/Bits.hpp \
	$(SRC_MISC_DIR)/ConicDomain.hpp \
	$(SRC_MISC_DIR)/CylinderDomain.hpp \
	$(SRC_MISC_DIR)/DomainClassifier.hpp \
	$(SRC_MISC_DIR)/EdgeIterator.hpp \
	$(SRC_MISC_DIR)/Exponent.hpp \
	$(SRC_MISC_DIR)/FileTokenizer.hpp \
	$(SRC_MISC_DIR)/Matrix3D.hpp \
	$(SRC_MISC_DIR)/MeshDomain1D.hpp \
	$(SRC_MISC_DIR)/MeshTransform.hpp \
	$(SRC_MISC_DIR)/MsqDebug.hpp \
	$(SRC_MISC_DIR)/MsqError.hpp \
	$(SRC_MISC_DIR)/MsqFPE.hpp \
	$(SRC_MISC_DIR)/MsqGeomPrim.hpp \
	$(SRC_MISC_DIR)/MsqHessian.hpp \
	$(SRC_MISC_DIR)/MsqInterrupt.hpp \
	$(SRC_MISC_DIR)/MsqMatrix.hpp \
	$(SRC_MISC_DIR)/MsqTimer.hpp \
	$(SRC_MISC_DIR)/PlanarDomain.hpp \
	$(SRC_MISC_DIR)/SimpleStats.hpp \
	$(SRC_MISC_DIR)/SphericalDomain.hpp \
	$(SRC_MISC_DIR)/SymMatrix3D.hpp \
	$(SRC_MISC_DIR)/Vector3D.hpp \
	$(SRC_MISC_DIR)/VtkTypeInfo.hpp \
	$(SRC_MISC_DIR)/XYPlanarDomain.hpp \
	$(SRC_MISC_DIR)/XYRectangle.hpp

SRC_OBJECTIVEFUNCTION_DIR = ObjectiveFunction
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_OBJECTIVEFUNCTION_DIR)
INCLUDED_MAKEFILES += $(SRC_OBJECTIVEFUNCTION_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_OBJECTIVEFUNCTION_DIR)/CompositeOFAdd.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/CompositeOFMultiply.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/CompositeOFScalarAdd.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/CompositeOFScalarMultiply.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/LInfTemplate.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/LPtoPTemplate.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/MaxTemplate.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/OFEvaluator.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/ObjectiveFunction.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/ObjectiveFunctionTemplate.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/PatchPowerMeanP.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/PMeanPTemplate.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/StdDevTemplate.cpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/VarianceTemplate.cpp

MSQ_HDRS += \
	$(SRC_OBJECTIVEFUNCTION_DIR)/CompositeOFAdd.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/CompositeOFMultiply.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/CompositeOFScalarAdd.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/CompositeOFScalarMultiply.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/LInfTemplate.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/LPtoPTemplate.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/MaxTemplate.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/OFEvaluator.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/ObjectiveFunction.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/ObjectiveFunctionTemplate.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/PatchPowerMeanP.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/PMeanPTemplate.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/StdDevTemplate.hpp \
	$(SRC_OBJECTIVEFUNCTION_DIR)/VarianceTemplate.hpp

SRC_QUALITYASSESSOR_DIR = QualityAssessor
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_QUALITYASSESSOR_DIR)
INCLUDED_MAKEFILES += $(SRC_QUALITYASSESSOR_DIR)/Makefile.am

MSQ_SRCS += $(SRC_QUALITYASSESSOR_DIR)/QualityAssessor.cpp 

MSQ_HDRS += $(SRC_QUALITYASSESSOR_DIR)/QualityAssessor.hpp 

SRC_QUALITYIMPROVER_DIR = QualityImprover
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_QUALITYIMPROVER_DIR)
INCLUDED_MAKEFILES += $(SRC_QUALITYIMPROVER_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_QUALITYIMPROVER_DIR)/PatchSetUser.cpp \
	$(SRC_QUALITYIMPROVER_DIR)/QualityImprover.cpp \
	$(SRC_QUALITYIMPROVER_DIR)/VertexMover.cpp

MSQ_HDRS += \
	$(SRC_QUALITYIMPROVER_DIR)/NullImprover.hpp \
	$(SRC_QUALITYIMPROVER_DIR)/PatchSetUser.hpp \
	$(SRC_QUALITYIMPROVER_DIR)/QualityImprover.hpp \
	$(SRC_QUALITYIMPROVER_DIR)/VertexMover.hpp 

SRC_OPTSOLVER_DIR = QualityImprover/OptSolvers
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_OPTSOLVER_DIR)
INCLUDED_MAKEFILES += $(SRC_OPTSOLVER_DIR)/Makefile.am

MSQ_SRCS += $(SRC_OPTSOLVER_DIR)/ConjugateGradient.cpp \
	    $(SRC_OPTSOLVER_DIR)/FeasibleNewton.cpp \
	    $(SRC_OPTSOLVER_DIR)/NonSmoothDescent.cpp \
	    $(SRC_OPTSOLVER_DIR)/QuasiNewton.cpp \
	    $(SRC_OPTSOLVER_DIR)/SteepestDescent.cpp \
	    $(SRC_OPTSOLVER_DIR)/NonGradient.cpp \
	    $(SRC_OPTSOLVER_DIR)/TrustRegion.cpp

MSQ_HDRS += $(SRC_OPTSOLVER_DIR)/ConjugateGradient.hpp \
	    $(SRC_OPTSOLVER_DIR)/FeasibleNewton.hpp \
	    $(SRC_OPTSOLVER_DIR)/NonSmoothDescent.hpp \
	    $(SRC_OPTSOLVER_DIR)/QuasiNewton.hpp \
	    $(SRC_OPTSOLVER_DIR)/SteepestDescent.hpp \
	    $(SRC_OPTSOLVER_DIR)/NonGradient.hpp \
	    $(SRC_OPTSOLVER_DIR)/TrustRegion.hpp 
SRC_RELAXATIONSMOOTHERS_DIR = QualityImprover/Relaxation
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_RELAXATIONSMOOTHERS_DIR)
INCLUDED_MAKEFILES += $(SRC_RELAXATIONSMOOTHERS_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_RELAXATIONSMOOTHERS_DIR)/Randomize.cpp \
	$(SRC_RELAXATIONSMOOTHERS_DIR)/RelaxationSmoother.cpp \
	$(SRC_RELAXATIONSMOOTHERS_DIR)/LaplacianSmoother.cpp \
	$(SRC_RELAXATIONSMOOTHERS_DIR)/SmartLaplacianSmoother.cpp 

MSQ_HDRS += \
	$(SRC_RELAXATIONSMOOTHERS_DIR)/Randomize.hpp \
	$(SRC_RELAXATIONSMOOTHERS_DIR)/RelaxationSmoother.hpp \
	$(SRC_RELAXATIONSMOOTHERS_DIR)/LaplacianSmoother.hpp \
	$(SRC_RELAXATIONSMOOTHERS_DIR)/SmartLaplacianSmoother.hpp 

SRC_QUALITYMETRIC_DIR = QualityMetric
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_QUALITYMETRIC_DIR) \
                -I$(top_srcdir)/src/$(SRC_QUALITYMETRIC_DIR)/Debug
INCLUDED_MAKEFILES += $(SRC_QUALITYMETRIC_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_QUALITYMETRIC_DIR)/Debug/CompareQM.cpp \
	$(SRC_QUALITYMETRIC_DIR)/Debug/NumericalQM.cpp \
	$(SRC_QUALITYMETRIC_DIR)/AddQualityMetric.cpp \
	$(SRC_QUALITYMETRIC_DIR)/AveragingQM.cpp \
	$(SRC_QUALITYMETRIC_DIR)/EdgeQM.cpp \
	$(SRC_QUALITYMETRIC_DIR)/ElemSampleQM.cpp \
	$(SRC_QUALITYMETRIC_DIR)/ElementMaxQM.cpp \
	$(SRC_QUALITYMETRIC_DIR)/ElementAvgQM.cpp \
	$(SRC_QUALITYMETRIC_DIR)/ElementPMeanP.cpp \
	$(SRC_QUALITYMETRIC_DIR)/ElementQM.cpp \
	$(SRC_QUALITYMETRIC_DIR)/MultiplyQualityMetric.cpp \
	$(SRC_QUALITYMETRIC_DIR)/PMeanPMetric.cpp \
	$(SRC_QUALITYMETRIC_DIR)/PowerQualityMetric.cpp \
	$(SRC_QUALITYMETRIC_DIR)/QualityMetric.cpp \
	$(SRC_QUALITYMETRIC_DIR)/ScalarAddQualityMetric.cpp \
	$(SRC_QUALITYMETRIC_DIR)/ScalarMultiplyQualityMetric.cpp \
	$(SRC_QUALITYMETRIC_DIR)/VertexMaxQM.cpp \
	$(SRC_QUALITYMETRIC_DIR)/VertexPMeanP.cpp \
	$(SRC_QUALITYMETRIC_DIR)/VertexQM.cpp 

MSQ_HDRS += \
	$(SRC_QUALITYMETRIC_DIR)/Debug/CompareQM.hpp \
	$(SRC_QUALITYMETRIC_DIR)/Debug/NumericalQM.hpp \
	$(SRC_QUALITYMETRIC_DIR)/AddQualityMetric.hpp \
	$(SRC_QUALITYMETRIC_DIR)/AveragingQM.hpp \
	$(SRC_QUALITYMETRIC_DIR)/EdgeQM.hpp \
	$(SRC_QUALITYMETRIC_DIR)/ElemSampleQM.hpp \
	$(SRC_QUALITYMETRIC_DIR)/ElementMaxQM.hpp \
	$(SRC_QUALITYMETRIC_DIR)/ElementAvgQM.hpp \
	$(SRC_QUALITYMETRIC_DIR)/ElementPMeanP.hpp \
	$(SRC_QUALITYMETRIC_DIR)/ElementQM.hpp \
	$(SRC_QUALITYMETRIC_DIR)/MultiplyQualityMetric.hpp \
	$(SRC_QUALITYMETRIC_DIR)/PMeanPMetric.hpp \
	$(SRC_QUALITYMETRIC_DIR)/PowerQualityMetric.hpp \
	$(SRC_QUALITYMETRIC_DIR)/QualityMetric.hpp \
	$(SRC_QUALITYMETRIC_DIR)/ScalarAddQualityMetric.hpp \
	$(SRC_QUALITYMETRIC_DIR)/ScalarMultiplyQualityMetric.hpp \
	$(SRC_QUALITYMETRIC_DIR)/VertexMaxQM.hpp \
	$(SRC_QUALITYMETRIC_DIR)/VertexPMeanP.hpp \
	$(SRC_QUALITYMETRIC_DIR)/VertexQM.hpp 

SRC_QUALITYMETRIC_SHAPE_DIR = QualityMetric/Shape
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_QUALITYMETRIC_SHAPE_DIR)
INCLUDED_MAKEFILES += $(SRC_QUALITYMETRIC_SHAPE_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/AspectRatioGammaQualityMetric.cpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/ConditionNumberQualityMetric.cpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/IdealWeightInverseMeanRatio.cpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/IdealWeightMeanRatio.cpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/VertexConditionNumberQualityMetric.cpp 

MSQ_HDRS += \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/AspectRatioGammaQualityMetric.hpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/ConditionNumberFunctions.hpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/ConditionNumberQualityMetric.hpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/IdealWeightInverseMeanRatio.hpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/IdealWeightMeanRatio.hpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/MeanRatioFunctions.hpp \
	$(SRC_QUALITYMETRIC_SHAPE_DIR)/VertexConditionNumberQualityMetric.hpp 

SRC_QUALITYMETRIC_SMOOTHNESS_DIR = QualityMetric/Smoothness
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_QUALITYMETRIC_SMOOTHNESS_DIR)
INCLUDED_MAKEFILES += $(SRC_QUALITYMETRIC_SMOOTHNESS_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_QUALITYMETRIC_SMOOTHNESS_DIR)/EdgeLengthQualityMetric.cpp \
	$(SRC_QUALITYMETRIC_SMOOTHNESS_DIR)/EdgeLengthRangeQualityMetric.cpp 

MSQ_HDRS += \
	$(SRC_QUALITYMETRIC_SMOOTHNESS_DIR)/EdgeLengthQualityMetric.hpp \
	$(SRC_QUALITYMETRIC_SMOOTHNESS_DIR)/EdgeLengthRangeQualityMetric.hpp 

SRC_QUALITYMETRIC_TMP_DIR = QualityMetric/TMP
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_QUALITYMETRIC_TMP_DIR)
INCLUDED_MAKEFILES += $(SRC_QUALITYMETRIC_TMP_DIR)/Makefile.am

MSQ_SRCS += $(SRC_QUALITYMETRIC_TMP_DIR)/AffineMapMetric.cpp \
            $(SRC_QUALITYMETRIC_TMP_DIR)/TargetMetricUtil.cpp \
            $(SRC_QUALITYMETRIC_TMP_DIR)/TargetMetricUtil.hpp \
            $(SRC_QUALITYMETRIC_TMP_DIR)/AWQualityMetric.cpp \
            $(SRC_QUALITYMETRIC_TMP_DIR)/TQualityMetric.cpp \
            $(SRC_QUALITYMETRIC_TMP_DIR)/TMPQualityMetric.cpp 

MSQ_HDRS += $(SRC_QUALITYMETRIC_TMP_DIR)/AffineMapMetric.hpp \
            $(SRC_QUALITYMETRIC_TMP_DIR)/AWQualityMetric.hpp \
            $(SRC_QUALITYMETRIC_TMP_DIR)/TQualityMetric.hpp \
            $(SRC_QUALITYMETRIC_TMP_DIR)/TMPQualityMetric.hpp

SRC_QUALITYMETRIC_UNTANGLE_DIR = QualityMetric/Untangle
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_QUALITYMETRIC_UNTANGLE_DIR)
INCLUDED_MAKEFILES += $(SRC_QUALITYMETRIC_UNTANGLE_DIR)/Makefile.am

MSQ_SRCS += $(SRC_QUALITYMETRIC_UNTANGLE_DIR)/UntangleBetaQualityMetric.cpp 

MSQ_HDRS += $(SRC_QUALITYMETRIC_UNTANGLE_DIR)/UntangleBetaQualityMetric.hpp 

SRC_QUALITYMETRIC_VOLUME_DIR = QualityMetric/Volume
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_QUALITYMETRIC_VOLUME_DIR)
INCLUDED_MAKEFILES += $(SRC_QUALITYMETRIC_VOLUME_DIR)/Makefile.am

MSQ_SRCS += $(SRC_QUALITYMETRIC_VOLUME_DIR)/EdgeLengthMetric.cpp \
	    $(SRC_QUALITYMETRIC_VOLUME_DIR)/LocalSizeQualityMetric.cpp \
	    $(SRC_QUALITYMETRIC_VOLUME_DIR)/SizeMetric.cpp

MSQ_HDRS += $(SRC_QUALITYMETRIC_VOLUME_DIR)/EdgeLengthMetric.hpp \
	    $(SRC_QUALITYMETRIC_VOLUME_DIR)/LocalSizeQualityMetric.hpp \
	    $(SRC_QUALITYMETRIC_VOLUME_DIR)/SizeMetric.hpp

SRC_TARGETCALCULATOR_DIR = TargetCalculator
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_TARGETCALCULATOR_DIR)
INCLUDED_MAKEFILES += $(SRC_TARGETCALCULATOR_DIR)/Makefile.am

MSQ_SRCS += \
	$(SRC_TARGETCALCULATOR_DIR)/CachingTargetCalculator.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/IdealShapeTarget.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/InverseMetricWeight.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/JacobianCalculator.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/LambdaTarget.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/LambdaConstant.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/LVQDTargetCalculator.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/MetricWeight.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/RefMeshTargetCalculator.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/ReferenceMesh.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/RefSizeTargetCalculator.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/RemainingWeight.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/TargetCalculator.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/TargetReader.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/TargetWriter.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/TetDihedralWeight.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/WeightCalculator.cpp \
	$(SRC_TARGETCALCULATOR_DIR)/WeightReader.cpp 

MSQ_HDRS += \
	$(SRC_TARGETCALCULATOR_DIR)/CachingTargetCalculator.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/IdealShapeTarget.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/InverseMetricWeight.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/JacobianCalculator.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/LambdaTarget.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/LambdaConstant.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/LVQDTargetCalculator.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/MetricWeight.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/RefMeshTargetCalculator.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/ReferenceMesh.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/RefSizeTargetCalculator.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/RemainingWeight.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/TargetCalculator.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/TargetReader.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/TargetWriter.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/TetDihedralWeight.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/WeightCalculator.hpp \
	$(SRC_TARGETCALCULATOR_DIR)/WeightReader.hpp 

SRC_TARGETMETRIC_DIR = TargetMetric
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_TARGETMETRIC_DIR)
INCLUDED_MAKEFILES += $(SRC_TARGETMETRIC_DIR)/Makefile.am

MSQ_SRCS += $(SRC_TARGETMETRIC_DIR)/AWMetric.cpp \
            $(SRC_TARGETMETRIC_DIR)/TMetric.cpp \
	    $(SRC_TARGETMETRIC_DIR)/TMPCommon.hpp \
	    $(SRC_TARGETMETRIC_DIR)/TMPDerivs.hpp \
            $(SRC_TARGETMETRIC_DIR)/TMetricBarrier.cpp \
            $(SRC_TARGETMETRIC_DIR)/TMetricNonBarrier.cpp \
            $(SRC_TARGETMETRIC_DIR)/AWMetricBarrier.cpp \
            $(SRC_TARGETMETRIC_DIR)/AWMetricNonBarrier.cpp 

MSQ_HDRS += $(SRC_TARGETMETRIC_DIR)/AWMetric.hpp \
            $(SRC_TARGETMETRIC_DIR)/TMetric.hpp \
            $(SRC_TARGETMETRIC_DIR)/TMetricBarrier.hpp \
            $(SRC_TARGETMETRIC_DIR)/TMetricNonBarrier.hpp \
            $(SRC_TARGETMETRIC_DIR)/AWMetricBarrier.hpp \
            $(SRC_TARGETMETRIC_DIR)/AWMetricNonBarrier.hpp 
	    

SRC_TARGETMETRIC_MISC_DIR = TargetMetric/Misc
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_TARGETMETRIC_MISC_DIR)
INCLUDED_MAKEFILES += $(SRC_TARGETMETRIC_MISC_DIR)/Makefile.am

MSQ_SRCS += $(SRC_TARGETMETRIC_MISC_DIR)/InvTransBarrier.cpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TOffset.cpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TMixed.cpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TPower2.cpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TScale.cpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TSquared.cpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TSum.cpp

MSQ_HDRS += $(SRC_TARGETMETRIC_MISC_DIR)/InvTransBarrier.hpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TOffset.hpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TMixed.hpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TPower2.hpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TScale.hpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TSquared.hpp \
	    $(SRC_TARGETMETRIC_MISC_DIR)/TSum.hpp 
SRC_TARGETMETRIC_SHAPE_DIR = TargetMetric/Shape
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_TARGETMETRIC_SHAPE_DIR)
INCLUDED_MAKEFILES += $(SRC_TARGETMETRIC_SHAPE_DIR)/Makefile.am

MSQ_SRCS += $(SRC_TARGETMETRIC_SHAPE_DIR)/AWShape2DB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/AWShape2DNB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/AWShape2DNB2.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TInverseMeanRatio.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TShape2DNB2.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TShape3DB2.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TShapeB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TShapeNB1.cpp 

MSQ_HDRS += $(SRC_TARGETMETRIC_SHAPE_DIR)/AWShape2DB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/AWShape2DNB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/AWShape2DNB2.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TInverseMeanRatio.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TShape2DNB2.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TShape3DB2.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TShapeB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_DIR)/TShapeNB1.hpp 
SRC_TARGETMETRIC_SHAPE_ORIENT_DIR = TargetMetric/ShapeOrient
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)
INCLUDED_MAKEFILES += $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/Makefile.am

MSQ_SRCS += $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/AWShapeOrientNB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/TShapeOrientB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/TShapeOrientB2.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/TShapeOrientNB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/TShapeOrientNB2.cpp

MSQ_HDRS += $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/AWShapeOrientNB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/TShapeOrientB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/TShapeOrientB2.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/TShapeOrientNB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPE_ORIENT_DIR)/TShapeOrientNB2.hpp
SRC_TARGETMETRIC_SHAPESIZE_DIR = TargetMetric/ShapeSize
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_TARGETMETRIC_SHAPESIZE_DIR)
INCLUDED_MAKEFILES += $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/Makefile.am

MSQ_SRCS += $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/AWShapeSizeB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize2DB2.cpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize2DNB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize2DNB2.cpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize3DB2.cpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize3DB4.cpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize3DNB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSizeB1.cpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSizeB3.cpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSizeNB3.cpp

MSQ_HDRS += $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/AWShapeSizeB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize2DB2.hpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize2DNB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize2DNB2.hpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize3DB2.hpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize3DB4.hpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSize3DNB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSizeB1.hpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSizeB3.hpp \
	    $(SRC_TARGETMETRIC_SHAPESIZE_DIR)/TShapeSizeNB3.hpp

SRC_TARGETMETRIC_SSO_DIR = TargetMetric/ShapeSizeOrient
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_TARGETMETRIC_SSO_DIR)
INCLUDED_MAKEFILES += $(SRC_TARGETMETRIC_SSO_DIR)/Makefile.am

MSQ_SRCS += $(SRC_TARGETMETRIC_SSO_DIR)/AWShapeSizeOrientNB1.cpp \
            $(SRC_TARGETMETRIC_SSO_DIR)/TShapeSizeOrientB1.cpp \
            $(SRC_TARGETMETRIC_SSO_DIR)/TShapeSizeOrientB2.cpp \
            $(SRC_TARGETMETRIC_SSO_DIR)/TShapeSizeOrientNB1.cpp
            
MSQ_HDRS += $(SRC_TARGETMETRIC_SSO_DIR)/AWShapeSizeOrientNB1.hpp \
            $(SRC_TARGETMETRIC_SSO_DIR)/TShapeSizeOrientB1.hpp \
            $(SRC_TARGETMETRIC_SSO_DIR)/TShapeSizeOrientB2.hpp \
            $(SRC_TARGETMETRIC_SSO_DIR)/TShapeSizeOrientNB1.hpp
SRC_TARGETMETRIC_SIZE_DIR = TargetMetric/Size
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_TARGETMETRIC_SIZE_DIR)
INCLUDED_MAKEFILES += $(SRC_TARGETMETRIC_MISC_DIR)/Makefile.am

MSQ_SRCS += $(SRC_TARGETMETRIC_SIZE_DIR)/AWSizeB1.cpp \
	    $(SRC_TARGETMETRIC_SIZE_DIR)/AWSizeNB1.cpp \
	    $(SRC_TARGETMETRIC_SIZE_DIR)/TSizeB1.cpp \
	    $(SRC_TARGETMETRIC_SIZE_DIR)/TSizeNB1.cpp \
	    $(SRC_TARGETMETRIC_SIZE_DIR)/TTau.cpp 

MSQ_HDRS += $(SRC_TARGETMETRIC_SIZE_DIR)/AWSizeB1.hpp \
	    $(SRC_TARGETMETRIC_SIZE_DIR)/AWSizeNB1.hpp \
	    $(SRC_TARGETMETRIC_SIZE_DIR)/TSizeB1.hpp \
	    $(SRC_TARGETMETRIC_SIZE_DIR)/TSizeNB1.hpp \
	    $(SRC_TARGETMETRIC_SIZE_DIR)/TTau.hpp 
SRC_TARGETMETRIC_UNTANGLE_DIR = TargetMetric/Untangle
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_TARGETMETRIC_UNTANGLE_DIR)
INCLUDED_MAKEFILES += $(SRC_TARGETMETRIC_UNTANGLE_DIR)/Makefile.am

MSQ_SRCS += $(SRC_TARGETMETRIC_UNTANGLE_DIR)/AWUntangleBeta.cpp \
            $(SRC_TARGETMETRIC_UNTANGLE_DIR)/TUntangle1.cpp \
            $(SRC_TARGETMETRIC_UNTANGLE_DIR)/TUntangleBeta.cpp \
	    $(SRC_TARGETMETRIC_UNTANGLE_DIR)/TUntangleMu.cpp

MSQ_HDRS += $(SRC_TARGETMETRIC_UNTANGLE_DIR)/AWUntangleBeta.hpp \
            $(SRC_TARGETMETRIC_UNTANGLE_DIR)/TUntangle1.hpp \
            $(SRC_TARGETMETRIC_UNTANGLE_DIR)/TUntangleBeta.hpp \
	    $(SRC_TARGETMETRIC_UNTANGLE_DIR)/TUntangleMu.hpp

