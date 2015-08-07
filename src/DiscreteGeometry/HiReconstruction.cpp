#include "moab/HiReconstruction.hpp"
#include "moab/Solvers.hpp"

#include <deque>
#include <iostream>

namespace moab
{
	HiReconstruction::HiReconstruction(Core *impl, ParallelComm *comm, EntityHandle meshIn, int minpnts, bool recwhole)
	 : mbImpl(impl),pcomm(comm),_mesh2rec(meshIn),_MINPNTS(minpnts){
	 	assert(NULL!=impl);
	 	ErrorCode error;

	 #ifdef MOAB_HAVE_MPI
	 	if(!pcomm){
	 		pcomm = moab::ParallelComm::get_pcomm(mbImpl,0);
	 	}
	 #endif

	 	error = initialize(recwhole);
	 	if(MB_SUCCESS!=error){
	 		std::cout << "Error initializing HiReconstruction\n" << std::endl;
	 		exit(1);
	 	}
	 }

	HiReconstruction::~HiReconstruction(){
	#ifdef MOAB_HAVE_AHF
	 	ahf = NULL;
	#else
	 	delete ahf;
	#endif
	}

	ErrorCode HiReconstruction::initialize(bool recwhole){
		ErrorCode error;

		ahf = new HalfFacetRep(mbImpl,pcomm,_mesh2rec);
		if(!ahf){
			return MB_MEMORY_ALLOCATION_FAILED;
		}

		error = ahf->initialize(); MB_CHK_ERR(error);
		if(HalfFacetRep::CURVE==ahf->thismeshtype){
			_dim = 1; _MAXPNTS = 13;
		}else if(HalfFacetRep::SURFACE==ahf->thismeshtype){
			_dim = 2; _MAXPNTS = 128;
		}else{
			MB_SET_ERR(MB_FAILURE,"Encounted a non-manifold mesh or a mesh with volume elements")
		}

		error = ahf->get_entity_ranges(_inverts,_inedges,_infaces,_incells); MB_CHK_ERR(error);

		//get locally hosted vertices by filtering pstatus
	#ifdef MOAB_HAVE_MPI
		error = pcomm->filter_pstatus(_inverts,PSTATUS_GHOST,PSTATUS_NOT,-1,&_verts2rec); MB_CHK_ERR(error);
	#else
		_verts2rec = _inverts;
	#endif
		_nv2rec = _verts2rec.size();

		if(recwhole){
			//compute normals(surface) or tangent vector(curve) for all locally hosted vertices
			if(2==_dim){
				compute_average_vertex_normals_surf();
			}else{
				compute_average_vertex_tangents_curve();
			}
			_hasderiv = true;
		}
		return error;
	}

	/***************************************************
	 *  User Interface for Reconstruction of Geometry  *
	 ***************************************************/

	ErrorCode reconstruct3D_surf_geom(int degree, bool interp, bool safeguard, bool reset){
		if(_hasfittings&&!reset){
			//This object has precomputed fitting results and user don't want to reset
			return MB_SUCCESS;
		}else{
			_initfittings = _hasfittings = false;
		}
	}

	ErrorCode reconstruct3D_surf_geom(size_t npts, int* degrees, bool* interps, bool safeguard, bool reset){

	}

	ErrorCode reconstruct3D_curve_geom(int degree, bool interp, bool safeguard, bool reset){

	}

	ErrorCode reconstruct3D_curve_geom(size_t npts, int* degrees, bool* interps, bool safeguard, bool reset){

	}

	ErrorCode polyfit3d_walf_surf_vertex(const EntityHandle vid, const bool interp, int degree, int minpnts, const bool safeguard, double* coords, int* degree_out, double* coeffs){

	}

	ErrorCode polyfit3d_walf_curve_vertex(const EntityHandle vid, const bool interp, int degree, int minpnts, const bool safeguard, double* coords, int* degree_out, double* coeffs){

	}

	/**************************************************************
	 *  User Interface for Evaluation via Reconstructed Geometry  *
	 **************************************************************/

	ErrorCode hiproj_walf_in_element(EntityHandle elem, const int nvpe, const int npts2fit, const double* naturalcoords2fit, double* newcoords){

	}

	ErrorCode hiproj_walf_around_vertex(EntityHandle vid, const int npts2fit, const double* coords2fit, double* hiproj_new){

	}

	void walf3d_surf_vertex_eval(const double* local_origin, const double* local_coords, const int local_deg, const double* local_coeffs, const bool interp, const int npts2fit, const double* coords2fit, double* hiproj_new){

	}

	void walf3d_curve_vertex_eval(const double* local_origin, const double* local_coords, const int local_deg, const double* local_coeffs, const bool interp, const int npts2fit, const double* coords2fit, double* hiproj_new){

	}

	/****************************************************************
	 *  Basic Internal Routines to initialize and set fitting data  *
	 ****************************************************************/

	 int estimate_num_rings(int degree, bool interp){

	 }

	 ErrorCode obtain_nring_ngbvs(const EntityHandle vid, int ring, const int minpnts, Range& ngbvs){

	 }

	 void initialize_surf_geom(const int degree){

	 }

	 void initialize_surf_geom(const size_t npts, const int* degrees){

	 }

	 void initialize_3Dcurve_geom(const int degree){

	 }

	 void initialize_3Dcurve_geom(const size_t npts, const int* degrees){
	 	
	 }
}//namespace moab