/*This unit test is for high order reconstruction capability, which could be used for mesh refinement*/
#include <iostream>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "../RefineMesh/moab/NestedRefine.hpp"
#include "../DiscreteGeometry/moab/Solvers.hpp"
#include "../DiscreteGeometry/moab/HiReconstruction.hpp"
#include "TestUtil.hpp"
#include <math.h>

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab_mpi.h"
#endif

using namespace moab;

#ifdef MOAB_HAVE_MPI
std::string read_options;
#endif

ErrorCode create_unitsq_tris(int n);
ErrorCode create_unitsq_quads(int n);
ErrorCode test_unitsq_tris();
ErrorCode test_unitsq_quads();

int main(int argc, char *argv[]){
#ifdef MOAB_HAVE_MPI
	MPI_Init(&argc,&argv);
	int nprocs,rank;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	char *infile;
	if(argc==1){
		
	}else if(argc==2){
		infile = argv[1];
	}else{
		std::cout << "usage: " << argv[0] << " <mesh file>\n";
		return 0;
	}
	Core moab;
	Interface* mbimpl=&moab;
	ParallelComm *pc=NULL;
	EntityHandle meshset;

	ErrorCode error;
	error = mbimpl->create_meshset(moab::MESHSET_SET,meshset); MB_CHK_ERR(error);
#ifdef MOAB_HAVE_MPI
	MPI_Comm comm = MPI_COMM_WORLD;
	EntityHandle partnset;
	error = mbimpl->create_meshset(moab::MESHSET_SET,partnset); MB_CHK_ERR(error);
	pc = moab::ParallelComm::get_pcomm(mbimpl,partnset,&comm);

	if(nprocs>1){
		read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";
		error = mbimpl->load_file(infile,&meshset,read_options.c_str()); MB_CHK_ERR(error);
	}else{
		error = mbimpl->load_file(infile,&meshset); MB_CHK_ERR(error);
	}
#else
	error = mbimpl->load_file(infile,&meshset); MB_CHK_ERR(error);
#endif

	int maxdeg = 6;
	//get ghost layers
#ifdef MOAB_HAVE_MPI
	if(nprocs>1){
		int nghostlayers = HiReconstruction::estimate_num_ghost_layers(maxdeg,true);

	}
#endif

	HiReconstruction hirec(dynamic_cast<Core*>(mbimpl),pc,meshset);

#ifdef MOAB_HAVE_MPI
	MPI_Finalize();
#endif
}

ErrorCode load_meshset_hirec(char* infile, Interface* mbimpl, EntityHandle& meshset, const int degree=0, MPI_Comm comm=MPI_COMM_WORLD, ParallelComm* pc=NULL){
	ErrorCode error;
	error = mbimpl->create_meshset(moab::MESHSET_SET,meshset); MB_CHK_ERR(error);
#ifdef MOAB_HAVE_MPI
	int nprocs,rank;
	MPI_Comm_size(comm,&nprocs);
	MPI_Comm_rank(comm,&rank);
	EntityHandle partnset;
	error = mbimpl->create_meshset(moab::MESHSET_SET,partnset); MB_CHK_ERR(error);
	pc = moab::ParallelComm::get_pcomm(mbimpl,partnset,&comm);

	if(nprocs>1){
		int nghlayers = degree>0?HiReconstruction::estimate_num_ghost_layers(degree,true):0;
		read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";
		error = mbimpl->load_file(infile,&meshset,read_options.c_str()); MB_CHK_ERR(error);
	}else{
		error = mbimpl->load_file(infile,&meshset); MB_CHK_ERR(error);
	}
#else
	error = mbimpl->load_file(infile,&meshset); MB_CHK_ERR(error);
#endif
	return error;
}

ErrorCode create_unitsq_tris(Interface *mbImpl, size_t n, std::vector<EntityHandle>& tris){
	if(n<2){
		MB_SET_ERR(MB_FAILURE,"n must be at least 2");
	}
	ErrorCode error;
	std::vector<EntityHandle> verts(n*n);
	size_t istr = tris.size();
	tris.resize(istr+2*(n-1)*(n-1));
	double istep = 1.0/(double) n;
	for(size_t i=0;i<n;++i){
		for(size_t j=0;j<n;++i){
			double coord[3] = {i*istep,j*istep,0};
			error = mbImpl->create_vertex(coord,verts[i*n+j]); MB_CHK_ERR(error);
		}
	}

	for(size_t ii=0;ii<n-1;++ii){
		for(size_t jj=0;jj<n-1;++jj){
			EntityHandle conn[3]={verts[ii*n+jj],verts[(ii+1)*n+jj+1],verts[(ii+1)*n+jj]};
			error = mbImpl->create_element(MBTRI,conn,3,tris[istr+ii*(n-1)+2*jj]); MB_CHK_ERR(error);
			conn[0] = verts[ii*n+jj]; conn[1] = verts[ii*n+jj+1]; conn[2] = verts[(ii+1)*n+jj+1];
			error = mbImpl->create_element(MBTRI,conn,3,tris[istr+ii*(n-1)+2*jj+1]); MB_CHK_ERR(error);
		}
	}
	return error;
}

ErrorCode create_unitsq_quads(Interface *mbImpl, size_t n, std::vector<EntityHandle>& quads){
	if(n<2){
		MB_SET_ERR(MB_FAILURE,"n must be at least 2");
	}
	ErrorCode error;
	std::vector<EntityHandle> verts(n*n);
	size_t istr = quads.size();
	quads.resize(istr+(n-1)*(n-1));
	double istep = 1.0/(double) n;
	for(size_t i=0;i<n;++i){
		for(size_t j=0;j<n;++i){
			double coord[3] = {i*istep,j*istep,0};
			error = mbImpl->create_vertex(coord,verts[i*n+j]); MB_CHK_ERR(error);
		}
	}

	for(size_t ii=0;ii<n-1;++ii){
		for(size_t jj=0;jj<n-1;++jj){
			EntityHandle conn[4]={verts[ii*n+jj],verts[ii*n+jj+1],verts[(ii+1)*n+jj+1],verts[(ii+1)*n+jj]};
			error = mbImpl->create_element(MBQUAD,conn,4,quads[istr+ii*(n-1)+jj]); MB_CHK_ERR(error);
		}
	}
	return error;
}

ErrorCode test_unitsq_tris(){
	ErrorCode error;
	for(size_t n=2;n<=8;++n){
		Core moab;
		Interface *mbImpl=&moab;
		std::vector<EntityHandle> tris;
		error = create_unitsq_tris(mbImpl,n,tris); MB_CHK_ERR(error);
		EntityHandle meshIn = 0;
		HiReconstruction hirec(dynamic_cast<Core*>(mbImpl),0,meshIn);
		for(int degree=1;degree<=6;++degree){
			//reconstruct geometry, interpolation
			hirec.reconstruct3D_surf_geom(degree, true, false, true);
			//test fitting result
			for(size_t itri=0;itri<tris.size();++itri){
				const int nvpe = 3;
				double naturalcoords2fit[nvpe] = {1.0/(double) nvpe,1.0/(double) nvpe,1.0/(double) nvpe}, newcoords[3];
				error = hirec.hiproj_walf_in_element(tris[itri],nvpe,1,naturalcoords2fit,newcoords); MB_CHK_ERR(error);
				std::vector<EntityHandle> conn;
				error = mbImpl->get_connectivity(&(tris[itri]),1,conn); MB_CHK_ERR(error);
				double coords[3*nvpe],linearcoords[3];
				error = mbImpl->get_coords(&(conn[0]),nvpe,coords); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,coords,naturalcoords2fit,linearcoords);
				assert(fabs(newcoords[0]-linearcoords[0])<1e-2);
				assert(fabs(newcoords[1]-linearcoords[1])<1e-2);
				assert(fabs(newcoords[2]-linearcoords[2])<1e-2);
			}

			//reconstruct geometry, least square fitting
			hirec.reconstruct3D_surf_geom(degree, false, false, true);
			//test fitting result
			for(size_t itri=0;itri<tris.size();++itri){
				const int nvpe = 3;
				double naturalcoords2fit[nvpe] = {1.0/(double) nvpe,1.0/(double) nvpe,1.0/(double) nvpe}, newcoords[3];
				error = hirec.hiproj_walf_in_element(tris[itri],nvpe,1,naturalcoords2fit,newcoords); MB_CHK_ERR(error);
				std::vector<EntityHandle> conn;
				error = mbImpl->get_connectivity(&(tris[itri]),1,conn); MB_CHK_ERR(error);
				double coords[3*nvpe],linearcoords[3];
				error = mbImpl->get_coords(&(conn[0]),nvpe,coords); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,coords,naturalcoords2fit,linearcoords);
				assert(fabs(newcoords[0]-linearcoords[0])<1e-2);
				assert(fabs(newcoords[1]-linearcoords[1])<1e-2);
				assert(fabs(newcoords[2]-linearcoords[2])<1e-2);
			}
		}
	}
	return error;
}

void compute_linear_coords(const int nvpe, double* elemcoords, double* naturals, double* linearcoords){
	assert(elemcoords&&linearcoords);
	for(int i=0;i<3;++i){
		linearcoords[i] = 0;
		for(int j=0;j<nvpe;++j){
			linearcoords[i] += naturals[j]*elemcoords[3*j+i];
		}
	}
}

ErrorCode test_unitsq_quads(){
	ErrorCode error;
	for(size_t n=2;n<=8;++n){
		Core moab;
		Interface *mbImpl=&moab;
		std::vector<EntityHandle> quads;
		error = create_unitsq_quads(mbImpl,n,quads); MB_CHK_ERR(error);
		EntityHandle meshIn = 0;
		HiReconstruction hirec(dynamic_cast<Core*>(mbImpl),0,meshIn);
		for(int degree=1;degree<=6;++degree){
			//reconstruct geometry, interpolation
			hirec.reconstruct3D_surf_geom(degree, true, false, true);
			//test fitting result
			for(size_t iquad=0;iquad<quads.size();++iquad){
				const int nvpe = 4; double w=1.0/(double) nvpe;
				double naturalcoords2fit[nvpe] = {w,w,w,w}, newcoords[3];
				error = hirec.hiproj_walf_in_element(quads[iquad],nvpe,1,naturalcoords2fit,newcoords); MB_CHK_ERR(error);
				std::vector<EntityHandle> conn;
				error = mbImpl->get_connectivity(&(quads[iquad]),1,conn); MB_CHK_ERR(error);
				double coords[3*nvpe],linearcoords[3];
				error = mbImpl->get_coords(&(conn[0]),nvpe,coords); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,coords,naturalcoords2fit,linearcoords);
				assert(fabs(newcoords[0]-linearcoords[0])<1e-2);
				assert(fabs(newcoords[1]-linearcoords[1])<1e-2);
				assert(fabs(newcoords[2]-linearcoords[2])<1e-2);
			}

			//reconstruct geometry, least square fitting
			hirec.reconstruct3D_surf_geom(degree, false, false, true);
			//test fitting result
			for(size_t iquad=0;iquad<quads.size();++iquad){
				const int nvpe = 4; double w=1.0/(double) nvpe;
				double naturalcoords2fit[nvpe] = {w,w,w,w}, newcoords[3];
				error = hirec.hiproj_walf_in_element(quads[iquad],nvpe,1,naturalcoords2fit,newcoords); MB_CHK_ERR(error);
				std::vector<EntityHandle> conn;
				error = mbImpl->get_connectivity(&(quads[iquad]),1,conn); MB_CHK_ERR(error);
				double coords[3*nvpe],linearcoords[3];
				error = mbImpl->get_coords(&(conn[0]),nvpe,coords); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,coords,naturalcoords2fit,linearcoords);
				assert(fabs(newcoords[0]-linearcoords[0])<1e-2);
				assert(fabs(newcoords[1]-linearcoords[1])<1e-2);
				assert(fabs(newcoords[2]-linearcoords[2])<1e-2);
			}
		}
	}
	return error;
}

ErrorCode test_unitsphere_tris(){
	//path to test files
#ifdef MESHDIR
	int nfiles = 2;
	char *filenames[] = {STRINGIFY(MESHDIR) "/sphere_tris_5.vtk", STRINGIFY(MESHDIR) "/sphere_tris_20.vtk"};
#else
#error Specify MESHDIR to compile test
#endif
	for(int ifile=0;ifile<nfiles;++ifile){
		//load file
		char *file = filenames[ifile];

		//initialize

		//reconstruction

		//fitting
	}
}

ErrorCode test_unitsphere_quads(){
	//path to test files
}

ErrorCode test_unitcircle(){

}

//parallel test