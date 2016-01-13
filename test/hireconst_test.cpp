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

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;

#ifdef MOAB_HAVE_MPI
std::string read_options;
#endif

ErrorCode load_meshset_hirec(const char* infile, Interface* mbimpl, EntityHandle& meshset, ParallelComm*& pc, const int degree=0, const int dim=2);
ErrorCode test_mesh(const char* infile,const int degree, const bool interp, const int dim);

void compute_linear_coords(const int nvpe, double* elemcoords, double* naturals, double* linearcoords);

ErrorCode create_unitsq_tris(Interface *mbImpl, size_t n, std::vector<EntityHandle>& tris);
ErrorCode create_unitsq_quads(Interface *mbImpl, size_t n, std::vector<EntityHandle>& quads);
ErrorCode test_unitsq_tris();
ErrorCode test_unitsq_quads();
ErrorCode test_unitsphere();
ErrorCode test_unitcircle();

ErrorCode exact_error_torus(const double R, const double r, const double center[3],  int npnts, double *pnt,  double &error_l1, double &error_l2, double &error_li );

int main(int argc, char *argv[]){
#ifdef MOAB_HAVE_MPI
	MPI_Init(&argc,&argv);
	int nprocs,rank;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	char *infile;
	int degree=2, dim=0;
	bool interp = false;
	ErrorCode error;
	if(argc==1){
		error = test_unitcircle(); MB_CHK_ERR(error);
		error = test_unitsq_tris(); MB_CHK_ERR(error);
		error = test_unitsq_quads(); MB_CHK_ERR(error);
		error = test_unitsphere(); MB_CHK_ERR(error);
		return 0;
	}else{
		infile = argv[1]; bool hasdim=false;
		for(int i=2;i<argc;++i){
			if(i+1!=argc){
				if(std::string(argv[i])=="-degree"){
					degree = atoi(argv[++i]);
				}else if(std::string(argv[i])=="-interp"){
					interp = atoi(argv[++i]);
				}else if(std::string(argv[i])=="-dim"){
					dim = atoi(argv[++i]); hasdim = true;
				}else{
					std::cout << argv[i] << std::endl;
					std::cout << "usage: " << argv[0] << " <mesh file> -degree <degree> -interp <0=least square, 1=interpolation> -dim <mesh dimension>" << std::endl;
					return 0;
				}
			}
		}
		if(!hasdim){
			std::cout << "Dimension of input mesh should be provided, positive and less than 3" << std::endl;
			return 0;
		}
		if(degree<=0||dim>2||dim<=0){
			std::cout << "Input degree should be positive number;\n";
			std::cout << "Input dimesion should be positive and less than 3;" << std::endl;
			return 0;
		}
		//std::cout << "Testing on " << infile << " with dimension " << dim << "\n";
		std::string opts = interp?"interpolation":"least square fitting";
		//std::cout << "High order reconstruction with degree " << degree << " " << opts << std::endl;
	}
	
	error = test_mesh(infile,degree,interp,dim); MB_CHK_ERR(error);

#ifdef MOAB_HAVE_MPI
	MPI_Finalize();
#endif
}

ErrorCode load_meshset_hirec(const char* infile, Interface* mbimpl, EntityHandle& meshset, ParallelComm*& pc, const int degree, const int dim){
	ErrorCode error;
	error = mbimpl->create_meshset(moab::MESHSET_SET,meshset); MB_CHK_ERR(error);
#ifdef MOAB_HAVE_MPI
	int nprocs,rank;
	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Comm_size(comm,&nprocs);
	MPI_Comm_rank(comm,&rank);
	EntityHandle partnset;
	error = mbimpl->create_meshset(moab::MESHSET_SET,partnset); MB_CHK_ERR(error);
	if(nprocs>1)
		pc = moab::ParallelComm::get_pcomm(mbimpl,partnset,&comm);

	if(nprocs>1){
		int nghlayers = degree>0?HiReconstruction::estimate_num_ghost_layers(degree,true):0;
		if(nghlayers){
			//get ghost layers
			if(dim==2){
				read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOST=2.0.";
			}else if(dim==1){
				read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOST=1.0.";
			}else{
				MB_SET_ERR(MB_FAILURE,"unsupported dimension");
			}
			read_options += ('0'+nghlayers);
		}else{
			read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";
		}
		error = mbimpl->load_file(infile,&meshset,read_options.c_str()); MB_CHK_ERR(error);
	}else{
		error = mbimpl->load_file(infile,&meshset); MB_CHK_ERR(error);
	}
#else
	error = mbimpl->load_file(infile,&meshset); MB_CHK_ERR(error);
#endif
	return error;
}

ErrorCode test_mesh(const char* infile,const int degree, const bool interp, const int dim){
	Core moab;
	Interface* mbimpl=&moab;
	ParallelComm *pc=NULL;
	EntityHandle meshset;

	ErrorCode error;
	error = load_meshset_hirec(infile,mbimpl,meshset,pc,degree,dim); MB_CHK_ERR(error);
	//initialize
	HiReconstruction hirec(dynamic_cast<Core*>(mbimpl),pc,meshset);
	Range elems;
	error = mbimpl->get_entities_by_dimension(meshset,dim,elems); MB_CHK_ERR(error);
	int nelems = elems.size();
	//std::cout << "Mesh has " << nelems << " elements" << std::endl;
	//reconstruction
	if(dim==2){
		//error = hirec.reconstruct3D_surf_geom(degree, interp, false); MB_CHK_ERR(error);
		error = hirec.reconstruct3D_surf_geom(degree, interp, true); MB_CHK_ERR(error);
	}else if(dim==1){
		error = hirec.reconstruct3D_curve_geom(degree, interp, false); MB_CHK_ERR(error);
	}	


	//fitting
	double mxdist=0, errl1=0, errl2=0, errli=0;

	double *pnts_proj = new double[elems.size()*3];

	for(Range::iterator ielem=elems.begin();ielem!=elems.end();++ielem){
		int nvpe; const EntityHandle* conn;
		error = mbimpl->get_connectivity(*ielem,conn,nvpe); MB_CHK_ERR(error);

		double w = 1.0/(double) nvpe;
		std::vector<double> naturalcoords2fit(nvpe,w);

		double newcoords[3],linearcoords[3];
		error = hirec.hiproj_walf_in_element(*ielem,nvpe,1,&(naturalcoords2fit[0]),newcoords); MB_CHK_ERR(error);

		pnts_proj[3*(*ielem-*elems.begin())] = newcoords[0];
		pnts_proj[3*(*ielem-*elems.begin())+1] = newcoords[1];
		pnts_proj[3*(*ielem-*elems.begin())+2] = newcoords[2];

		std::vector<double> coords(3*nvpe);
		error = mbimpl->get_coords(conn,nvpe,&(coords[0])); MB_CHK_ERR(error);
		compute_linear_coords(nvpe,&(coords[0]),&(naturalcoords2fit[0]),linearcoords);
		mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));	
	}

	double center[3]={0,0,0};
	double R=1, r=0.3;

	error = exact_error_torus(R, r, center, (int)elems.size(), pnts_proj, errl1, errl2, errli);MB_CHK_ERR(error);

	//std::cout << "Maximum projection lift is " << mxdist << std::endl;
	std::cout<<"Errors using exact torus: L1 = "<<errl1<<", L2 = "<<errl2<<", Linf = "<<errli<<std::endl;

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
	double istep = 1.0/(double) (n-1);
	for(size_t j=0;j<n;++j){
		for(size_t i=0;i<n;++i){
			double coord[3] = {i*istep,j*istep,0};
			EntityHandle temp;
			error = mbImpl->create_vertex(coord,temp); MB_CHK_ERR(error);
			verts[j*n+i] = temp;
		}
	}

	for(size_t jj=0;jj<n-1;++jj){
		for(size_t ii=0;ii<n-1;++ii){
			EntityHandle conn[3]={verts[jj*n+ii],verts[(jj+1)*n+ii+1],verts[(jj+1)*n+ii]};
			error = mbImpl->create_element(MBTRI,conn,3,tris[istr+2*jj*(n-1)+2*ii]); MB_CHK_ERR(error);
			conn[0] = verts[jj*n+ii]; conn[1] = verts[jj*n+ii+1]; conn[2] = verts[(jj+1)*n+ii+1];
			error = mbImpl->create_element(MBTRI,conn,3,tris[istr+2*jj*(n-1)+2*ii+1]); MB_CHK_ERR(error);
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
	double istep = 1.0/(double) (n-1);
	for(size_t j=0;j<n;++j){
		for(size_t i=0;i<n;++i){
			double coord[3] = {i*istep,j*istep,0};
			error = mbImpl->create_vertex(coord,verts[j*n+i]); MB_CHK_ERR(error);
		}
	}

	for(size_t jj=0;jj<n-1;++jj){
		for(size_t ii=0;ii<n-1;++ii){
			EntityHandle conn[4]={verts[jj*n+ii],verts[jj*n+ii+1],verts[(jj+1)*n+ii+1],verts[(jj+1)*n+ii]};
			error = mbImpl->create_element(MBQUAD,conn,4,quads[istr+jj*(n-1)+ii]); MB_CHK_ERR(error);
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
			double mxdist=0;
			for(size_t itri=0;itri<tris.size();++itri){
				const int nvpe = 3;
				double naturalcoords2fit[nvpe] = {1.0/(double) nvpe,1.0/(double) nvpe,1.0/(double) nvpe}, newcoords[3];
				error = hirec.hiproj_walf_in_element(tris[itri],nvpe,1,naturalcoords2fit,newcoords); MB_CHK_ERR(error);
				std::vector<EntityHandle> conn;
				error = mbImpl->get_connectivity(&(tris[itri]),1,conn); MB_CHK_ERR(error);
				double coords[3*nvpe],linearcoords[3];
				error = mbImpl->get_coords(&(conn[0]),nvpe,coords); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,coords,naturalcoords2fit,linearcoords);
				mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));
			}
			std::cout << "triangulated unit square n= " << n << " degree= " << degree << " interpolation:\n";
			std::cout << "maximum projection lift is " << mxdist << std::endl;
			mxdist = 0;
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
				mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));
			}
			std::cout << "unit square n= " << n << " degree= " << degree << " least square:\n";
			std::cout << "maximum projection lift is " << mxdist << std::endl;
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
			double mxdist=0;
			for(size_t iquad=0;iquad<quads.size();++iquad){
				const int nvpe = 4; double w=1.0/(double) nvpe;
				double naturalcoords2fit[nvpe] = {w,w,w,w}, newcoords[3];
				error = hirec.hiproj_walf_in_element(quads[iquad],nvpe,1,naturalcoords2fit,newcoords); MB_CHK_ERR(error);
				std::vector<EntityHandle> conn;
				error = mbImpl->get_connectivity(&(quads[iquad]),1,conn); MB_CHK_ERR(error);
				double coords[3*nvpe],linearcoords[3];
				error = mbImpl->get_coords(&(conn[0]),nvpe,coords); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,coords,naturalcoords2fit,linearcoords);
				mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));
			}
			std::cout << "quadrilateral unit square n= " << n << " degree= " << degree << " interpolation:\n";
			std::cout << "maximum projection lift is " << mxdist << std::endl;
			mxdist = 0;
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
				mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));
			}
			std::cout << "quadrilateral unit square n= " << n << " degree= " << degree << " least square:\n";
			std::cout << "maximum projection lift is " << mxdist << std::endl;
		}
	}
	return error;
}

ErrorCode test_unitsphere(){
	//path to test files
#ifdef MESHDIR
	int nfiles = 4;
	char *filenames[] = {STRINGIFY(MESHDIR) "/sphere_tris_5.vtk", STRINGIFY(MESHDIR) "/sphere_tris_20.vtk", STRINGIFY(MESHDIR) "/sphere_quads_5.vtk", STRINGIFY(MESHDIR) "/sphere_quads_20.vtk"};
#else
#error Specify MESHDIR to compile test
#endif
	ErrorCode error;
	int maxdeg = 6;
	for(int ifile=0;ifile<nfiles;++ifile){
		Core moab;
		Interface* mbimpl=&moab;
		ParallelComm *pc=NULL;
		EntityHandle meshset;
		//load file
		error = load_meshset_hirec(filenames[ifile],mbimpl,meshset,pc,maxdeg); MB_CHK_ERR(error);
		//initialize
		HiReconstruction hirec(&moab,pc,meshset);
		Range elems;
		error = mbimpl->get_entities_by_dimension(meshset,2,elems); MB_CHK_ERR(error);
		//reconstruction
		for(int degree=1;degree<=maxdeg;++degree){
			hirec.reconstruct3D_surf_geom(degree, true, false, true);
			//fitting
			double mxdist=0,mxerr=0;
			for(Range::iterator ielem=elems.begin();ielem!=elems.end();++ielem){
				int nvpe; const EntityHandle* conn;
				error = mbimpl->get_connectivity(*ielem,conn,nvpe); MB_CHK_ERR(error);
				double w = 1.0/(double) nvpe;
				std::vector<double> naturalcoords2fit(nvpe,w);
				double newcoords[3],linearcoords[3];
				error = hirec.hiproj_walf_in_element(*ielem,nvpe,1,&(naturalcoords2fit[0]),newcoords); MB_CHK_ERR(error);
				std::vector<double> coords(3*nvpe);
				error = mbimpl->get_coords(conn,nvpe,&(coords[0])); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,&(coords[0]),&(naturalcoords2fit[0]),linearcoords);
				mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));
				mxerr = std::max(mxerr,fabs(Solvers::vec_2norm(3,newcoords)-1));
			}
			std::cout << filenames[ifile] << ": unit sphere" << " degree= " << degree << " interpolation:\n";
			std::cout << "maximum projection lift is " << mxdist << ", maximum error is " << mxerr << std::endl;
			mxdist = 0; mxerr = 0;
			hirec.reconstruct3D_surf_geom(degree, false, false, true);
			//fitting
			for(Range::iterator ielem=elems.begin();ielem!=elems.end();++ielem){
				int nvpe; const EntityHandle* conn;
				error = mbimpl->get_connectivity(*ielem,conn,nvpe); MB_CHK_ERR(error);
				double w = 1.0/(double) nvpe;
				std::vector<double> naturalcoords2fit(nvpe,w);
				double newcoords[3],linearcoords[3];
				error = hirec.hiproj_walf_in_element(*ielem,nvpe,1,&(naturalcoords2fit[0]),newcoords); MB_CHK_ERR(error);
				std::vector<double> coords(3*nvpe);
				error = mbimpl->get_coords(conn,nvpe,&(coords[0])); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,&(coords[0]),&(naturalcoords2fit[0]),linearcoords);
				mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));
				mxerr = std::max(mxerr,fabs(Solvers::vec_2norm(3,newcoords)-1));
			}
			std::cout << filenames[ifile] << ": unit sphere" << " degree= " << degree << " least square:\n";
			std::cout << "maximum projection lift is " << mxdist << ", maximum error is " << mxerr << std::endl;
		}
	}
	return error;
}

ErrorCode test_unitcircle(){
	//path to test files
#ifdef MESHDIR
	int nfiles = 4;
	char *filenames[] = {STRINGIFY(MESHDIR) "/circle_3.vtk", STRINGIFY(MESHDIR) "/circle_4.vtk", STRINGIFY(MESHDIR) "/circle_10.vtk", STRINGIFY(MESHDIR) "/circle_20.vtk"};
#else
#error Specify MESHDIR to compile test
#endif
	ErrorCode error;
	int maxdeg = 6;
	for(int ifile=0;ifile<nfiles;++ifile){
		Core moab;
		Interface* mbimpl=&moab;
		ParallelComm *pc=0;
		EntityHandle meshset;
		int dim=1;
		//load file
		error = load_meshset_hirec(filenames[ifile],mbimpl,meshset,pc,maxdeg,dim); MB_CHK_ERR(error);
		//initialize
		HiReconstruction hirec(&moab,pc,meshset);
		Range edges;
		error = mbimpl->get_entities_by_dimension(meshset,dim,edges); MB_CHK_ERR(error);
		//reconstruction
		for(int degree=1;degree<=maxdeg;++degree){
			hirec.reconstruct3D_curve_geom(degree, true, false, true);
			//fitting
			double mxdist=0,mxerr=0;
			for(Range::iterator iedge=edges.begin();iedge!=edges.end();++iedge){
				int nvpe; const EntityHandle* conn;
				error = mbimpl->get_connectivity(*iedge,conn,nvpe); MB_CHK_ERR(error);
				double w = 1.0/(double) nvpe;
				std::vector<double> naturalcoords2fit(nvpe,w);
				double newcoords[3],linearcoords[3];
				error = hirec.hiproj_walf_in_element(*iedge,nvpe,1,&(naturalcoords2fit[0]),newcoords); MB_CHK_ERR(error);
				std::vector<double> coords(3*nvpe);
				error = mbimpl->get_coords(conn,nvpe,&(coords[0])); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,&(coords[0]),&(naturalcoords2fit[0]),linearcoords);
				mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));
				mxerr = std::max(mxerr,fabs(Solvers::vec_2norm(3,newcoords)-1));
			}
			std::cout << filenames[ifile] << ": unit circle" << " degree= " << degree << " interpolation:\n";
			std::cout << "maximum projection lift is " << mxdist << ", maximum error is " << mxerr << std::endl;
			mxdist = 0; mxerr = 0;
			hirec.reconstruct3D_curve_geom(degree, false, false, true);
			//fitting
			for(Range::iterator iedge=edges.begin();iedge!=edges.end();++iedge){
				int nvpe; const EntityHandle* conn;
				error = mbimpl->get_connectivity(*iedge,conn,nvpe); MB_CHK_ERR(error);
				double w = 1.0/(double) nvpe;
				std::vector<double> naturalcoords2fit(nvpe,w);
				double newcoords[3],linearcoords[3];
				error = hirec.hiproj_walf_in_element(*iedge,nvpe,1,&(naturalcoords2fit[0]),newcoords); MB_CHK_ERR(error);
				std::vector<double> coords(3*nvpe);
				error = mbimpl->get_coords(conn,nvpe,&(coords[0])); MB_CHK_ERR(error);
				compute_linear_coords(nvpe,&(coords[0]),&(naturalcoords2fit[0]),linearcoords);
				mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));
				mxerr = std::max(mxerr,fabs(Solvers::vec_2norm(3,newcoords)-1));
			}
			std::cout << filenames[ifile] << ": unit circle" << " degree= " << degree << " least square:\n";
			std::cout << "maximum projection lift is " << mxdist << ", maximum error is " << mxerr << std::endl;
		}
	}
	return error;
}

ErrorCode exact_error_torus(const double R, const double r, const double center[], int npnts, double *pnts, double &error_l1, double &error_l2, double &error_li)
{
  error_l1 = 0;
  error_l2 = 0;
  error_li = 0;

  double x=0, y=0, z=0;
  double error_single=0;

  for (int i = 0; i< npnts; i++)
    {
      x = pnts[3*i]-center[0]; y = pnts[3*i+1]-center[1]; z = pnts[3*i+2]-center[2];

      error_single = fabs (sqrt(pow(sqrt(pow(x,2)+pow(y,2)) - R, 2)+pow(z,2)) - r);
      error_l1=error_l1+error_single;
      error_l2=error_l2+error_single*error_single;

      if (error_li<error_single)
        error_li = error_single;
    }
  error_l1 = error_l1/(double)npnts;
  error_l2 = sqrt(error_l2/(double)npnts);

  return MB_SUCCESS;
}
