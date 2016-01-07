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

void usage(){
	std::cout << "usage: mpirun -np <number of processors> ./hireconst_test_parallel <mesh file> -degree <degree> -interp <0=least square, 1=interpolation> -dim <mesh dimension>" << std::endl;
}

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
		usage();
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
				#ifdef MOAB_HAVE_MPI
					if(0==rank){
						usage();
					}
				#else
					usage();
				#endif
					return 0;
				}
			}
		}
		if(!hasdim){
		#ifdef MOAB_HAVE_MPI
			if(0==rank){
				std::cout << "Dimension of input mesh should be provided, positive and less than 3" << std::endl;
			}
		#else
			std::cout << "Dimension of input mesh should be provided, positive and less than 3" << std::endl;
		#endif
			return 0;
		}
		if(degree<=0||dim>2||dim<=0){
		#ifdef MOAB_HAVE_MPI
			if(0==rank){
				std::cout << "Input degree should be positive number;\n";
				std::cout << "Input dimesion should be positive and less than 3;" << std::endl;
			}
		#else
			std::cout << "Input degree should be positive number;\n";
			std::cout << "Input dimesion should be positive and less than 3;" << std::endl;
		#endif
			return 0;
		}
	#ifdef MOAB_HAVE_MPI
		if(0==rank){
			std::cout << "Testing on " << infile << " with dimension " << dim << "\n";
			std::string opts = interp?"interpolation":"least square fitting";
			std::cout << "High order reconstruction with degree " << degree << " " << opts << std::endl;
		}
	#else
		std::cout << "Testing on " << infile << " with dimension " << dim << "\n";
		std::string opts = interp?"interpolation":"least square fitting";
		std::cout << "High order reconstruction with degree " << degree << " " << opts << std::endl;
	#endif
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
				read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=2.0.";
			}else if(dim==1){
				read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=1.0.";
			}else{
				MB_SET_ERR(MB_FAILURE,"unsupported dimension");
			}
			read_options += (char)('0'+nghlayers);
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
#ifdef MOAB_HAVE_MPI
	int nprocs,rank;
	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Comm_size(comm,&nprocs);
	MPI_Comm_rank(comm,&rank);
#endif

	ErrorCode error;
	//mesh will be loaded and communicator pc will be updated
	error = load_meshset_hirec(infile,mbimpl,meshset,pc,degree,dim); MB_CHK_ERR(error);
	//initialize
	HiReconstruction hirec(dynamic_cast<Core*>(mbimpl),pc,meshset);
	Range elems,elems_owned;
	error = mbimpl->get_entities_by_dimension(meshset,dim,elems); MB_CHK_ERR(error);
	int nelems = elems.size();

#ifdef MOAB_HAVE_MPI
	if(pc){
		error = pc->filter_pstatus(elems,PSTATUS_GHOST,PSTATUS_NOT,-1,&elems_owned); MB_CHK_ERR(error);
	}else{
		elems_owned = elems;
	}
#endif

#ifdef MOAB_HAVE_MPI
	std::cout << "Mesh has " << nelems << " elements on Processor " << rank << " in total;";
	std::cout << elems_owned.size() << " of which are locally owned elements" << std::endl;
#else
	std::cout << "Mesh has " << nelems << " elements" << std::endl;
#endif
	//reconstruction
	if(dim==2){
		error = hirec.reconstruct3D_surf_geom(degree, interp, false); MB_CHK_ERR(error);
	}else if(dim==1){
		error = hirec.reconstruct3D_curve_geom(degree, interp, false); MB_CHK_ERR(error);
	}
#ifdef MOAB_HAVE_MPI
	std::cout << "HiRec has been done on Processor " << rank << std::endl;
#else
	std::cout << "HiRec has been done " << std::endl;
#endif
	//fitting
	double mxdist=0;
	for(Range::iterator ielem=elems_owned.begin();ielem!=elems_owned.end();++ielem){
		int nvpe; const EntityHandle* conn;
		error = mbimpl->get_connectivity(*ielem,conn,nvpe); MB_CHK_ERR(error);
		double w = 1.0/(double) nvpe;
		std::vector<double> naturalcoords2fit(nvpe,w);
		double newcoords[3],linearcoords[3];
		error = hirec.hiproj_walf_in_element(*ielem,nvpe,1,&(naturalcoords2fit[0]),newcoords); 
		if(MB_FAILURE==error) continue;
		std::vector<double> coords(3*nvpe);
		error = mbimpl->get_coords(conn,nvpe,&(coords[0])); MB_CHK_ERR(error);
		compute_linear_coords(nvpe,&(coords[0]),&(naturalcoords2fit[0]),linearcoords);
		mxdist = std::max(mxdist,Solvers::vec_distance(3,newcoords,linearcoords));
	/*#ifdef MOAB_HAVE_MPI
		std::cout << "Error on element " << *ielem << " is " << Solvers::vec_distance(3,newcoords,linearcoords) << "on Processor " << rank << std::endl;
	#else
		std::cout << "Error on element " << *ielem << " is " << Solvers::vec_distance(3,newcoords,linearcoords) << std::endl;
	#endif*/
	}
#ifdef MOAB_HAVE_MPI
	std::cout << "Maximum projection lift is " << mxdist << " on Processor " << rank << std::endl;
#else
	std::cout << "Maximum projection lift is " << mxdist << std::endl;
#endif
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