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

int main(int argc, char *argv[]){
#ifdef MOAB_HAVE_MPI
	MPI_Init(&argc,&argv);
	int nprocs,rank;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	if(argc!=2){
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
		error = mbimpl->load_file(argv[1],&meshset,read_options.c_str()); MB_CHK_ERR(error);
	}else{
		error = mbimpl->load_file(argv[1],&meshset); MB_CHK_ERR(error);
	}
#else
	error = mbimpl->load_file(argv[1],&meshset); MB_CHK_ERR(error);
#endif

	HiReconstruction hirec(&moab,pc,meshset);

#ifdef MOAB_HAVE_MPI
	MPI_Finalize();
#endif
}