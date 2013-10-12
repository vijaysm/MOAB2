/*
 * MPASpartitioner.cpp
 *  simple mpas reader that should work in parallel to partition the element space
 *  will link against zoltan
 *  Created on: Oct 8
 */


//
//
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "MBZoltan.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "CslamUtils.hpp"
#include "moab/ParallelComm.hpp"

#include "mpi.h"
#include "pnetcdf.h"

#define NCFUNC(func) ncmpi_ ## func

//! Collective I/O mode get
#define NCFUNCAG(func) ncmpi_get ## func ## _all

//! Independent I/O mode get
#define NCFUNCG(func) ncmpi_get ## func

//! Nonblocking get (request aggregation), used so far only for ucd mesh
#define NCFUNCREQG(func) ncmpi_iget ## func

#define NCDF_SIZE MPI_Offset
#define NCDF_DIFF MPI_Offset


#include <algorithm>
#include <string>
#include <assert.h>

#include <cmath>


using namespace moab;

int ierr, rank=0;

// Error routines for use with MOAB API
#define CHKERR(CODE, MSG)                                 \
  do {                                                    \
    if (MB_SUCCESS != (CODE)) {                           \
      std::string errstr;  mbi->get_last_error(errstr);   \
      std::cerr << errstr << std::endl;                   \
      std::cerr << MSG << std::endl;                      \
      MPI_Finalize();                                     \
    }                                                     \
  } while(false)

// Error routines for use with MPI API
#define MPICHKERR(CODE, MSG)                              \
  do {                                                    \
    if (0 != CODE) {                                      \
      std::cerr << MSG << std::endl;                      \
      MPI_Finalize();                                     \
    }                                                     \
  } while(false)

#define dbgprint(MSG)                                \
  do {                                              \
      if (!rank) std::cerr << MSG << std::endl;     \
  } while(false)

#define dbgprintall(MSG)                                      \
  do {                                                        \
      std::cerr << "[" << rank << "]: " << MSG << std::endl;  \
  } while(false)


#define INS_ID(stringvar, prefix, id) \
          sprintf(stringvar, prefix, id)

#define GET_DIM(ncdim, name, val)\
    {                            \
    int gdfail = NCFUNC(inq_dimid)(ncFile, name, &ncdim);          \
    if (NC_NOERR == gdfail) {                                             \
      MPI_Offset tmp_val;                                                   \
      gdfail = NCFUNC(inq_dimlen)(ncFile, ncdim, &tmp_val);                        \
      if (NC_NOERR != gdfail) {                                           \
        std::cout<<"couldn't get dimension length for" << name<< " \n"; \
        return 1;                                              \
      }                                                                 \
      else val = tmp_val;                                               \
    } else val = 0;}

#define GET_DIMB(ncdim, name, varname, id, val) \
          INS_ID(name, varname, id); \
          GET_DIM(ncdim, name, val);

#define GET_VAR(name, id, dims) \
    {                           \
    id = -1;\
    int gvfail = NCFUNC(inq_varid)( ncFile, name, &id);   \
    if (NC_NOERR == gvfail) {       \
    int ndims;\
    gvfail = NCFUNC(inq_varndims)( ncFile, id, &ndims);\
    if (NC_NOERR == gvfail) {\
    dims.resize(ndims);    \
    gvfail = NCFUNC(inq_vardimid)( ncFile, id, &dims[0]);}}}

#define GET_1D_INT_VAR(name, id, vals) \
    {GET_VAR(name, id, vals);  \
  if (-1 != id) {\
    NCDF_SIZE ntmp;\
    int ivfail = NCFUNC(inq_dimlen)(ncFile, vals[0], &ntmp);\
    vals.resize(ntmp);\
    NCDF_SIZE ntmp1 = 0;                                                           \
    ivfail = NCFUNC(get_vara_int_all)(ncFile, id, &ntmp1, &ntmp, &vals[0]);\
    if (NC_NOERR != ivfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable " <<name <<"\n";\
      return 1;}}}

#define GET_1D_DBL_VAR(name, id, vals) \
    {std::vector<int> dum_dims;        \
  GET_VAR(name, id, dum_dims);\
  if (-1 != id) {\
    std::cout << " var:" << name << " id: " << id << "\n"; \
    NCDF_SIZE ntmp;\
    int dvfail = NCFUNC(inq_dimlen)(ncFile, dum_dims[0], &ntmp);\
    vals.resize(ntmp);\
    NCDF_SIZE ntmp1 = 0;               \
    dvfail = NCFUNC(get_vara_double_all)( ncFile, id, &ntmp1, &ntmp, &vals[0]);\
    if (NC_NOERR != dvfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";\
      return 1;}}}

#define GET_1D_DBL_VAR_RANGE(name, id, startIndex, endIndex, vals) \
    {std::vector<int> dum_dims;        \
  GET_VAR(name, id, dum_dims);\
  if (-1 != id) {\
    if (rank==0) std::cout << " var:" << name << " id: " << id << "\n"; \
    NCDF_SIZE ntmp;\
    int dvfail = NCFUNC(inq_dimlen)(ncFile, dum_dims[0], &ntmp);\
    if (rank==0) std::cout << " prepare read var " << name << " start:" << startIndex << " end:"<< endIndex <<" Actual size:" << \
       ntmp << "\n" ; \
    if (startIndex > (int)ntmp || endIndex > (int) ntmp || startIndex>=endIndex) \
      { std::cout << "bad input \n"; return 1;} \
    vals.resize(endIndex-startIndex);\
    NCDF_SIZE ntmp1 = startIndex;    NCDF_SIZE count = endIndex-startIndex;           \
    dvfail = NCFUNC(get_vara_double_all)( ncFile, id, &ntmp1, &count, &vals[0]);\
    if (NC_NOERR != dvfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";\
      return 1;}}}

#define FILL_1D_DBL_VAR(name, id, vals) \
    {std::vector<int> dum_dims;        \
  GET_VAR(name, id, dum_dims);\
  if (-1 != id) {\
    MPI_Offset ntmp;\
    int dvfail = NCFUNC(inq_dimlen)(ncFile, dum_dims[0], &ntmp);\
    std::cout <<" var:" << name << " ntmp:" << ntmp << "\n"; \
    MPI_Offset ntmp1 = 0;                                                           \
    dvfail = NCFUNC(get_vara_double_all)( ncFile, id, &ntmp1, &ntmp, vals);\
    if (NC_NOERR != dvfail) {\
      std::cout<<"FILL_1D_DBL_VAR:: Problem getting variable "<< name<<"\n";\
      return 1;}}}

/*
int get_2d_flt_var(int ncFile, const char * name, int index, std::vector<float> & vals)
{
  int id;
  std::vector<int> dum_dims;
  GET_VAR(name, id, dum_dims);
  if (-1 != id) {
    size_t ntmp;
  int dvfail = nc_inq_dimlen(ncFile, dum_dims[1], &ntmp);
  vals.resize(ntmp);
  size_t ntmp1[2] = {index, 0};

  int ntimes;
  dvfail = nc_inq_dimlen(ncFile, dum_dims[0], &ntimes);
  if (index>=ntimes) {
    std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";
    return 1;
  }
  size_t count[2] ={1, ntmp};
  dvfail = nc_get_vara_float(ncFile, id, ntmp1, count, &vals[0]);
  if (NC_NOERR != dvfail) {
    std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";
    return 1;
  }
  return 0;
}
*/
/* get the variable along an index */
#define GET_2D_FLT_VAR(name, id, index, vals) \
    {std::vector<int> dum_dims;        \
  GET_VAR(name, id, dum_dims);\
  if (-1 != id) {\
    MPI_Offset ntmp, ntmp2;\
    int dvfail = NCFUNC(inq_dimlen)( ncFile, dum_dims[1], &ntmp);\
    dvfail = NCFUNC(inq_dimlen)(ncFile, dum_dims[0], &ntmp2);\
    vals.resize(ntmp);\
    MPI_Offset ntmp1[2] = {index, 0}; \
    if (index>=ntmp2) { \
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n"; \
      return 1; \
    } \
    MPI_Offset count[2] ={1, ntmp}; \
    dvfail = NCFUNC(get_vara_float)( ncFile, id, ntmp1, count, &vals[0]);\
    if (NC_NOERR != dvfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";\
      return 1;}}}

int main(int argc, char ** argv)
{

  MPI_Comm comm = MPI_COMM_WORLD;
  ierr = MPI_Init(&argc, &argv);
  MPICHKERR(ierr, "MPI Init  failed");
  ierr = MPI_Comm_rank(comm, &rank);
  MPICHKERR(ierr, "rank failed");


  const char *data_file = "/home/iulian/source/myMoabRepo/MeshFiles/unittest/io/mpasx1.642.t.2.nc";
  if (argc>1)
  {
    data_file = argv[1];
  }

  if (0==rank)
    std::cout<<" processing file: " << data_file << "\n";
  Interface * mb = new Core;
  ErrorCode rval;

  ParallelComm* pcomm = new ParallelComm(mb, MPI_COMM_WORLD);
  int nprocs = pcomm->proc_config().proc_size();

  MBZoltan * mbZTool = new MBZoltan(mb, false, argc, argv);

  int ncFile, temp_dim;        // netcdf/exodus file

  int fail = NCFUNC(open)(comm, data_file, 0, MPI_INFO_NULL, &ncFile);
  if (NC_NOWRITE != fail) {
    std::cout<<"ReadNCDF:: problem opening Netcdf/Exodus II "<<data_file <<"\n";
    return 1;
  }
  int maxEdges;
  GET_DIM(temp_dim, "maxEdges", maxEdges);

  int nVertices;
  GET_DIM(temp_dim, "nVertices", nVertices);

  int nCells;
  GET_DIM(temp_dim, "nCells", nCells);

  if (rank==0)
  {
    std::cout << "maxEdges:" << maxEdges << "\n";
    std::cout << "nVertices:" << nVertices << "\n";
    std::cout << "nCells:" << nCells << "\n";
  }
  // see if the mesh from metadata makes sense
  // create polygons with the connectivity from conn array

  // find out how many element centers do we have to read per processor
  int numLocalCells = nCells/nprocs;

  int leftOver= nCells%nprocs;
  if (rank<leftOver)
    numLocalCells++;
  // start index
  int startCellIndex = rank*numLocalCells;
  if (rank>=leftOver)
    startCellIndex+=leftOver;
  int endCellIndex = startCellIndex + numLocalCells; // exclusive
  // so now, every processor will read xCell, yCell, zCell for cells with indices:
  // [startCellIndex, endCellIndex)
  // their global ids will be [startCellIndex+1, endCellIndex+1) // global id starts at 1 :)
  // read 3 coordinates for x, y, z, separately

  std::vector<double> x, y, z;
  GET_1D_DBL_VAR_RANGE("xCell", temp_dim, startCellIndex, endCellIndex, x) ;
  GET_1D_DBL_VAR_RANGE("yCell", temp_dim, startCellIndex, endCellIndex, y) ;
  GET_1D_DBL_VAR_RANGE("zCell", temp_dim, startCellIndex, endCellIndex, z) ;

  Range localGIDs;
  rval = mbZTool->repartition(x, y, z, startCellIndex+1, "RCB", localGIDs );
  if (rval !=MB_SUCCESS)
    std::cout <<" error in partitioning\n";

  // gather all psizes on root for printing
  int lsize[2]={ (int)localGIDs.psize(), (int)localGIDs.size()};
  /*ierr = MPI_Gather(
    void *sendbuf,
    int sendcnt,
    MPI_Datatype sendtype,
    void *recvbuf,
    int recvcnt,
    MPI_Datatype recvtype,
    int root,
    MPI_Comm comm
  );*/
  int * sizes=0;
  if (rank==0)
    sizes = new int [2* nprocs];
  ierr = MPI_Gather( lsize, 2, MPI_INT, sizes, 2, MPI_INT, 0, comm);

  if (rank==0)
  {
    for (int i=0; i<nprocs; i++)
      std::cout << " on rank:" << i <<"  psize:" << sizes[2*i] << " size: " << sizes[2*i+1] << "\n";
    delete [] sizes;
    std::cout<<" localGIDs on rank 0\n";
    localGIDs.print(0);
  }
  // form tuples (toProc, GlobalId)
  MPI_Finalize();

  return 0;
}
