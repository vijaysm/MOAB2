/** @example GenLargeMesh.cpp \n
 * \brief Create a large structured mesh, partitioned \n
 * <b>To run</b>: mpiexec -np 2 GenLargeMesh \n
 *
 *  It shows how to load create a mesh on the fly, on multiple
 *  processors, block wise
 *  Each processor will create its version of a block mesh, partitioned
 *  as AxBxC blocks. Each block will be with blockSize^3 hexahedrons, and will get a
 *  different PARALLEL_PARTITION tag
 *  When -t option is used, instead of each hex, we are creating 6 tetrahedrons
 *
 *  The number of tasks will be MxNxK, and it must match the mpi size
 *  Each task will generate its mesh at location (m,n,k)
 *
 *  By default M=1, N=1, K=1, so by default it should be launched on 1 proc
 *  By default, blockSize is 4, and A=2, B=2, C=2, so each task will generate locally
 *    blockSize^3 x A x B x C hexahedrons (value = 64x8 = 512 hexas, in 8 partitions)
 *   (if -t, multiple by 6 for total number of cells/tets)
 *  The total number of partitions will be A*B*C*M*N*K (default 8)
 *
 *  Each part in partition will get a proper tag
 *
 *  The vertices will get a proper global id, which will be used to resolve the
 *  shared entities
 *  The output will be written in parallel, and we will try sizes as big as we can
 *  (up to a billion vertices, because we use int for global ids)
 *
 *  Within each partition, the hexas will be numbered contiguously, and also the
 *  vertices; The global id will be determined easily, for a vertex, but the entity
 *  handle space will be more interesting to handle, within a partition (we want
 *  contiguous handles within a partition)
 *
 *  We may or may not use ScdInterface, because we want control over global id and
 *  entity handle within a partition
 *
 *
 *  to run: ./GenLargeMesh
 *
 *  when launched on more procs, you have to make sure
 *  num procs = M*N*K
 *
 *  so you can launch with
 *  mpiexec -np 8 ./GenLargeMesh -M 2 -N 2 -K 2
 *
 *  We also added -q option; it works now only for hexa mesh, it will generate
 *  quadratic hex27 elements
 */

#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/CN.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/MergeMesh.hpp"
#include <time.h>

#include <iostream>
#include <vector>

using namespace moab;
using namespace std;

#define CHECKE(message) if(rval != MB_SUCCESS) { cout << (message) << "\n"; MPI_Finalize(); return 1; }
int main(int argc, char **argv)
{
  int A=2, B=2, C=2, M=1, N=1, K=1;
  int blockSize = 4;

  bool newMergeMethod=false;
  bool quadratic=false;
  bool keep_skins=false;
  bool tetra = false;
  bool adjEnts = false;

  MPI_Init(&argc, &argv);

  ProgOptions opts;

  opts.addOpt<int>(std::string("blockSize,b"),
      std::string("Block size of mesh (default=4)"), &blockSize);
  opts.addOpt<int>(std::string("xproc,M"),
      std::string("Number of processors in x dir (default=1)"), &M);
  opts.addOpt<int>(std::string("yproc,N"),
      std::string("Number of processors in y dir (default=1)"), &N);
  opts.addOpt<int>(std::string("zproc,K"),
      std::string("Number of processors in z dir (default=1)"), &K);

  opts.addOpt<int>(std::string("xblocks,A"),
      std::string("Number of blocks on a task in x dir (default=2)"), &A);
  opts.addOpt<int>(std::string("yblocks,B"),
      std::string("Number of blocks on a task in y dir (default=2)"), &B);
  opts.addOpt<int>(std::string("zblocks,C"),
      std::string("Number of blocks on a task in x dir (default=2)"), &C);

  opts.addOpt<void>("newMerge,w", "use new merging method",
          &newMergeMethod);

  opts.addOpt<void>("quadratic,q", "use hex 27 elements",
      &quadratic);

  opts.addOpt<void>("keep_skins,k", "keep skins with shared entities",
        &keep_skins);

  opts.addOpt<void>("tetrahedrons,t", "generate tetrahedrons ",
        &tetra);

  opts.addOpt<void>("faces_edges,f", "create all faces and edges", &adjEnts);

  vector<string> intTagNames;
  string firstIntTag;
  opts.addOpt<string>(std::string("int_tag_vert,i"), string("add integer tag on vertices"), &firstIntTag);

  vector<string> doubleTagNames;
  string firstDoubleTag;
  opts.addOpt<string>(std::string("double_tag_cell,d"), string("add double tag on cells"), &firstDoubleTag);

  string outFileName="test1.h5m";
  opts.addOpt<std::string> ("outFile,o", "Specify the output file name string (default test1.h5m)", &outFileName );

  opts.parseCommandLine(argc, argv);

  opts.getOptAllArgs("int_tag_vert,i", intTagNames);
  opts.getOptAllArgs("double_tag_cell,d", doubleTagNames);

  Interface* mb = new Core;
  if (NULL == mb)
  {
    MPI_Finalize();
    return 1;
  }
  ReadUtilIface *iface;
  ErrorCode rval = mb->query_interface(iface);
  CHECKE("Can't get reader interface");

  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  if (M*N*K != size)
  {
    if (rank==0)
      cout <<"M*N*K=" << M*N*K << "  != size = "<< size << "\n";
    MPI_Finalize();
    return 1;
  }

  if (adjEnts)
    keep_skins = true; // do not delete anything

  // determine m, n, k for processor rank
  int m,n,k;
  k = rank/(M*N);
  int leftover = rank%(M*N);
  n = leftover/M;
  m = leftover%M;
  if (rank==size-1)
    cout << "m, n, k for last rank:" << m << " " << n << " " << k << "\n";

  // so there are a total of M * A * blockSize elements in x direction (so M * A * blockSize + 1 verts in x direction)
  // so there are a total of N * B * blockSize elements in y direction (so N * B * blockSize + 1 verts in y direction)
  // so there are a total of K * C * blockSize elements in z direction (so K * C * blockSize + 1 verts in z direction)

  // there are ( M * A blockSize )      *  ( N * B * blockSize)      * (K * C * blockSize )    hexas
  // there are ( M * A * blockSize + 1) *  ( N * B * blockSize + 1 ) * (K * C * blockSize + 1) vertices
  // x is the first dimension that varies

  clock_t tt = clock();

  // used for nodes increments
  int q = (quadratic)? 2 : 1;
  // used for element increments
  int factor = (tetra)? 6 : 1;

  int NX = (q * M * A * blockSize + 1);
  int NY = (q * N * B * blockSize + 1);
  int nex = M * A * blockSize; // number of elements in x direction, used for global id on element
  int ney = N * B * blockSize; // number of elements in y direction  ....
  // int NZ = ( K * C * blockSize + 1); not used
  int blockSize1 = q*blockSize + 1;// used for vertices

  //int xstride = 1;
  int ystride = blockSize1;

  int zstride = blockSize1 * blockSize1;
  // generate the block at (a, b, c); it will represent a partition , it will get a partition tag

  Tag global_id_tag;
  mb->tag_get_handle("GLOBAL_ID", 1, MB_TYPE_INTEGER,
               global_id_tag);

  Tag part_tag;
  int dum_id=-1;
  mb->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER,
                 part_tag, MB_TAG_CREAT|MB_TAG_SPARSE, &dum_id);

  // create tags on vertices and cells, look in the list of options
  vector<Tag> intTags(intTagNames.size());
  vector<Tag> doubleTags(doubleTagNames.size());
  for (size_t i=0; i<intTagNames.size(); i++)
  {
    rval = mb->tag_get_handle(intTagNames[i].c_str(), 1, MB_TYPE_INTEGER, intTags[i],
        MB_TAG_CREAT|MB_TAG_DENSE, &dum_id);
    CHECKE("Can't create integer tag.");
  }
  double defval=0.;
  for (size_t i=0; i<doubleTagNames.size(); i++)
  {
    rval = mb->tag_get_handle(doubleTagNames[i].c_str(), 1, MB_TYPE_DOUBLE, doubleTags[i],
        MB_TAG_CREAT|MB_TAG_DENSE, &defval);
    CHECKE("Can't create double tag.");
  }
  for (int a=0; a<A; a++)
  {
    for (int b=0; b<B; b++)
    {
      for (int c=0; c<C; c++)
      {
        // we will generate (q*block+1)^3 vertices, and block^3 hexas; q is 1 for linear, 2 for quadratic
        // the global id of the vertices will come from m, n, k, a, b, c
        // x will vary from  m*A*q*block + a*q*block to m*A*q*block+(a+1)*q*block etc;
        int num_nodes = blockSize1*blockSize1*blockSize1;

        vector<double*> arrays;
        EntityHandle startv;
        rval = iface->get_node_coords(3, num_nodes, 0, startv, arrays);
        CHECKE("Can't get node coords.");

        // will start with the lower corner:
        int x = m*A*q*blockSize + a*q*blockSize;
        int y = n*B*q*blockSize + b*q*blockSize;
        int z = k*C*q*blockSize + c*q*blockSize;
        int ix=0;
        vector<int>  gids(num_nodes);
        Range verts(startv, startv+num_nodes-1);
        for (int kk=0; kk<blockSize1; kk++)
        {
          for (int jj=0; jj<blockSize1; jj++)
          {
            for (int ii=0; ii<blockSize1; ii++)
            {
              arrays[0][ix] = (double)x+ii;
              arrays[1][ix] = (double)y+jj;
              arrays[2][ix] = (double)z+kk;
              gids[ix] = 1 + (x+ii) + (y+jj) * NX + (z+kk) * (NX*NY) ;
              // set int tags, some nice values?
              EntityHandle v = startv + ix;
              for (size_t i=0; i<intTags.size(); i++)
              {
                int valv=gids[ix]/2+3 + i*1000;
                mb->tag_set_data(intTags[i], &v, 1, &valv);
              }
              ix++;
            }
          }
        }
        mb->tag_set_data(global_id_tag, verts, &gids[0]);
        int num_hexas = (blockSize)*(blockSize)*(blockSize);
        int num_el = num_hexas * factor;

        EntityHandle starte; // connectivity
        EntityHandle * conn;
        int num_v_per_elem = 8;
        if (quadratic)
        {
          num_v_per_elem = 27;
          rval = iface->get_element_connect(num_el, 27, MBHEX, 0, starte, conn);
        }
        else if (tetra)
        {
          num_v_per_elem = 4;
          rval = iface->get_element_connect(num_el, 4, MBTET, 0, starte, conn);
        }
        else
          rval = iface->get_element_connect(num_el, 8, MBHEX, 0, starte, conn);
        CHECKE("Can't get element  connectivity.");

        Range cells(starte, starte+num_el-1); // should be elements
        // fill  cells
        ix=0;
        // identify the elements at the lower corner, for their global ids
        int xe = m*A*blockSize + a*blockSize;
        int ye = n*B*blockSize + b*blockSize;
        int ze = k*C*blockSize + c*blockSize;
        gids.resize(num_el);
        int ie=0; // index now in the elements, for global ids
        for (int kk=0; kk<blockSize; kk++)
        {
          for (int jj=0; jj<blockSize; jj++)
          {
            for (int ii=0; ii<blockSize; ii++)
            {
              EntityHandle corner=startv + q * ii + q * jj * ystride + q * kk * zstride;
              gids[ie] = 1 + ((xe+ii) + (ye+jj) * nex + (ze+kk) * (nex*ney))*factor ; // 6 more for tetra

              EntityHandle eh = starte + ie;
              for (size_t i=0; i<doubleTags.size(); i++)
              {
                double valv=gids[ie]/30. + i*5000.;
                mb->tag_set_data(doubleTags[i], &eh, 1, &valv);
              }
              ie++;
              if (quadratic)
              {

  //                    4   ----- 19   -----  7
  //                .   |                 .   |
  //            16         25         18      |
  //         .          |          .          |
  //      5   ----- 17   -----  6             |
  //      |            12       | 23         15
  //      |                     |             |
  //      |     20      |  26   |     22      |
  //      |                     |             |
  //     13         21  |      14             |
  //      |             0   ----- 11   -----  3
  //      |         .           |         .
  //      |      8         24   |     10
  //      |  .                  |  .
  //      1   -----  9   -----  2
  //
                conn[ix]=    corner;
                conn[ix+1]=  corner + 2;
                conn[ix+2]=  corner + 2 + 2 * ystride;
                conn[ix+3]=  corner +     2 * ystride;
                conn[ix+4]=  corner                   + 2 * zstride;
                conn[ix+5]=  corner + 2               + 2 * zstride;
                conn[ix+6]=  corner + 2 + 2 * ystride + 2 * zstride;
                conn[ix+7]=  corner +     2 * ystride + 2 * zstride;
                conn[ix+8]=  corner + 1  ;                                         // 0-1
                conn[ix+9]=  corner + 2 +     ystride ;                            // 1-2
                conn[ix+10]= corner + 1 + 2 * ystride ;                            // 2-3
                conn[ix+11]= corner +         ystride;                             // 3-0
                conn[ix+12]= corner +                       zstride;               // 0-4
                conn[ix+13]= corner + 2 +                   zstride;               // 1-5
                conn[ix+14]= corner + 2 + 2 * ystride +     zstride;               // 2-6
                conn[ix+15]= corner +     2 * ystride +     zstride;               // 3-7
                conn[ix+16]= corner + 1 +               2 * zstride;               // 4-5
                conn[ix+17]= corner + 2 +     ystride + 2 * zstride;               //5-6
                conn[ix+18]= corner + 1 + 2 * ystride + 2 * zstride;               // 6-7
                conn[ix+19]= corner +         ystride + 2 * zstride;               // 4-7
                conn[ix+20]= corner + 1 +                   zstride;               // 0154
                conn[ix+21]= corner + 2 +     ystride +     zstride;               // 1265
                conn[ix+22]= corner + 1 + 2 * ystride +     zstride;               // 2376
                conn[ix+23]= corner +         ystride +     zstride;               // 0374
                conn[ix+24]= corner + 1 +     ystride              ;               // 0123
                conn[ix+25]= corner + 1 +     ystride + 2 * zstride;               // 4567
                conn[ix+26]= corner + 1 +     ystride +     zstride;               // center
                ix+=27;
              }
              else if (tetra)
              {

                //        E      H
                //     F     G
                //
                //        A     D
                //     B     C
                EntityHandle AA = corner;
                EntityHandle BB = corner + 1;
                EntityHandle CC = corner + 1 + ystride;
                EntityHandle D = corner +     ystride;
                EntityHandle E = corner +               zstride;
                EntityHandle F = corner + 1 +           zstride;
                EntityHandle G = corner + 1 + ystride + zstride;
                EntityHandle H = corner +     ystride + zstride;

                // tet EDHG
                conn[ix]    = E;
                conn[ix+1]  = D;
                conn[ix+2]  = H;
                conn[ix+3]  = G;

                // tet ABCF
                conn[ix+4]  = AA;
                conn[ix+5]  = BB;
                conn[ix+6]  = CC;
                conn[ix+7]  = F;

                // tet ADEF
                conn[ix+8]  = AA;
                conn[ix+9]  = D;
                conn[ix+10] = E;
                conn[ix+11] = F;

                // tet CGDF
                conn[ix+12] = CC;
                conn[ix+13] = G;
                conn[ix+14] = D;
                conn[ix+15] = F;

                // tet ACDF
                conn[ix+16] = AA;
                conn[ix+17] = CC;
                conn[ix+18] = D;
                conn[ix+19] = F;

                // tet DGEF
                conn[ix+20] = D;
                conn[ix+21] = G;
                conn[ix+22] = E;
                conn[ix+23] = F;
                ix+=24;
                for (int ff=0; ff<factor-1; ff++)
                {
                  gids[ie] = gids[ie-1]+1 ; // 6 more for tetra

                  eh = starte + ie;
                  for (size_t i=0; i<doubleTags.size(); i++)
                  {
                    double valv=gids[ie]/30. + i*5000.;
                    mb->tag_set_data(doubleTags[i], &eh, 1, &valv);
                  }
                  ie++;
                }
              }
              else // linear hex
              {
                conn[ix]=  corner;
                conn[ix+1]=corner + 1;
                conn[ix+2]=corner + 1 + ystride;
                conn[ix+3]=corner +     ystride;
                conn[ix+4]=corner +               zstride;
                conn[ix+5]=corner + 1 +           zstride;
                conn[ix+6]=corner + 1 + ystride + zstride;
                conn[ix+7]=corner +     ystride + zstride;
                ix+=8;
              }
            }
          }
        }
        EntityHandle part_set;
        rval = mb->create_meshset(MESHSET_SET, part_set);
        CHECKE("Can't create mesh set.");
        rval = mb->add_entities(part_set, cells);
        CHECKE("Can't add entities to set.");
        // if needed, add all edges and faces
        if (adjEnts)
        {
          // we need to update adjacencies now, because some elements are new
          rval = iface->update_adjacencies(starte, num_el, num_v_per_elem, conn);
          CHECKE("Can't update adjacencies");
          // generate all adj entities dimension 1 and 2 (edges and faces/ tri or qua
          Range edges, faces;
          rval = mb->get_adjacencies(cells, 1, true, edges,
              Interface::UNION);
          CHECKE("Can't get edges");
          rval = mb->get_adjacencies(cells, 2, true, faces,
              Interface::UNION);
          CHECKE("Can't get faces");
          rval = mb->add_entities(part_set, edges);
          CHECKE("Can't add edges to partition set.");
          rval = mb->add_entities(part_set, faces);
          CHECKE("Can't add faces to partition set.");
        }
        rval = mb->tag_set_data(global_id_tag, cells, &gids[0]);
        CHECKE("Can't set global ids to elements.");
        int part_num= a +  m*A + (b + n*B)*(M*A) + (c+k*C)*(M*A * N*B);
        rval = mb->tag_set_data(part_tag, &part_set, 1, &part_num);
        CHECKE("Can't set part tag on set");



      }
    }
  }

  if (0==rank)
  {
      std::cout << "generate local mesh:  "
            << (clock() - tt) / (double) CLOCKS_PER_SEC << " seconds" << std::endl;
      tt = clock();
  }

  /*// before merge locally
  rval = mb->write_file("test0.h5m", 0, ";;PARALLEL=WRITE_PART");
  CHECKE("Can't write in parallel, before merging");*/
  // after the mesh is generated on each proc, merge the vertices
  MergeMesh mm(mb);
  Range all3dcells;
  rval = mb->get_entities_by_dimension(0, 3, all3dcells);

  CHECKE("Can't get all 3d cells  elements.");

  Range verts;
  rval = mb->get_entities_by_dimension(0, 0, verts);

  CHECKE("Can't get all vertices.");
  if (A*B*C!=1) //  merge needed
  {
    if (newMergeMethod)
      rval = mm.merge_using_integer_tag( verts, global_id_tag);
    else
      rval = mm.merge_entities(all3dcells, 0.0001);
    CHECKE("Can't merge");
    if (0==rank)
    {
       std::cout << "merge locally:  "
             << (clock() - tt) / (double) CLOCKS_PER_SEC << " seconds" << std::endl;
       tt = clock();
    }
  }
  if (size>1)
  {
    ParallelComm* pcomm = ParallelComm::get_pcomm(mb, 0);
    if (pcomm==NULL)
    {
      pcomm = new ParallelComm( mb, MPI_COMM_WORLD );
    }
    rval = pcomm->resolve_shared_ents( 0, all3dcells, 3, 0 );
    CHECKE("Can't resolve shared ents");

    if (0==rank)
    {
       std::cout << "resolve shared entities:  "
             << (clock() - tt) / (double) CLOCKS_PER_SEC << " seconds" << std::endl;
       tt = clock();
    }

    if (!keep_skins) // default is to delete the 1- and 2-dimensional entities
    {
      // delete all quads and edges
      Range toDelete;
      rval =  mb->get_entities_by_dimension(0, 1, toDelete);
      CHECKE("Can't get edges");

      rval = mb->get_entities_by_dimension(0, 2, toDelete);
      CHECKE("Can't get faces");

      rval = pcomm->delete_entities(toDelete) ;
      CHECKE("Can't delete entities")

      if (0==rank)
      {
        std::cout << "delete edges and faces, and correct sharedEnts:  "
              << (clock() - tt) / (double) CLOCKS_PER_SEC << " seconds" << std::endl;
        tt = clock();
      }
    }
  }
  rval = mb->write_file(outFileName.c_str(), 0, ";;PARALLEL=WRITE_PART");
  CHECKE("Can't write in parallel");

  if (0==rank)
  {
    std::cout << "write file " << outFileName << " in "
          << (clock() - tt) / (double) CLOCKS_PER_SEC << " seconds" << std::endl;
    tt = clock();
  }

  MPI_Finalize();
  return 0;
}


