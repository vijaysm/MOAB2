#include "TestUtil.hpp"
#include "ElemUtil.hpp"
#include <iostream>

using namespace moab;

void test_tet();
void test_hex();
void test_spectral_hex();
void test_spectral_quad();

int main()
{
  int rval = 0;
  rval += RUN_TEST(test_tet);
  rval += RUN_TEST(test_hex);
  rval += RUN_TEST(test_spectral_hex);
  rval += RUN_TEST(test_spectral_quad);
  return rval;
}

void test_tet() {
  // test constructor on one simple tet
  std::vector<CartVect> verts;
  verts.push_back(CartVect(0.,0.,0.2));
  verts.push_back(CartVect(1.1,0.,0.));
  verts.push_back(CartVect(0.,1.13,0.));
  verts.push_back(CartVect(0.,0.,1.));
  moab::Element::LinearTet tet(verts);
  CartVect p(0.25, 0.25, 0.25);
  CartVect res = tet.ievaluate(p);
  CartVect p1 = tet.evaluate(res);
  double d = (p1-p).length_squared();
  CHECK_REAL_EQUAL(0., d, 1.e-10);

  return;
}// test_tet()

void test_hex() {
  moab::Element::LinearHex hex;
}// test_hex()
#include "moab/Core.hpp"
#include "moab/Range.hpp"

void test_spectral_hex()
{
  // first load a model that has spectral elements
  moab::Core *mb = new moab::Core();
  std::string meshFile = STRINGIFY(SRCDIR) "/spectral.h5m";
  moab::ErrorCode rval = mb->load_mesh(meshFile.c_str());
  if (moab::MB_SUCCESS != rval) return ;

  // get the ent set with SEM_DIMS tag
  moab::Range spectral_sets;
  moab::Tag  sem_tag;
  rval = mb->tag_get_handle("SEM_DIMS", 3, moab::MB_TYPE_INTEGER, sem_tag);
  if (moab::MB_SUCCESS != rval)
  {
    std::cout<< "can't find tag, no spectral set\n";
    return ;
  }
  rval = mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &sem_tag, NULL, 1, spectral_sets);
  if (moab::MB_SUCCESS != rval || spectral_sets.empty())
  {
    std::cout<< "can't get sem set\n";
    return ;
  }

  moab::Range ents;

  int sem_dims[3];
  moab::EntityHandle firstSemSet=spectral_sets[0];
  rval = mb->tag_get_data(sem_tag, &firstSemSet, 1, (void*)sem_dims);
  if (moab::MB_SUCCESS != rval) return ;

  rval = mb->get_entities_by_dimension(firstSemSet, 3, ents);
  if (moab::MB_SUCCESS != rval) return ;
  std::cout << "Found " << ents.size() << " " << 3 << "-dimensional entities:" << std::endl;
  
  if (sem_dims[0]!=sem_dims[1] || sem_dims[0] != sem_dims[2])
  {
    std::cout << " dimensions are different. bail out\n";
    return;
  }

  // get the SEM_X ...tags  
  moab::Tag xm1Tag, ym1Tag, zm1Tag;
  int ntot = sem_dims[0]*sem_dims[1]*sem_dims[2];
  rval = mb->tag_get_handle("SEM_X", ntot, moab::MB_TYPE_DOUBLE, xm1Tag); 
  if (moab::MB_SUCCESS != rval) 
  {
     std::cout << "can't get xm1tag \n";
     return;
  }
  rval = mb->tag_get_handle("SEM_Y", ntot, moab::MB_TYPE_DOUBLE, ym1Tag); 
  if (moab::MB_SUCCESS != rval) 
  {
     std::cout << "can't get ym1tag \n";
     return;
  }
  rval = mb->tag_get_handle("SEM_Z", ntot, moab::MB_TYPE_DOUBLE, zm1Tag); 
  if (moab::MB_SUCCESS != rval) 
  {
     std::cout << "can't get zm1tag \n";
     return;
  }
  moab::Tag velTag;
  
  rval = mb->tag_get_handle("VX", ntot, moab::MB_TYPE_DOUBLE, velTag); 
  if (moab::MB_SUCCESS != rval) 
  {
     std::cout << "can't get veltag \n";
     return;
  }
  moab::Element::SpectralHex specHex(sem_dims[0] );

 // compute the data for some elements 
  for (moab::Range::iterator rit=ents.begin(); rit!=ents.end(); rit++)
  { 
  // get the tag pointers to the internal storage for xm1, to not copy the values
     moab::EntityHandle eh= *rit;
     const double * xval;
     const double * yval;
     const double * zval;
     rval = mb-> tag_get_by_ptr(xm1Tag, &eh, 1,(const void **) &xval );
     if (moab::MB_SUCCESS != rval)
     {
       std::cout << "can't get xm1 values \n";
       return;
     }
     rval = mb-> tag_get_by_ptr(ym1Tag, &eh, 1, (const void **)&yval );
     if (moab::MB_SUCCESS != rval)
     {
       std::cout << "can't get ym1 values \n";
       return;
     }
     rval = mb-> tag_get_by_ptr(zm1Tag, &eh, 1, (const void **)&zval );
     if (moab::MB_SUCCESS != rval)
     {
       std::cout << "can't get zm1 values \n";
       return;
     }
     if (rit==ents.begin())
     {
        std::cout << " xm1 for first element: \n";
        for (int i=0; i< ntot; i++)
          std::cout << " " << xval[i] ; 
        std::cout << "\n";
     }
     specHex.set_gl_points((double*)xval, (double*)yval, (double*)zval);
     // first evaluate a point, then inverse it to see if we get the same thing
     moab::CartVect rst(0.1, -0.1, 0.5);
     moab::CartVect pos = specHex.evaluate(rst);
     moab::CartVect inverse = specHex.ievaluate(pos);
     std::cout << "difference" << rst-inverse << "\n";
     Matrix3 jac=specHex.jacobian(rst);
     std::cout<< "jacobian: \n" << jac << " \n determinant: " << jac.determinant() << "\n";
     // evaluate vx at some point
     const double * vx;
     rval = mb-> tag_get_by_ptr(velTag, &eh, 1, (const void **)&vx );
     if (moab::MB_SUCCESS != rval)
     {
       std::cout << "can't get vel values \n";
       return;
     }
     double vel1 = specHex.evaluate_scalar_field(rst, vx);
     std::cout << "velocity: " << vel1 << "\n";
     // compute integral over vx:
     double integral = specHex.integrate_scalar_field(vx);
     std::cout << "integral over vx: " << integral << "\n";

  }
  std::cout << "success...\n";
  
  return;
}
void test_spectral_quad()
{
  // first load a model that has spectral elements
  moab::Core *mb = new moab::Core();
  // use the grid on Sphere from mbcslam folder
  std::string meshFile = STRINGIFY(SRCDIR) "/../mbcslam/eulerHomme.vtk";
  moab::ErrorCode rval = mb->load_mesh(meshFile.c_str());
  if (moab::MB_SUCCESS != rval) return ;

  // for each element, compute the gl points and project them on sphere
  // the radius is the same as the one from intersection test on sphere
  //double R = 6. * sqrt(3.) / 2; // input



  moab::Range ents;


  rval = mb->get_entities_by_type(0, moab::MBQUAD, ents);// get all quads
  if (moab::MB_SUCCESS != rval) return ;
  std::cout << "Found " << ents.size() << " " << 2 << "-dimensional entities:" << std::endl;

//
  int NP =4; // test this....
  moab::Element::SpectralQuad specQuad(NP);

 // compute the gl points for some elements
  for (moab::Range::iterator rit = ents.begin(); rit != ents.end(); rit++)
  {

    const moab::EntityHandle * conn4 = NULL;
    int num_nodes = 0;
    rval = mb->get_connectivity(*rit, conn4, num_nodes);
    if (moab::MB_SUCCESS != rval)
    {
      std::cout << "can't get conectivity for quad \n";
      return;
    }
    assert(num_nodes==4);

    std::vector<CartVect> verts(4);
    rval = mb->get_coords(conn4, num_nodes, &(verts[0][0]));
    if (moab::MB_SUCCESS != rval)
    {
      std::cout << "can't get coords for quad \n";
      return;
    }


     specQuad.set_vertices(verts);
     specQuad.compute_gl_positions();
     // do something with the gl positions, project them on a sphere, and create another mesh?
     if (rit==ents.begin())
      {
         std::cout << " gl points for first element: \n";
         int size;
         double * xyz[3];
         specQuad.get_gl_points(xyz[0], xyz[1], xyz[2], size);
         for (int i=0; i<size; i++)
           std::cout << xyz[0][i] << " " <<  xyz[1][i] << " " <<  xyz[2][i] << "\n";
      }

     // project them on a sphere, and create another mesh with it?
  }
  std::cout << "success...\n";

  return;
}
