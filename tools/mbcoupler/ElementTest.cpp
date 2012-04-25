#include "TestUtil.hpp"
#include "ElemUtil.hpp"
#include <iostream>

using namespace moab;

void test_tet();
void test_hex();
void test_spectral_hex();

int main()
{
  int rval = 0;
  rval += RUN_TEST(test_tet);
  rval += RUN_TEST(test_hex);
  rval += RUN_TEST(test_spectral_hex);
  return rval;
}

void test_tet() {
  moab::Element::LinearTet tet;
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
  moab::ErrorCode rval = mb->load_mesh("spectral.h5m");
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

  }
  std::cout << "success...\n";
  
  return;
}

