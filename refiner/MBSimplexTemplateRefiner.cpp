#include "MBSimplexTemplateRefiner.hpp"

#include "MBInterface.hpp"

#include <iostream>

/// Construct a template refiner.
MBSimplexTemplateRefiner::MBSimplexTemplateRefiner( MBInterface* mesh )
  : MBEntityRefiner( mesh )
{
}

/// Empty destructor for good form.
MBSimplexTemplateRefiner::~MBSimplexTemplateRefiner()
{
}

bool MBSimplexTemplateRefiner::refine_entity( MBEntityHandle entity )
{
  this->reset_heap_pointers();
  bool rval = true;
  switch ( this->mesh->type_from_handle( entity ) )
    {
    case MBVERTEX:
      this->refine_0_simplex();
      rval = false;
      break;
    case MBEDGE:
      rval = this->refine_1_simplex();
      break;
    case MBTRI:
      rval = this->refine_2_simplex();
      break;
    case MBQUAD:
      std::cerr << "Quadrilaterals not handled yet\n";
      rval = false;
      break;
    case MBPOLYGON:
      std::cerr << "Polygons not handled yet\n";
      rval = false;
      break;
    case MBTET:
      rval = this->refine_3_simplex();
      break;
    case MBPYRAMID:
      std::cerr << "Pyramid schemes not handled yet\n";
      rval = false;
      break;
    case MBPRISM:
      std::cerr << "Prisms not handled yet\n";
      rval = false;
      break;
    case MBKNIFE:
      std::cerr << "Knifahedra not handled yet\n";
      rval = false;
      break;
    case MBHEX:
      std::cerr << "Hexahedra not handled yet\n";
      rval = false;
      break;
    case MBPOLYHEDRON:
      std::cerr << "Polyhedra not handled yet\n";
      rval = false;
      break;
    case MBENTITYSET:
      std::cerr <<
        "How should entity sets be handled? We might iterate over its entries or skip it."
        "Must coordinate with MBMeshRefiner's loop over entities...\n";
      rval = false;
      break;
    case MBMAXTYPE:
      rval = false;
      break;
    }
  return rval;
}

/**\fn unsigned long MBSimplexTemplateRefiner::get_heap_size_bound( int max_recursions ) const
  *\brief Bound on the number of new vertices used to allocate the heap.
  *
  * This bound is based on a hexahedron that is divided into 48 tetrahedra (a point is added to
  * the center of each face so that compatible boundaries are guaranteed on neighboring hexahedra),
  * each of which has 4 edges.
  */

/**\brief "Refine" a vertex by passing it through to the output.
  *
  * FIXME: There is some question as to whether this should pass vertices through
  * since there does not appear to be a distinction between vertices as points
  * in space and vertices as degrees-of-freedom in a mesh (i.e. a vertex that is
  * treated as a lumped-parameter model).
  */
void MBSimplexTemplateRefiner::refine_0_simplex()
{
}

/**\brief Refine an edge.
  */
bool MBSimplexTemplateRefiner::refine_1_simplex()
{
  return true;
}

/**\brief Refine a triangle.
  */
bool MBSimplexTemplateRefiner::refine_2_simplex()
{
  return true;
}

/**\brief Refine a tetrahedron.
  */
bool MBSimplexTemplateRefiner::refine_3_simplex()
{
  return true;
}


