#include "MBSimplexTemplateRefiner.hpp"

#include "MBEdgeSizeEvaluator.hpp"
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
  const MBEntityHandle* conn;
  int num_nodes;
  if ( this->mesh->get_connectivity( entity, conn, num_nodes ) != MB_SUCCESS )
    {
    return false;
    }
  std::vector<double> entity_coords;
  entity_coords.resize( 6 * num_nodes );
  // Have to make num_nodes calls to get_coords() because we need xyz interleaved with rst coords.
  for ( int n = 0; n < num_nodes; ++n )
    {
    if ( this->mesh->get_coords( &conn[n], 1, &entity_coords[3 * n + 3] ) != MB_SUCCESS )
      {
      return false;
      }
    }
  // Still need to get tags.

  switch ( this->mesh->type_from_handle( entity ) )
    {
    case MBVERTEX:
      this->refine_0_simplex( &entity_coords[0], 0 ); // FIXME
      rval = false;
      break;
    case MBEDGE:
      rval = this->refine_1_simplex( 0,  &entity_coords[0], 0,  &entity_coords[6], 0 ); // FIXME
      break;
    case MBTRI:
      rval = this->refine_2_simplex( 0 ,   0, 0,   0, 0,   0, 0 ); // FIXME
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
void MBSimplexTemplateRefiner::refine_0_simplex( const double* v0, const void* t0 )
{
  (*this->output_functor)( v0, t0 );
  (*this->output_functor)( MBVERTEX );
}

/**\brief Refine an edge.
  */
bool MBSimplexTemplateRefiner::refine_1_simplex(
  int max_depth, const double* v0, const void* t0, const double* v1, const void* t1 )
{
  bool edge_code = false;

  double* midptc;
  void* midptt;

  if ( max_depth-- > 0 )
    {
    midptc = this->heap_coord_storage();
    midptt = this->heap_tag_storage();
    int i;
    // make valgrind happy
    //vtkstd::fill( midpt0, midpt0 + 6, 0. );
    for ( i = 0; i < 6; i++ )
      midptc[i] = ( v0[i] + v1[i] ) / 2.;

    this->edge_size_evaluator->evaluate_tags_at_midpoint( v0, t0, midptc, midptt, v1, t1 );
    edge_code = this->edge_size_evaluator->evaluate_edge( v0, t0, midptc, midptt, v1, t1 );
    }

  switch ( edge_code )
    {
    // No edges to subdivide
  case 0:
    (*this->output_functor)( v0, t0 );
    (*this->output_functor)( v1, t1 );
    (*this->output_functor)( MBEDGE );
    break ;

    // One edge to subdivide
  case 1:
    this->refine_1_simplex( max_depth, v0, t0, midptc, midptt );
    this->refine_1_simplex( max_depth, midptc, midptt, v1, t1 );
    break;
    }

  return edge_code;
}

/**\brief Refine a triangle.
  */
bool MBSimplexTemplateRefiner::refine_2_simplex(
  int max_depth, const double* v0, const void* t0, const double* v1, const void* t1, const double* v2, const void* t2 )
{
  int edgeCode = 0;

  double* midpt0c;
  double* midpt1c;
  double* midpt2c;
  void* midpt0t;
  void* midpt1t;
  void* midpt2t;

  if ( max_depth-- > 0 )
    {
    int i;
    for ( i = 0; i < 3; ++i )
      {

      }
    }

  return true;
}

/**\brief Refine a tetrahedron.
  */
bool MBSimplexTemplateRefiner::refine_3_simplex()
{
  return true;
}


