#include "MBSimplexTemplateTagAssigner.hpp"

#include "MBEdgeSizeEvaluator.hpp"
#include "MBInterface.hpp"
#include "MBSimplexTemplateRefiner.hpp"

#include <vector>

#include <math.h>

using namespace std;

/// Construct a template tag assigner.
MBSimplexTemplateTagAssigner::MBSimplexTemplateTagAssigner( MBSimplexTemplateRefiner* r )
{
  this->mesh_refiner = r;
  this->edge_size_evaluator = r->get_edge_size_evaluator();
}

/// Empty destructor for good form.
MBSimplexTemplateTagAssigner::~MBSimplexTemplateTagAssigner()
{
}

/**\brief Given endpoint coordinates and tag values plus midpoint coordinates, compute midpoint tag values.
  *
  * Normally, this function will be invoked by the MBEntityRefiner before evaluate_edge is called.
  * However, if evaluate_edge() changes the parametric coordinates of the midpoint,
  * it should call evaluate_tags_at_midpoint() again to update any tag values;
  * that is why this function is a member of MBEdgeSizeEvaluator and not MBEntityRefiner.
  *
  * @param[in] c0 Pointer to endpoint 0 coordinates. The parametric coordinates (3) are followed by world coordinates (3).
  * @param[in] t0 Pointer to endpoint 0 tag values.
  * @param[in] cm Pointer to midpoint coordinates. The parametric coordinates (3) are followed by world coordinates (3).
  * @param[out] tm Pointer to midpoint tag values.
  * @param[in] c1 Pointer to endpoint 1 coordinates. The parametric coordinates (3) are followed by world coordinates (3).
  * @param[in] t1 Pointer to endpoint 1 tag values.
  */
void MBSimplexTemplateTagAssigner::operator () (
  const double* c0, const void* t0, MBEntityHandle h0,
  const double* cm, void* tm, 
  const double* c1, const void* t1, MBEntityHandle h1 )
{
  double c0m_squared = 0.;
  double c01_squared = 0.;
  for ( int i = 0; i < 3; ++i )
    {
    double tmp = cm[i] - c0[i];
    c0m_squared += tmp * tmp;
    tmp = c1[i] - c0[i];
    c01_squared += tmp * tmp;
    }
  double lambda = sqrt( c0m_squared / c01_squared );
  double one_minus_lambda = 1. - lambda;

  MBDataType data_type;
  int tag_size;
  int num_components;
  int num_tags = this->edge_size_evaluator->get_number_of_vertex_tags();
  MBTag tag_handle;
  int tag_offset;
  for ( int i = 0; i < num_tags; ++i )
    {
    this->edge_size_evaluator->get_vertex_tag( i, tag_handle, tag_offset );
    this->mesh_refiner->get_mesh()->tag_get_data_type( tag_handle, data_type );
    this->mesh_refiner->get_mesh()->tag_get_size( tag_handle, tag_size );
    
    switch ( data_type )
      {
      case MB_TYPE_DOUBLE:
        {
        num_components = tag_size / sizeof( double );
        double* t0i = (double*) ( (char*)t0 + tag_offset );
        double* tmi = (double*) ( (char*)tm + tag_offset );
        double* t1i = (double*) ( (char*)t1 + tag_offset );
        for ( int i = 0; i < num_components; ++ i )
          tmi[i] = one_minus_lambda * t0i[i] + lambda * t1i[i];
        }
        break;
      default:
        memcpy( (char*)tm + tag_offset, (char*)( h0 < h1 ? t0 : t1 ) + tag_offset, tag_size );
        break;
      }
    }
}

void MBSimplexTemplateTagAssigner::operator () ( const void* t0,
                                                 const void* t1,
                                                 const void* t2,
                                                 void* tp )
{
  (void)t0;
  (void)t1;
  (void)t2;
  (void)tp;
}
