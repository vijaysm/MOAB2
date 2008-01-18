#include "MBSimplexTemplateRefiner.hpp"

#include "MBEdgeSizeEvaluator.hpp"
#include "MBInterface.hpp"

#include <iostream>
#include <stack>

// Static arrays holding parametric coordinates of element vertices
static double MBVertexParametric[] = { 0., 0., 0. };
static double MBEdgeParametric[]   = { 0., 0., 0.,   1., 0., 0. };
static double MBTriParametric[]    = { 0., 0., 0.,   1., 0., 0.,   0., 1., 0. };
static double MBTetParametric[]    = { 0., 0., 0.,   1., 0., 0.,   0., 1., 0.,   0., 0., 1. };

#ifdef MB_DEBUG_TESSELLATOR
#  define MB_TESSELLATOR_INCR_CASE_COUNT(cs) this->case_counts[cs]++
#  define MB_TESSELLATOR_INCR_SUBCASE_COUNT(cs,sc) this->subcase_counts[cs][sc]++
#else // MB_DEBUG_TESSELLATOR
#  define MB_TESSELLATOR_INCR_CASE_COUNT(cs)
#  define MB_TESSELLATOR_INCR_SUBCASE_COUNT(cs,sc)
#endif // MB_DEBUG_TESSELLATOR

/// Construct a template refiner.
MBSimplexTemplateRefiner::MBSimplexTemplateRefiner( MBInterface* mesh )
  : MBEntityRefiner( mesh )
{
  this->tag_assigner = new MBSimplexTemplateTagAssigner( this );
  this->corner_coords.resize( 6 * 8 ); // Hex has 8 verts w/ 6 coordinates each
  this->corner_tags.resize( 8 ); // Hex has 8 verts (this is a pointer, not the actual tag data)
}

/// Empty destructor for good form.
MBSimplexTemplateRefiner::~MBSimplexTemplateRefiner()
{
  delete this->tag_assigner;
}

/**\brief Stream a single mesh entity through the refiner.
  */
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
  this->corner_coords.resize( 6 * num_nodes );
  this->corner_tags.resize( num_nodes );

  MBEntityType etyp = this->mesh->type_from_handle( entity );
  // Have to make num_nodes calls to get_coords() because we need xyz interleaved with rst coords.
  MBTag tag_handle;
  int tag_offset;
  void* tag_data;
  for ( int n = 0; n < num_nodes; ++ n )
    {
    if ( this->mesh->get_coords( &conn[n], 1, &corner_coords[6 * n + 3] ) != MB_SUCCESS )
      {
      return false;
      }
    tag_data = this->heap_tag_storage();
    for ( int i = 0; i < this->edge_size_evaluator->get_number_of_vertex_tags(); ++ i )
      {
      this->edge_size_evaluator->get_vertex_tag( i, tag_handle, tag_offset );
      if ( this->mesh->tag_get_data( tag_handle, &conn[n], 1, (char*)tag_data + tag_offset ) != MB_SUCCESS )
        {
        return false;
        }
      }
    this->corner_tags[n] = tag_data;
    }

  switch ( etyp )
    {
    case MBVERTEX:
      this->assign_parametric_coordinates( 1, MBVertexParametric, &this->corner_coords[0] );
      this->refine_0_simplex( &this->corner_coords[0], this->corner_tags[0] );
      rval = false;
      break;
    case MBEDGE:
      this->assign_parametric_coordinates( 2, MBEdgeParametric, &this->corner_coords[0] );
      rval = this->refine_1_simplex( this->maximum_number_of_subdivisions,
        &this->corner_coords[0], this->corner_tags[0],  &this->corner_coords[6], this->corner_tags[1] );
      break;
    case MBTRI:
      this->assign_parametric_coordinates( 3, MBTriParametric, &this->corner_coords[0] );
      rval = this->refine_2_simplex( this->maximum_number_of_subdivisions, 7,
        &this->corner_coords[ 0], this->corner_tags[0],
        &this->corner_coords[ 6], this->corner_tags[1],
        &this->corner_coords[12], this->corner_tags[2] );
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
      this->assign_parametric_coordinates( 4, MBTetParametric, &this->corner_coords[0] );
      rval = this->refine_3_simplex( 0, 0, 0, 0, 0, 0, 0, 0, 0 ); // FIXME
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

/**\brief Set the function object used to decide which tag values an edge or face midpoint is assigned.
  *
  * This will change the tag assigner's edge size evaluator to match the refiner's.
  * @param[in] ta The new tag assigner. This must be non-NULL.
  * @retval True if the tag assigner was changed and false otherwise.
  */
bool MBSimplexTemplateRefiner::set_tag_assigner( MBSimplexTemplateTagAssigner* ta )
{ 
  if ( ! ta || ta == this->tag_assigner )
    return false;

  this->tag_assigner = ta; 
  this->tag_assigner->set_edge_size_evaluator( this->edge_size_evaluator );
  return true;
}


/**\brief Set the function object used to decide whether an edge is subdivided or not.
  */
bool MBSimplexTemplateRefiner::set_edge_size_evaluator( MBEdgeSizeEvaluator* es ) 
{ 
  if ( this->MBEntityRefiner::set_edge_size_evaluator( es ) )
    {
    this->tag_assigner->set_edge_size_evaluator( es );
    return true;
    }
  return false;
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
  int i0, i1;

  if ( max_depth-- > 0 )
    {
    midptc = this->heap_coord_storage();
    midptt = this->heap_tag_storage();
    int i;
    // make valgrind happy
    //vtkstd::fill( midpt0, midpt0 + 6, 0. );
    for ( i = 0; i < 6; i++ )
      midptc[i] = ( v0[i] + v1[i] ) / 2.;

    (*this->tag_assigner)( v0, t0, i0, midptc, midptt, v1, t1, i1 );
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
  int max_depth, int move, const double* v0, const void* t0, const double* v1, const void* t1, const double* v2, const void* t2 )
{
  int edge_code = 0;

  double* midpt0c;
  double* midpt1c;
  double* midpt2c;
  void* midpt0t;
  void* midpt1t;
  void* midpt2t;
  int i0, i1, i2;

  if ( max_depth-- > 0 )
    {
    int i;
    midpt0c = this->heap_coord_storage();
    midpt1c = this->heap_coord_storage();
    midpt2c = this->heap_coord_storage();
    midpt0t = this->heap_tag_storage();
    midpt1t = this->heap_tag_storage();
    midpt2t = this->heap_tag_storage();
    for ( i = 0; i < 6; ++i )
      {
      midpt0c[i] = ( v0[i] + v1[i] ) / 2.;
      midpt1c[i] = ( v1[i] + v2[i] ) / 2.;
      midpt2c[i] = ( v2[i] + v0[i] ) / 2.;
      }
    (*this->tag_assigner)( v0, t0, i0, midpt0c, midpt0t, v1, t1, i1 );
    (*this->tag_assigner)( v1, t1, i1, midpt1c, midpt1t, v2, t2, i2 );
    (*this->tag_assigner)( v2, t2, i2, midpt2c, midpt2t, v0, t0, i0 );
    if ( ( move & 1 ) && this->edge_size_evaluator->evaluate_edge( v0, t0, midpt0c, midpt0t, v1, t1 ) )
      edge_code += 1;
    if ( ( move & 2 ) && this->edge_size_evaluator->evaluate_edge( v1, t1, midpt1c, midpt1t, v2, t2 ) )
      edge_code += 2;
    if ( ( move & 4 ) && this->edge_size_evaluator->evaluate_edge( v2, t2, midpt2c, midpt2t, v0, t0 ) )
      edge_code += 4;
    }

  switch ( edge_code )
    {
    // No edges to subdivide
  case 0:
    (*this->output_functor)( v0, t0 );
    (*this->output_functor)( v1, t1 );
    (*this->output_functor)( v2, t2 );
    (*this->output_functor)( MBTRI );
    break ;

    // One edge to subdivide
  case 1:
    this->refine_2_simplex( max_depth, move | 2, v0, t0, midpt0c, midpt0t, v2, t2 );
    this->refine_2_simplex( max_depth, move | 4, midpt0c, midpt0t, v1, t1, v2, t2 );
    break;
  case 2:
    this->refine_2_simplex( max_depth, move | 4, v0, t0, v1, t1, midpt1c, midpt1t );
    this->refine_2_simplex( max_depth, move | 1, v0, t0, midpt1c, midpt1t, v2, t2 );
    break;
  case 4:
    this->refine_2_simplex( max_depth, move | 2, v0, t0, v1, t1, midpt2c, midpt2t );
    this->refine_2_simplex( max_depth, move | 1, midpt2c, midpt2t, v1, t1, v2, t2 );
    break;

    // Two edges to subdivide
  case 3:
    this->refine_2_simplex( max_depth, move | 4, midpt0c, midpt0t, v1, t1, midpt1c, midpt1t );
    if ( this->compare_Hopf_cross_string_dist( v2, midpt0c, v0, midpt1c ) )
      {
      this->refine_2_simplex( max_depth, move | 5, midpt0c, midpt0t, midpt1c, midpt1t,   v2,      t2    );
      this->refine_2_simplex( max_depth, move | 2,   v0,      t0,    midpt0c, midpt0t,   v2,      t2    );
      }
    else                                         
      {
      this->refine_2_simplex( max_depth, move | 6,   v0,      t0,    midpt0c, midpt0t, midpt1c, midpt1t );
      this->refine_2_simplex( max_depth, move | 1,   v0,      t0,    midpt1c, midpt1t,   v2,      t2    );
      }
    break;
  case 5:
    this->refine_2_simplex( max_depth, move | 2, v0, t0, midpt0c, midpt0t, midpt2c, midpt2t );
    if ( this->compare_Hopf_cross_string_dist( v2, midpt0c, v1, midpt2c ) )
      {
      this->refine_2_simplex( max_depth, move | 4, midpt0c, midpt0t,   v1,      t1,       v2,      t2   );
      this->refine_2_simplex( max_depth, move | 3, midpt2c, midpt2t, midpt0c, midpt0t,    v2,      t2   );
      }
    else
      {
      this->refine_2_simplex( max_depth, move | 6, midpt0c, midpt0t,   v1,      t1,    midpt2c, midpt2t );
      this->refine_2_simplex( max_depth, move | 1, midpt2c, midpt2t,   v1,      t1,       v2,      t2   );
      }
    break;
  case 6:
    this->refine_2_simplex( max_depth, move | 1, midpt2c, midpt2t, midpt1c, midpt1t, v2, t2 );
    if ( this->compare_Hopf_cross_string_dist( v0, midpt1c, v1, midpt2c ) )
      {
      this->refine_2_simplex( max_depth, move | 3,   v0,      t0,    midpt1c, midpt1t, midpt2c, midpt2t );
      this->refine_2_simplex( max_depth, move | 4,   v0,      t0,      v1,      t1,    midpt1c, midpt1t );
      }
    else
      {
      this->refine_2_simplex( max_depth, move | 2,   v0,      t0,      v1,      t1,    midpt2c, midpt2t );
      this->refine_2_simplex( max_depth, move | 5, midpt2c, midpt2t,   v1,      t1,    midpt1c, midpt1t );
      }
    break;

    // Three edges to subdivide
  case 7:
    this->refine_2_simplex( max_depth,        7, midpt0c, midpt0t, midpt1c, midpt1t, midpt2c, midpt2t );
    this->refine_2_simplex( max_depth, move | 2,   v0   ,   t0   , midpt0c, midpt0t, midpt2c, midpt2t );
    this->refine_2_simplex( max_depth, move | 4, midpt0c, midpt0t,   v1   ,   t1   , midpt1c, midpt1t );
    this->refine_2_simplex( max_depth, move | 1, midpt2c, midpt2t, midpt1c, midpt1t,   v2   ,   t2    );
    break;
    }

  return true;
}

/**\brief Refine a tetrahedron.
  */
bool MBSimplexTemplateRefiner::refine_3_simplex( int max_depth,
                                                 double* v0, void* t0, 
                                                 double* v1, void* t1, 
                                                 double* v2, void* t2,
                                                 double* v3, void* t3 )
{
  bool edge_code = false;

  double* midpt0c;
  double* midpt1c;
  double* midpt2c;
  double* midpt3c;
  double* midpt4c;
  double* midpt5c;

  void* midpt0t;
  void* midpt1t;
  void* midpt2t;
  void* midpt3t;
  void* midpt4t;
  void* midpt5t;

  int i0, i1, i2, i3;

  if ( max_depth-- > 0 )
    {
    midpt0c = this->heap_coord_storage();
    midpt1c = this->heap_coord_storage();
    midpt2c = this->heap_coord_storage();
    midpt3c = this->heap_coord_storage();
    midpt4c = this->heap_coord_storage();
    midpt5c = this->heap_coord_storage();

    midpt0t = this->heap_tag_storage();
    midpt1t = this->heap_tag_storage();
    midpt2t = this->heap_tag_storage();
    midpt3t = this->heap_tag_storage();
    midpt4t = this->heap_tag_storage();
    midpt5t = this->heap_tag_storage();

    for ( int i = 0; i < 6; ++ i )
      {
      midpt0c[i] = ( v0[i] + v1[i] ) * .5;
      midpt1c[i] = ( v1[i] + v2[i] ) * .5;
      midpt2c[i] = ( v2[i] + v0[i] ) * .5;
      midpt3c[i] = ( v0[i] + v3[i] ) * .5;
      midpt4c[i] = ( v1[i] + v3[i] ) * .5;
      midpt5c[i] = ( v2[i] + v3[i] ) * .5;
      }

    (*this->tag_assigner)( v0, t0, i0, midpt0c, midpt0t, v1, t1, i1 );
    (*this->tag_assigner)( v1, t1, i1, midpt1c, midpt1t, v2, t2, i2 );
    (*this->tag_assigner)( v2, t2, i2, midpt2c, midpt2t, v0, t0, i0 );
    (*this->tag_assigner)( v0, t0, i0, midpt3c, midpt3t, v3, t3, i3 );
    (*this->tag_assigner)( v1, t1, i1, midpt4c, midpt4t, v3, t3, i3 );
    (*this->tag_assigner)( v2, t2, i2, midpt5c, midpt5t, v3, t3, i3 );

    if ( this->edge_size_evaluator->evaluate_edge( v0, t0, midpt0c, midpt0t, v1, t1 ) )
      edge_code |=  1;
    if ( this->edge_size_evaluator->evaluate_edge( v1, t1, midpt1c, midpt1t, v2, t2 ) )
      edge_code |=  2;
    if ( this->edge_size_evaluator->evaluate_edge( v2, t2, midpt2c, midpt2t, v0, t0 ) )
      edge_code |=  4;
    if ( this->edge_size_evaluator->evaluate_edge( v0, t0, midpt3c, midpt3t, v3, t3 ) )
      edge_code |=  8;
    if ( this->edge_size_evaluator->evaluate_edge( v1, t1, midpt4c, midpt4t, v3, t3 ) )
      edge_code |= 16;
    if ( this->edge_size_evaluator->evaluate_edge( v2, t2, midpt5c, midpt5t, v3, t3 ) )
      edge_code |= 32;
    }

  double edge_length2[6];
  edge_length2[0]
    = edge_length2[1]
    = edge_length2[2]
    = edge_length2[3]
    = edge_length2[4]
    = edge_length2[5]
    = 0;

  for ( int c = 0; c < 3; ++ c )
    {
    double tmp;
    tmp = v1[c] - v0[c];
    edge_length2[0] += tmp * tmp;
    tmp = v2[c] - v1[c];
    edge_length2[1] += tmp * tmp;
    tmp = v2[c] - v0[c];
    edge_length2[2] += tmp * tmp;
    tmp = v3[c] - v0[c];
    edge_length2[3] += tmp * tmp;
    tmp = v3[c] - v1[c];
    edge_length2[4] += tmp * tmp;
    tmp = v3[c] - v2[c];
    edge_length2[5] += tmp * tmp;
    }

  if ( ! edge_code )
    {
    // No edges to subdivide
    (*this->output_functor)( v0, t0 );
    (*this->output_functor)( v1, t1 );
    (*this->output_functor)( v2, t2 );
    (*this->output_functor)( v3, t3 );
    (*this->output_functor)( MBTET );

    return false;
    }

  double* facept0c;
  double* facept1c;
  double* facept2c;
  double* facept3c;
  facept0c = this->heap_coord_storage();
  facept1c = this->heap_coord_storage();
  facept2c = this->heap_coord_storage();
  facept3c = this->heap_coord_storage();
  double* vertex_coords[14] = { v0, v1, v2, v3, 
				midpt0c, midpt1c, midpt2c, 
				midpt3c, midpt4c, midpt5c,
                                facept0c, facept1c, facept2c, facept3c };

  void* facept0t;
  void* facept1t;
  void* facept2t;
  void* facept3t;
  facept0t = this->heap_tag_storage();
  facept1t = this->heap_tag_storage();
  facept2t = this->heap_tag_storage();
  facept3t = this->heap_tag_storage();
  void* vertex_tags[14] = { t0, t1, t2, t3, 
                            midpt0t, midpt1t, midpt2t, 
                            midpt3t, midpt4t, midpt5t,
                            facept0t, facept1t, facept2t, facept3t };
  
  // Generate tetrahedra that are compatible except when edge
  // lengths are equal on indeterminately subdivided faces.
  double* permuted_coords[14];
  void* permuted_tags[14];
  double permlen[6]; // permuted edge lengths
  int C = MBSimplexTemplateRefiner::template_index[edge_code][0];
  int P = MBSimplexTemplateRefiner::template_index[edge_code][1];
  
  // 1. Permute the tetrahedron into our canonical configuration
  for ( int i = 0; i < 4; ++ i )
    {
    permuted_coords[i] = vertex_coords[MBSimplexTemplateRefiner::permutations_from_index[P][i]];
    permuted_tags[i] = vertex_tags[MBSimplexTemplateRefiner::permutations_from_index[P][i]];
    }

  for ( int i = 4 ; i < 10; ++ i )
    {
    // permute edges too
    permuted_coords[i] = vertex_coords[MBSimplexTemplateRefiner::permutations_from_index[P][i]];
    permuted_tags[i] = vertex_tags[MBSimplexTemplateRefiner::permutations_from_index[P][i]];
    permlen[i-4]  = edge_length2[MBSimplexTemplateRefiner::permutations_from_index[P][i] - 4];
    }
  // Add our local (heap) storage for face point coordinates to the list.
  permuted_coords[10] = facept0c;
  permuted_coords[11] = facept1c;
  permuted_coords[12] = facept2c;
  permuted_coords[13] = facept3c;
  permuted_tags[10] = facept0t;
  permuted_tags[11] = facept1t;
  permuted_tags[12] = facept2t;
  permuted_tags[13] = facept3t;

  int comparison_bits;
  std::stack<int*> output_tets;
  std::stack<int*> output_perm;
  std::stack<int>  output_sign;

  // cout << "Case " << C << "  Permutation " << P << endl;
  // 2. Generate tetrahedra based on the configuration.
  //    Note that case 0 is handled above (edgeCode == 0).
  
  switch ( C )
    {
    case 1: // Ruprecht-Müller Case 1
      MB_TESSELLATOR_INCR_CASE_COUNT(0);
      output_tets.push( MBSimplexTemplateRefiner::templates + 0 );
      output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
      output_sign.push( 1 );
      MB_TESSELLATOR_INCR_SUBCASE_COUNT(0,0);
      break;
    case 2: // Ruprecht-Müller Case 2a
      comparison_bits = 
        ( permlen[0] <= permlen[1] ? 1 : 0 ) | ( permlen[0] >= permlen[1] ? 2 : 0 ) |
        0;
      if ( ( comparison_bits & 3 ) == 3 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[10][i] = ( permuted_coords[0][i] + permuted_coords[2][i] ) * .375 + permuted_coords[1][i] * .25;
          }
	(*this->tag_assigner)( t0, t2, t1, permuted_tags[10] );
        }
      MB_TESSELLATOR_INCR_CASE_COUNT(1);
      output_tets.push( MBSimplexTemplateRefiner::templates + 9 );
      output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
      output_sign.push( 1 );
      MB_TESSELLATOR_INCR_SUBCASE_COUNT(1,0);
      switch ( comparison_bits )
        {
        case 2: // 0>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 14 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(1,1);
          break;
        case 1: // 1>0
          output_tets.push( MBSimplexTemplateRefiner::templates + 14 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[13] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(1,2);
          break;
        case 3: // 0=1
          output_tets.push( MBSimplexTemplateRefiner::templates + 23 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(1,3);
          break;
        }
      break;
    case 3: // Ruprecht-Müller Case 2b
      MB_TESSELLATOR_INCR_CASE_COUNT(2);
      output_tets.push( MBSimplexTemplateRefiner::templates + 40 );
      output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
      output_sign.push( 1 );
      MB_TESSELLATOR_INCR_SUBCASE_COUNT(2,0);
      break;
    case 4: // Ruprecht-Müller Case 3a
      comparison_bits = 
        ( permlen[0] <= permlen[3] ? 1 : 0 ) | ( permlen[0] >= permlen[3] ? 2 : 0 ) |
        ( permlen[2] <= permlen[3] ? 4 : 0 ) | ( permlen[2] >= permlen[3] ? 8 : 0 ) |
        ( permlen[0] <= permlen[2] ? 16 : 0 ) | ( permlen[0] >= permlen[2] ? 32 : 0 ) |
        0;
      if ( ( comparison_bits & 3 ) == 3 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[11][i] = ( permuted_coords[1][i] + permuted_coords[3][i] ) * .375 + permuted_coords[0][i] * .25;
          }
	(*this->tag_assigner)( t1, t3, t0, permuted_tags[11] );
        }
      if ( ( comparison_bits & 12 ) == 12 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[13][i] = ( permuted_coords[2][i] + permuted_coords[3][i] ) * .375 + permuted_coords[0][i] * .25;
          }
	(*this->tag_assigner)( t2, t3, t0, permuted_tags[13] );
        }
      if ( ( comparison_bits & 48 ) == 48 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[10][i] = ( permuted_coords[1][i] + permuted_coords[2][i] ) * .375 + permuted_coords[0][i] * .25;
          }
	(*this->tag_assigner)( t1, t2, t0, permuted_tags[10] );
        }
      MB_TESSELLATOR_INCR_CASE_COUNT(3);
      output_tets.push( MBSimplexTemplateRefiner::templates + 57 );
      output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
      output_sign.push( 1 );
      MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,0);
      switch ( comparison_bits )
        {
        case 42: // 0>2>3<0
          output_tets.push( MBSimplexTemplateRefiner::templates + 62 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,1);
          break;
        case 25: // 2>3>0<2
          output_tets.push( MBSimplexTemplateRefiner::templates + 62 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[11] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,2);
          break;
        case 37: // 3>0>2<3
          output_tets.push( MBSimplexTemplateRefiner::templates + 62 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[3] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,3);
          break;
        case 21: // 3>2>0<3
          output_tets.push( MBSimplexTemplateRefiner::templates + 62 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[22] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,4);
          break;
        case 26: // 2>0>3<2
          output_tets.push( MBSimplexTemplateRefiner::templates + 62 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[12] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,5);
          break;
        case 38: // 0>3>2<0
          output_tets.push( MBSimplexTemplateRefiner::templates + 62 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,6);
          break;
        case 58: // 0=2>3<0
          output_tets.push( MBSimplexTemplateRefiner::templates + 75 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,7);
          break;
        case 29: // 2=3>0<2
          output_tets.push( MBSimplexTemplateRefiner::templates + 75 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[11] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,8);
          break;
        case 39: // 0=3>2<0
          output_tets.push( MBSimplexTemplateRefiner::templates + 75 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[3] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,9);
          break;
        case 53: // 3>0=2<3
          output_tets.push( MBSimplexTemplateRefiner::templates + 96 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,10);
          break;
        case 46: // 0>2=3<0
          output_tets.push( MBSimplexTemplateRefiner::templates + 96 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[11] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,11);
          break;
        case 27: // 2>0=3<2
          output_tets.push( MBSimplexTemplateRefiner::templates + 96 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[3] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,12);
          break;
        case 63: // 0=2=3=0
          output_tets.push( MBSimplexTemplateRefiner::templates + 117 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,13);
          break;
        }
      break;
    case 5: // Ruprecht-Müller Case 3b
      MB_TESSELLATOR_INCR_CASE_COUNT(4);
      output_tets.push( MBSimplexTemplateRefiner::templates + 162 );
      output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
      output_sign.push( 1 );
      MB_TESSELLATOR_INCR_SUBCASE_COUNT(4,0);
      break;
    case 6: // Ruprecht-Müller Case 3c
      comparison_bits = 
        ( permlen[0] <= permlen[1] ? 1 : 0 ) | ( permlen[0] >= permlen[1] ? 2 : 0 ) |
        ( permlen[0] <= permlen[3] ? 4 : 0 ) | ( permlen[0] >= permlen[3] ? 8 : 0 ) |
        0;
      if ( ( comparison_bits & 3 ) == 3 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[10][i] = ( permuted_coords[0][i] + permuted_coords[2][i] ) * .375 + permuted_coords[1][i] * .25;
          }
	(*this->tag_assigner)( t0, t2, t1, permuted_tags[10] );
        }
      if ( ( comparison_bits & 12 ) == 12 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[11][i] = ( permuted_coords[1][i] + permuted_coords[3][i] ) * .375 + permuted_coords[0][i] * .25;
          }
	(*this->tag_assigner)( t1, t3, t0, permuted_tags[11] );
        }
      MB_TESSELLATOR_INCR_CASE_COUNT(5);
      switch ( comparison_bits )
        {
        case 10: // 0>1,0>3
          output_tets.push( MBSimplexTemplateRefiner::templates + 179 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,0);
          break;
        case 5: // 1>0,3>0
          output_tets.push( MBSimplexTemplateRefiner::templates + 200 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,1);
          break;
        case 6: // 0>1,3>0
          output_tets.push( MBSimplexTemplateRefiner::templates + 221 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,2);
          break;
        case 9: // 1>0,0>3
          output_tets.push( MBSimplexTemplateRefiner::templates + 242 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,3);
          break;
        case 11: // 0=1,0>3
          output_tets.push( MBSimplexTemplateRefiner::templates + 263 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,4);
          break;
        case 14: // 0=3,0>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 263 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,5);
          break;
        case 7: // 3>0,0=1
          output_tets.push( MBSimplexTemplateRefiner::templates + 292 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,6);
          break;
        case 13: // 1>0,0=3
          output_tets.push( MBSimplexTemplateRefiner::templates + 292 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,7);
          break;
        case 15: // 0=1,0=3
          output_tets.push( MBSimplexTemplateRefiner::templates + 321 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,8);
          break;
        }
      break;
    case 7: // Ruprecht-Müller Case 3d
      comparison_bits = 
        ( permlen[0] <= permlen[2] ? 1 : 0 ) | ( permlen[0] >= permlen[2] ? 2 : 0 ) |
        ( permlen[0] <= permlen[4] ? 4 : 0 ) | ( permlen[0] >= permlen[4] ? 8 : 0 ) |
        0;
      if ( ( comparison_bits & 3 ) == 3 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[10][i] = ( permuted_coords[1][i] + permuted_coords[2][i] ) * .375 + permuted_coords[0][i] * .25;
          }
	(*this->tag_assigner)( t1, t2, t0, permuted_tags[10] );
        }
      if ( ( comparison_bits & 12 ) == 12 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[11][i] = ( permuted_coords[0][i] + permuted_coords[3][i] ) * .375 + permuted_coords[1][i] * .25;
          }
	(*this->tag_assigner)( t0, t3, t1, permuted_tags[11] );
        }
      MB_TESSELLATOR_INCR_CASE_COUNT(6);
      switch ( comparison_bits )
        {
        case 10: // 0>4,0>2
          output_tets.push( MBSimplexTemplateRefiner::templates + 362 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,0);
          break;
        case 5: // 4>0,2>0
          output_tets.push( MBSimplexTemplateRefiner::templates + 383 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,1);
          break;
        case 9: // 0>4,2>0
          output_tets.push( MBSimplexTemplateRefiner::templates + 404 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,2);
          break;
        case 6: // 4>0,0>2
          output_tets.push( MBSimplexTemplateRefiner::templates + 425 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,3);
          break;
        case 14: // 0=4,0>2
          output_tets.push( MBSimplexTemplateRefiner::templates + 446 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,4);
          break;
        case 11: // 0=2,0>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 446 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,5);
          break;
        case 13: // 2>0,0=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 475 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,6);
          break;
        case 7: // 4>0,0=2
          output_tets.push( MBSimplexTemplateRefiner::templates + 475 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,7);
          break;
        case 15: // 0=4,0=2
          output_tets.push( MBSimplexTemplateRefiner::templates + 504 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,8);
          break;
        }
      break;
    case 8: // Ruprecht-Müller Case 4a
      comparison_bits = 
        ( permlen[4] <= permlen[5] ? 1 : 0 ) | ( permlen[4] >= permlen[5] ? 2 : 0 ) |
        ( permlen[3] <= permlen[4] ? 4 : 0 ) | ( permlen[3] >= permlen[4] ? 8 : 0 ) |
        0;
      if ( ( comparison_bits & 3 ) == 3 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[12][i] = ( permuted_coords[1][i] + permuted_coords[2][i] ) * .375 + permuted_coords[3][i] * .25;
          }
	(*this->tag_assigner)( t1, t2, t3, permuted_tags[12] );
        }
      if ( ( comparison_bits & 12 ) == 12 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[11][i] = ( permuted_coords[0][i] + permuted_coords[1][i] ) * .375 + permuted_coords[3][i] * .25;
          }
	(*this->tag_assigner)( t0, t1, t3, permuted_tags[11] );
        }
      MB_TESSELLATOR_INCR_CASE_COUNT(7);
      output_tets.push( MBSimplexTemplateRefiner::templates + 545 );
      output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
      output_sign.push( 1 );
      MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,0);
      switch ( comparison_bits )
        {
        case 5: // 5>4>3
          output_tets.push( MBSimplexTemplateRefiner::templates + 554 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,1);
          break;
        case 10: // 3>4>5
          output_tets.push( MBSimplexTemplateRefiner::templates + 554 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[13] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,2);
          break;
        case 6: // 3<4>5
          output_tets.push( MBSimplexTemplateRefiner::templates + 571 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,3);
          break;
        case 9: // 3>4<5
          output_tets.push( MBSimplexTemplateRefiner::templates + 588 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,4);
          break;
        case 14: // 3=4>5
          output_tets.push( MBSimplexTemplateRefiner::templates + 605 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,5);
          break;
        case 7: // 4=5,4>3
          output_tets.push( MBSimplexTemplateRefiner::templates + 605 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[13] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,6);
          break;
        case 13: // 5>4,3=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 630 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,7);
          break;
        case 11: // 3>4=5
          output_tets.push( MBSimplexTemplateRefiner::templates + 630 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[13] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,8);
          break;
        case 15: // 3=4=5
          output_tets.push( MBSimplexTemplateRefiner::templates + 655 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,9);
          break;
        }
      break;
    case 9: // Ruprecht-Müller Case 4b
      comparison_bits = 
        ( permlen[1] <= permlen[2] ? 1 : 0 ) | ( permlen[1] >= permlen[2] ? 2 : 0 ) |
        ( permlen[2] <= permlen[3] ? 4 : 0 ) | ( permlen[2] >= permlen[3] ? 8 : 0 ) |
        ( permlen[3] <= permlen[4] ? 16 : 0 ) | ( permlen[3] >= permlen[4] ? 32 : 0 ) |
        ( permlen[1] <= permlen[4] ? 64 : 0 ) | ( permlen[1] >= permlen[4] ? 128 : 0 ) |
        0;
      if ( ( comparison_bits & 3 ) == 3 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[10][i] = ( permuted_coords[1][i] + permuted_coords[0][i] ) * .375 + permuted_coords[2][i] * .25;
          }
	(*this->tag_assigner)( t1, t0, t2, permuted_tags[10] );
        }
      if ( ( comparison_bits & 12 ) == 12 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[13][i] = ( permuted_coords[2][i] + permuted_coords[3][i] ) * .375 + permuted_coords[0][i] * .25;
          }
	(*this->tag_assigner)( t2, t3, t0, permuted_tags[13] );
        }
      if ( ( comparison_bits & 48 ) == 48 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[11][i] = ( permuted_coords[0][i] + permuted_coords[1][i] ) * .375 + permuted_coords[3][i] * .25;
          }
	(*this->tag_assigner)( t0, t1, t3, permuted_tags[11] );
        }
      if ( ( comparison_bits & 192 ) == 192 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[12][i] = ( permuted_coords[2][i] + permuted_coords[3][i] ) * .375 + permuted_coords[1][i] * .25;
          }
	(*this->tag_assigner)( t2, t3, t1, permuted_tags[12] );
        }
      MB_TESSELLATOR_INCR_CASE_COUNT(8);
      switch ( comparison_bits )
        {
        case 85: // 2>1,3>2,4>3,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 688 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,0);
          break;
        case 102: // 1>2,3>2,3>4,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 688 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,1);
          break;
        case 170: // 1>2,2>3,3>4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 688 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,2);
          break;
        case 153: // 2>1,2>3,4>3,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 688 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,3);
          break;
        case 90: // 1>2,2>3,4>3,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 688 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[9] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,4);
          break;
        case 105: // 2>1,2>3,3>4,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 688 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[7] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,5);
          break;
        case 165: // 2>1,3>2,3>4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 688 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[19] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,6);
          break;
        case 150: // 1>2,3>2,4>3,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 688 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[23] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,7);
          break;
        case 101: // 2>1,3>2,3>4,4>1
          {
          int alternates[] = { 713, 738, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 0, 1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,8);
          break;
        case 86: // 1>2,3>2,4>3,4>1
          {
          int alternates[] = {713, 738, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 14, -1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,9);
          break;
        case 154: // 1>2,2>3,4>3,1>4
          {
          int alternates[] = {713, 738, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 5, 1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,10);
          break;
        case 169: // 2>1,2>3,3>4,1>4
          {
          int alternates[] = {713, 738, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 15, -1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,11);
          break;
        case 89: // 2>1,2>3,4>3,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 763 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,12);
          break;
        case 166: // 1>2,3>2,3>4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 763 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,13);
          break;
        case 103: // 1=2,3>2,3>4,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 788 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,14);
          break;
        case 87: // 1=2,3>2,4>3,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 788 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,15);
          break;
        case 185: // 2>1,2>3,3=4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 788 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,16);
          break;
        case 186: // 1>2,2>3,3=4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 788 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,17);
          break;
        case 158: // 1>2,2=3,4>3,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 788 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[9] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,18);
          break;
        case 229: // 2>1,3>2,3>4,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 788 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[7] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,19);
          break;
        case 233: // 2>1,2>3,3>4,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 788 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[19] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,20);
          break;
        case 94: // 1>2,2=3,4>3,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 788 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[23] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,21);
          break;
        case 155: // 1=2,2>3,4>3,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 825 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,22);
          break;
        case 171: // 1=2,2>3,3>4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 825 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,23);
          break;
        case 118: // 1>2,3>2,3=4,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 825 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,24);
          break;
        case 117: // 2>1,3>2,3=4,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 825 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,25);
          break;
        case 109: // 2>1,2=3,3>4,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 825 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[9] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,26);
          break;
        case 218: // 1>2,2>3,4>3,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 825 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[7] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,27);
          break;
        case 214: // 1>2,3>2,4>3,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 825 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[19] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,28);
          break;
        case 173: // 2>1,2=3,3>4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 825 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[23] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,29);
          break;
        case 91: // 1=2,2>3,4>3,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 862 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,30);
          break;
        case 167: // 1=2,3>2,3>4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 862 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,31);
          break;
        case 182: // 1>2,3>2,3=4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 862 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,32);
          break;
        case 121: // 2>1,2>3,3=4,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 862 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,33);
          break;
        case 93: // 2>1,2=3,4>3,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 862 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[9] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,34);
          break;
        case 217: // 2>1,2>3,4>3,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 862 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[7] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,35);
          break;
        case 230: // 1>2,3>2,3>4,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 862 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[19] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,36);
          break;
        case 174: // 1>2,2=3,3>4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 862 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[23] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,37);
          break;
        case 119: // 1=2,3>2,3=4,4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 899 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,38);
          break;
        case 187: // 1=2>3=4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 899 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,39);
          break;
        case 222: // 1>2,2=3,4>3,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 899 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[9] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,40);
          break;
        case 237: // 2>1,2=3,3>4,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 899 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[7] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,41);
          break;
        case 95: // 4>1=2=3,4>3
          output_tets.push( MBSimplexTemplateRefiner::templates + 944 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,42);
          break;
        case 231: // 1=2,3>2,3>4,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 944 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,43);
          break;
        case 190: // 1>2=3=4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 944 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,44);
          break;
        case 249: // 2>1,2>3,3=4,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 944 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,45);
          break;
        case 175: // 1=2=3>4,1>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 993 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,46);
          break;
        case 219: // 1=2>3,4>3,1=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 993 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,47);
          break;
        case 125: // 2>1,2=3=4>1
          output_tets.push( MBSimplexTemplateRefiner::templates + 993 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,48);
          break;
        case 246: // 1>2,3>2,3=4=1
          output_tets.push( MBSimplexTemplateRefiner::templates + 993 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,49);
          break;
        case 255: // 1=2=3=4=1
          output_tets.push( MBSimplexTemplateRefiner::templates + 1042 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,50);
          break;
        }
      break;
    case 10: // Ruprecht-Müller Case 5
      comparison_bits = 
        ( permlen[1] <= permlen[2] ? 1 : 0 ) | ( permlen[1] >= permlen[2] ? 2 : 0 ) |
        ( permlen[3] <= permlen[4] ? 4 : 0 ) | ( permlen[3] >= permlen[4] ? 8 : 0 ) |
        0;
      if ( ( comparison_bits & 3 ) == 3 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[10][i] = ( permuted_coords[1][i] + permuted_coords[0][i] ) * .375 + permuted_coords[2][i] * .25;
          }
	(*this->tag_assigner)( t1, t0, t2, permuted_tags[10] );
        }
      if ( ( comparison_bits & 12 ) == 12 )
        {
        // Compute face point and tag
        for ( int i = 0; i < 6; ++ i )
          {
          permuted_coords[11][i] = ( permuted_coords[0][i] + permuted_coords[1][i] ) * .375 + permuted_coords[3][i] * .25;
          }
	(*this->tag_assigner)( t0, t1, t3, permuted_tags[11] );
        }
      MB_TESSELLATOR_INCR_CASE_COUNT(9);
      output_tets.push( MBSimplexTemplateRefiner::templates + 1107 );
      output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
      output_sign.push( 1 );
      MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,0);
      switch ( comparison_bits )
        {
        case 10: // 1>2,3>4
          output_tets.push( MBSimplexTemplateRefiner::templates + 1116 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,1);
          break;
        case 5: // 2>1,4>3
          output_tets.push( MBSimplexTemplateRefiner::templates + 1116 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,2);
          break;
        case 6: // 1>2,4>3
          {
          int alternates[] = { 1137, 1158, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 0, 1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,3);
          break;
        case 9: // 2>1,3>4
          {
          int alternates[] = {1137, 1158, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 14, -1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,4);
          break;
        case 11: // 1=2,3>4
          {
          int alternates[] = { 1179, 1212, 1245, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 0, 1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,5);
          break;
        case 7: // 1=2,4>3
          {
          int alternates[] = {1179, 1212, 1245, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 14, -1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[14] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,6);
          break;
        case 14: // 3=4,1>2
          {
          int alternates[] = {1179, 1212, 1245, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 5, 1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[5] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,7);
          break;
        case 13: // 3=4,2>1
          {
          int alternates[] = {1179, 1212, 1245, -1 };
          output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 15, -1 ) );
          }
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[15] );
          output_sign.push( -1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,8);
          break;
        case 15: // 1=2,3=4
          output_tets.push( MBSimplexTemplateRefiner::templates + 1278 );
          output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
          output_sign.push( 1 );
          MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,9);
          break;
        }
      break;
    case 11: // Ruprecht-Müller Case 6
      MB_TESSELLATOR_INCR_CASE_COUNT(10);
      output_tets.push( MBSimplexTemplateRefiner::templates + 1319 );
      output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
      output_sign.push( 1 );
      MB_TESSELLATOR_INCR_SUBCASE_COUNT(10,0);
        {
        int alternates[] = { 1336, 1353, 1370, -1 };
        output_tets.push( MBSimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 0, 1 ) );
        }
        output_perm.push( MBSimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(10,1);
        break;
    }

  int* tets;
  int  ntets;
  int* perm;
  int  sgn;
#ifdef MB_DEBUG_TESSELLATOR
  if ( output_tets.empty() )
    {
    cout << "Argh! Case " << C << " Perm " << P << " has no output!" << endl;
    }
#endif // MB_DEBUG_TESSELLATOR
  while ( ! output_tets.empty() )
    {
    tets = output_tets.top();
    ntets = *tets;
    tets++;
    perm = output_perm.top();
    sgn = output_sign.top();

    output_tets.pop();
    output_perm.pop();
    output_sign.pop();

    int t;
    if ( sgn > 0 )
      {
      for ( t = 0; t < ntets; ++t )
        {
        this->refine_3_simplex( max_depth,
                                permuted_coords[perm[tets[0]]], 
                                permuted_tags[perm[tets[0]]], 
                                permuted_coords[perm[tets[1]]],
                                permuted_tags[perm[tets[1]]], 
                                permuted_coords[perm[tets[2]]], 
                                permuted_tags[perm[tets[2]]], 
                                permuted_coords[perm[tets[3]]],
                                permuted_tags[perm[tets[3]]] );
        tets += 4;
        }
      }
    else
      {
      // we have an inverted tet... reverse the first 2 vertices
      // so the orientation is positive.
      for ( t = 0; t < ntets; ++t )
        {
        this->refine_3_simplex( max_depth,
                                permuted_coords[perm[tets[1]]], 
                                permuted_tags[perm[tets[1]]], 
                                permuted_coords[perm[tets[0]]],
                                permuted_tags[perm[tets[0]]], 
                                permuted_coords[perm[tets[2]]], 
                                permuted_tags[perm[tets[2]]], 
                                permuted_coords[perm[tets[3]]],
                                permuted_tags[perm[tets[3]]] );
        tets += 4;
        }
      }
    }

  return true;
}

/**\brief This is used by refine_entity to assign parametric coordinates to corners of each element.
  */
void MBSimplexTemplateRefiner::assign_parametric_coordinates( int num_nodes, const double* src, double* tgt )
{
  for ( int i = 0; i < num_nodes; ++i, src +=3, tgt += 6 )
    for ( int j = 0; j < 3; ++j )
      tgt[j] = src[j];
}

/**\brief Returns true if || a0a1 || < || b0b1 ||
  *
  * We use this to test which triangulation has the best
  * aspect ratio when there are 2 to choose from.
  */
bool MBSimplexTemplateRefiner::compare_Hopf_cross_string_dist(
  const double* a0, const double* a1, const double* b0, const double* b1 )
{   
  double sq_mag_a = 0.;
  double sq_mag_b = 0.;
  for ( int i = 0; i < 3; ++ i )
    {
    double tmp;
    tmp = a0[i] - a1[i];
    sq_mag_a += tmp*tmp;
    tmp = b0[i] - b1[i];
    sq_mag_b += tmp*tmp;
    } 
  return sq_mag_a < sq_mag_b;
}

/*
 * The array below is indexed by the edge code for a tetrahedron.
 * Looking up a row with a tet's edge code will return C and P.
 * C is a configuration number and P is a permutation index. 
 *
 * C is based on the case number from Ruprecht and
 * Müller's (1998) paper on adaptive tetrahedra. (The case
 * numbers are shown to the left of the row in the column
 * labeled case. The only difference is that we introduce
 * a case 3d which is part of case 3c in the paper.)
 *
 * P is an index into the permutations_from_index array below,
 * and is used to transform the current tetrahedron into
 * the canonical configuration associated with C.
 *
 * The 6-digit binary number to the left (which is shown in
 * the horribly UNconventional LSB->MSB order) is the edge
 * code for the row. The 6 digits correspond to the 6 edges
 * of the tetrahedron; a '0' implies no subdivision while
 * a '1' implies subdivision should occur. The ordering of
 * the bits is
 *
 * Edge 0-1, Edge 1-2, Edge 2-0, Edge 0-3, Edge 1-3, Edge 2-3,
 *
 * where the numbers are vertices of the tetrahedron 0-1-2-3.
 * Note that Tet 0-1-2-3 must be positive (i.e., the plane
 * specified by Triangle 0-1-2 must have a normal pointing
 * towards vertex 3, and Triangle 0-1-2's normal must be
 * calculated using the cross-product (Edge 0-1) x (Edge 0-2)).
 *
 * ===========
 * References:
 * (Ruprect and Müller, 1998) A Scheme for Edge-based Adaptive
 *   Tetrahedron Subdivision, Mathematical Visualization (eds.
 *   Hege and Polthier), pp. 61--70. Springer-Verlag. 1998.
 */
int MBSimplexTemplateRefiner::template_index[64][2] =
{
  /*      code case      C    P */
  /* 000000  0  0  */ {  0,   0 },
  /* 100000  1  1  */ {  1,   0 },
  /* 010000  2  1  */ {  1,   1 },
  /* 110000  3  2a */ {  2,   0 },
  /* 001000  4  1  */ {  1,   2 },
  /* 101000  5  2a */ {  2,   2 },
  /* 011000  6  2a */ {  2,   1 },
  /* 111000  7  3b */ {  5,  11 },
  /* 000100  8  1  */ {  1,  10 },
  /* 100100  9  2a */ {  2,   5 },
  /* 010100 10  2b */ {  3,   1 },
  /* 110100 11  3c */ {  6,   0 },
  /* 001100 12  2a */ {  2,  10 },
  /* 101100 13  3a */ {  4,   0 },
  /* 011100 14  3d */ {  7,   2 },
  /* 111100 15  4a */ {  8,   6 },
  /* 000010 16  1  */ {  1,   6 },
  /* 100010 17  2a */ {  2,   4 },
  /* 010010 18  2a */ {  2,   8 },
  /* 110010 19  3a */ {  4,   1 },
  /* 001010 20  2b */ {  3,   2 },
  /* 101010 21  3d */ {  7,   0 },
  /* 011010 22  3c */ {  6,   1 },
  /* 111010 23  4a */ {  8,   9 },
  /* 000110 24  2a */ {  2,   3 },
  /* 100110 25  3b */ {  5,   0 },
  /* 010110 26  3d */ {  7,   4 },
  /* 110110 27  4a */ {  8,  11 },
  /* 001110 28  3c */ {  6,  10 },
  /* 101110 29  4a */ {  8,   7 },
  /* 011110 30  4b */ {  9,   0 },
  /* 111110 31  5  */ { 10,   7 },
  /* 000001 32  1  */ {  1,   7 },
  /* 100001 33  2b */ {  3,   0 },
  /* 010001 34  2a */ {  2,   7 },
  /* 110001 35  3d */ {  7,   1 },
  /* 001001 36  2a */ {  2,  11 },
  /* 101001 37  3c */ {  6,   2 },
  /* 011001 38  3a */ {  4,   2 },
  /* 111001 39  4a */ {  8,   3 },
  /* 000101 40  2a */ {  2,   9 },
  /* 100101 41  3d */ {  7,  10 },
  /* 010101 42  3c */ {  6,   7 },
  /* 110101 43  4b */ {  9,   2 },
  /* 001101 44  3b */ {  5,   7 },
  /* 101101 45  4a */ {  8,   8 },
  /* 011101 46  4a */ {  8,   4 },
  /* 111101 47  5  */ { 10,   6 },
  /* 000011 48  2a */ {  2,   6 },
  /* 100011 49  3c */ {  6,   4 },
  /* 010011 50  3b */ {  5,   1 },
  /* 110011 51  4a */ {  8,  10 },
  /* 001011 52  3d */ {  7,   7 },
  /* 101011 53  4b */ {  9,   1 },
  /* 011011 54  4a */ {  8,   5 },
  /* 111011 55  5  */ { 10,  10 },
  /* 000111 56  3a */ {  4,  10 },
  /* 100111 57  4a */ {  8,   1 },
  /* 010111 58  4a */ {  8,   2 },
  /* 110111 59  5  */ { 10,   2 },
  /* 001111 60  4a */ {  8,   0 },
  /* 101111 61  5  */ { 10,   1 },
  /* 011111 62  5  */ { 10,   0 },
  /* 111111 63  6  */ { 11,   0 },
};


/* Does this mean anything? If so, then you are either 
 * superstitious or much more clever than I (or both?).
 */
/* permutation index, P:  0  1  2  3  4  5  6  7  8  9 10 11 */
/* number of references: 12  9  9  3  4  2  5  6  2  3  7  2 */


/*
 * The array below is a list of all the _positive_
 * permutations of Tetrahedron 0-1-2-3. Given a
 * permutation index, it returns a row of 14 values:
 * these are the vertex numbers of the permuted
 * tetrahedron. The first 4 values are the permuted
 * corner indices, the next 6 values are the
 * permuted edge midpoint indices, and the final
 * entries reference mid-face points inserted
 * to maintain a compatible tetrahedralization.
 *
 * There are 24 entries, 6 for each of the 4 faces of
 * the tetrahedron.
 */
int MBSimplexTemplateRefiner::permutations_from_index[24][14] =
{
  /* corners      midpoints          face points   */
  /* POSITIVE ARRANGEMENTS                         */
  { 0, 1, 2, 3,   4, 5, 6, 7, 8, 9,  10, 11, 12, 13 }, /* Face 0-1-2 */
  { 1, 2, 0, 3,   5, 6, 4, 8, 9, 7,  10, 12, 13, 11 },
  { 2, 0, 1, 3,   6, 4, 5, 9, 7, 8,  10, 13, 11, 12 },

  { 0, 3, 1, 2,   7, 8, 4, 6, 9, 5,  11, 13, 12, 10 }, /* Face 0-3-1 */
  { 3, 1, 0, 2,   8, 4, 7, 9, 5, 6,  11, 12, 10, 13 },
  { 1, 0, 3, 2,   4, 7, 8, 5, 6, 9,  11, 10, 13, 12 },

  { 1, 3, 2, 0,   8, 9, 5, 4, 7, 6,  12, 11, 13, 10 }, /* Face 1-3-2 */
  { 3, 2, 1, 0,   9, 5, 8, 7, 6, 4,  12, 13, 10, 11 },
  { 2, 1, 3, 0,   5, 8, 9, 6, 4, 7,  12, 10, 11, 13 },

  { 2, 3, 0, 1,   9, 7, 6, 5, 8, 4,  13, 12, 11, 10 }, /* Face 2-3-0 */
  { 3, 0, 2, 1,   7, 6, 9, 8, 4, 5,  13, 11, 10, 12 },
  { 0, 2, 3, 1,   6, 9, 7, 4, 5, 8,  13, 10, 12, 11 },

  /* NEGATIVE ARRANGEMENTS                         */
  { 0, 2, 1, 3,   6, 5, 4, 7, 9, 8,  10, 13, 12, 11 }, /* Face 0-1-2 */
  { 2, 1, 0, 3,   5, 4, 6, 9, 8, 7,  10, 12, 11, 13 },
  { 1, 0, 2, 3,   4, 6, 5, 8, 7, 9,  10, 11, 13, 12 },

  { 0, 1, 3, 2,   4, 8, 7, 6, 5, 9,  11, 10, 12, 13 }, /* Face 0-3-1 */
  { 1, 3, 0, 2,   8, 7, 4, 5, 9, 6,  11, 12, 13, 10 },
  { 3, 0, 1, 2,   7, 4, 8, 9, 6, 5,  11, 13, 10, 12 },

  { 1, 2, 3, 0,   5, 9, 8, 4, 6, 7,  12, 10, 13, 11 }, /* Face 1-3-2 */
  { 2, 3, 1, 0,   9, 8, 5, 6, 7, 4,  12, 13, 11, 10 },
  { 3, 1, 2, 0,   8, 5, 9, 7, 4, 6,  12, 11, 10, 13 },

  { 2, 0, 3, 1,   6, 7, 9, 5, 4, 8,  13, 10, 11, 12 }, /* Face 2-3-0 */
  { 0, 3, 2, 1,   7, 9, 6, 4, 8, 5,  13, 11, 12, 10 },
  { 3, 2, 0, 1,   9, 6, 7, 8, 5, 4,  13, 12, 10, 11 }
};

/*
 * Below is a list of output tetrahedra. The array is
 * generated by TessellatorGenerator.py
 * which also generates the code that references it.
 * Each set of tetrahedra begins with a single integer
 * that is the number of tetrahedra for that particular
 * case. It is followed by 5 integers for each output
 * tetrahedron; the first four numbers on each row are
 * indices of the output tetrahedron. The final number
 * is a bit vector specifying which edges of the
 * tetrahedron are internal to the parent tetrahedron
 * being decomposed.
 *
 * Multiple lists of output tetrahedra may be
 * combined to create the tessellation of a single
 * input tetrahedron.
 */

int MBSimplexTemplateRefiner::templates[] = 
{
// case 1_0
   2,
   0,  4,  2,  3,
   4,  1,  2,  3,

// case 2a_0
   1,
   3,  4,  5,  1,

// case 2a, 0>1
   2,
   0,  4,  2,  3,
   4,  5,  2,  3,

// case 2a, 0=1
   4,
  10,  3,  0,  4,
  10,  3,  4,  5,
  10,  3,  5,  2,
  10,  3,  2,  0,

// case 2b_0
   4,
   0,  4,  9,  3,
   4,  1,  9,  3,
   0,  4,  2,  9,
   4,  1,  2,  9,

// case 3a_0
   1,
   4,  7,  6,  0,

// case 3a, 0>2>3<0
   3,
   1,  3,  2,  4,
   4,  6,  3,  2,
   4,  6,  7,  3,

// case 3a, 0=2>3<0
   5,
   4,  6,  7,  3,
  10,  1,  2,  3,
  10,  2,  6,  3,
  10,  6,  4,  3,
  10,  4,  1,  3,

// case 3a, 3>0=2<3
   5,
   1,  3,  2,  7,
  10,  1,  2,  7,
  10,  2,  6,  7,
  10,  6,  4,  7,
  10,  4,  1,  7,

// case 3a, 0=2=3=0
  11,
   2,  6, 10, 13,
   3,  7, 13, 11,
   4,  1, 10, 11,
  11,  6, 10,  4,
  11,  6, 13, 10,
  11,  6,  7, 13,
  11,  6,  4,  7,
   2, 10, 11, 13,
   1, 10, 11,  2,
   2, 11,  3, 13,
   3,  2,  1, 11,

// case 3b_0
   4,
   0,  7,  4,  2,
   4,  7,  8,  2,
   4,  8,  1,  2,
   7,  3,  8,  2,

// case 3c, 0>1,0>3
   5,
   4,  2,  7,  5,
   4,  2,  0,  7,
   4,  3,  1,  5,
   4,  3,  5,  7,
   3,  5,  7,  2,

// case 3c, 1>0,3>0
   5,
   0,  5,  2,  7,
   0,  5,  7,  4,
   7,  1,  4,  5,
   7,  1,  5,  3,
   3,  5,  7,  2,

// case 3c, 0>1,3>0
   5,
   4,  2,  7,  5,
   4,  2,  0,  7,
   7,  1,  4,  5,
   7,  1,  5,  3,
   3,  5,  7,  2,

// case 3c, 1>0,0>3
   5,
   0,  5,  2,  7,
   0,  5,  7,  4,
   4,  3,  1,  5,
   4,  3,  5,  7,
   3,  5,  7,  2,

// case 3c, 0=1,0>3
   7,
   4,  1,  5,  3,
  10,  0,  4,  7,
  10,  2,  0,  7,
  10,  7,  4,  3,
  10,  2,  7,  3,
  10,  5,  2,  3,
  10,  4,  5,  3,

// case 3c, 3>0,0=1
   7,
   7,  1,  5,  3,
   7,  5,  2,  3,
  10,  0,  4,  7,
  10,  2,  0,  7,
  10,  5,  2,  7,
  10,  4,  5,  7,
   1,  5,  4,  7,

// case 3c, 0=1,0=3
  10,
   4,  1,  5, 11,
  11,  1,  5,  3,
  10,  0,  4,  7,
  10,  2,  0,  7,
  10,  5,  2,  3,
  10,  2,  7,  3,
  10,  7,  4, 11,
  10,  7, 11,  3,
  10,  4,  5, 11,
  10, 11,  5,  3,

// case 3d, 0>4,0>2
   5,
   4,  3,  6,  0,
   4,  3,  8,  6,
   4,  2,  8,  1,
   4,  2,  6,  8,
   2,  3,  6,  8,

// case 3d, 4>0,2>0
   5,
   8,  0,  6,  4,
   8,  0,  3,  6,
   6,  1,  8,  4,
   6,  1,  2,  8,
   2,  3,  6,  8,

// case 3d, 0>4,2>0
   5,
   4,  3,  6,  0,
   4,  3,  8,  6,
   6,  1,  8,  4,
   6,  1,  2,  8,
   2,  3,  6,  8,

// case 3d, 4>0,0>2
   5,
   8,  0,  6,  4,
   8,  0,  3,  6,
   4,  2,  8,  1,
   4,  2,  6,  8,
   2,  3,  6,  8,

// case 3d, 0=4,0>2
   7,
   4,  1,  2,  8,
  11,  4,  0,  6,
  11,  0,  3,  6,
  11,  2,  4,  6,
  11,  3,  2,  6,
  11,  3,  8,  2,
  11,  8,  4,  2,

// case 3d, 2>0,0=4
   7,
   6,  2,  8,  1,
   6,  8,  2,  3,
  11,  4,  0,  6,
  11,  0,  3,  6,
   8, 11,  3,  6,
   8,  4, 11,  6,
   1,  6,  4,  8,

// case 3d, 0=4,0=2
  10,
   4,  1, 10,  8,
  10,  2,  8,  1,
  11,  4,  0,  6,
  11,  0,  3,  6,
  11,  3,  8,  2,
  11,  3,  2,  6,
  11, 10,  4,  6,
  11, 10,  6,  2,
   8,  4, 11, 10,
  11, 10,  2,  8,

// case 4a_0
   2,
   7,  8,  9,  3,
   7,  9,  8,  6,

// case 4a, 5>4>3
   4,
   8,  0,  6,  1,
   8,  0,  7,  6,
   9,  1,  6,  2,
   9,  1,  8,  6,

// case 4a, 3<4>5
   4,
   8,  0,  6,  1,
   8,  0,  7,  6,
   8,  2,  6,  9,
   8,  2,  1,  6,

// case 4a, 3>4<5
   4,
   6,  9,  8,  1,
   6,  9,  1,  2,
   6,  7,  0,  1,
   6,  7,  1,  8,

// case 4a, 3=4>5
   6,
   6,  7,  0, 11,
   6,  0,  1, 11,
   6,  7, 11,  8,
   6, 11,  1,  8,
   1,  2,  6,  8,
   2,  6,  8,  9,

// case 4a, 5>4,3=4
   6,
   6,  7,  0, 11,
   6,  0,  1, 11,
   6,  7, 11,  8,
   6, 11,  1,  8,
   1,  2,  6,  9,
   1,  6,  8,  9,

// case 4a, 3=4=5
   8,
   6,  7,  0, 11,
   6,  0,  1, 11,
   6,  7, 11,  8,
   6, 11,  1,  8,
   6,  1,  2, 12,
   6,  2,  9, 12,
   6,  9,  8, 12,
   6,  8,  1, 12,

// case 4b, 2>1,3>2,4>3,4>1
   6,
   6,  8,  1,  5,
   6,  8,  0,  1,
   6,  8,  7,  0,
   6,  8,  2,  7,
   7,  8,  2,  3,
   6,  8,  5,  2,

// case 4b, 2>1,3>2,3>4,4>1
   6,
   6,  8,  1,  5,
   6,  8,  7,  1,
   6,  7,  0,  1,
   8,  7,  3,  2,
   6,  8,  5,  2,
   6,  8,  2,  7,

// case 4b, 2>1,3>2,3>4,4>1, a
   6,
   7,  8,  1,  5,
   6,  5,  7,  1,
   6,  7,  0,  1,
   8,  7,  3,  2,
   7,  8,  5,  2,
   6,  5,  2,  7,

// case 4b, 2>1,2>3,4>3,4>1
   6,
   6,  8,  5,  2,
   6,  8,  2,  3,
   6,  8,  3,  7,
   6,  8,  7,  0,
   6,  8,  0,  1,
   6,  8,  1,  5,

// case 4b, 1=2,3>2,3>4,4>1
   9,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1,  7,
  10,  7,  1,  8,
   6,  7, 10,  8,
   6, 10,  5,  8,
   6,  2,  7,  8,
   6,  5,  2,  8,
   7,  8,  2,  3,

// case 4b, 1=2,2>3,4>3,1>4
   9,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1,  8,
  10,  7,  0,  8,
   6,  7, 10,  8,
   6, 10,  5,  8,
   6,  3,  7,  8,
   6,  5,  3,  8,
   6,  5,  2,  3,

// case 4b, 1=2,2>3,4>3,4>1
   9,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1,  8,
  10,  7,  0,  8,
   6,  7, 10,  8,
   6, 10,  5,  8,
   6,  3,  7,  8,
   6,  5,  2,  8,
   6,  2,  3,  8,

// case 4b, 1=2,3>2,3=4,4>1
  11,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1, 11,
  10, 11,  1,  8,
  10,  0, 11,  7,
  10,  7, 11,  8,
   6,  7, 10,  8,
   6, 10,  5,  8,
   6,  2,  7,  8,
   6,  5,  2,  8,
   7,  8,  2,  3,

// case 4b, 4>1=2=3,4>3
  12,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1,  8,
  10,  7,  0,  8,
  13,  6,  2,  5,
  13,  3,  7,  8,
  13,  2,  3,  8,
  13,  2,  8,  5,
   6,  7, 10,  8,
   6, 10,  5,  8,
   6, 13,  7,  8,
   6,  5, 13,  8,

// case 4b, 1=2=3>4,1>4
  12,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1,  7,
  10,  7,  1,  8,
  13,  6,  2,  5,
  13,  3,  7,  8,
  13,  2,  3,  5,
  13,  3,  8,  5,
   6,  7, 10,  8,
   6, 10,  5,  8,
   6, 13,  7,  8,
   6,  5, 13,  8,

// case 4b, 1=2=3=4=1
  16,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1, 11,
  10, 11,  1,  8,
  10,  0, 11,  7,
  10,  7, 11,  8,
  13,  6,  2,  5,
  13,  3,  7,  8,
  13,  2,  3, 12,
  13,  2, 12,  5,
  13, 12,  3,  8,
  13, 12,  5,  8,
   6,  7, 10,  8,
   6, 10,  5,  8,
   6,  5, 13,  8,
   6, 13,  7,  8,

// case 5_0
   2,
   7,  8,  9,  3,
   6,  5,  2,  9,

// case 5, 1>2,3>4
   5,
   5,  7,  1,  8,
   5,  7,  0,  1,
   5,  7,  6,  0,
   5,  7,  9,  6,
   5,  7,  8,  9,

// case 5, 1>2,4>3
   5,
   0,  5,  6,  7,
   0,  5,  7,  8,
   0,  5,  8,  1,
   5,  7,  9,  6,
   5,  7,  8,  9,

// case 5, 1>2,4>3, a
   5,
   0,  5,  6,  8,
   0,  6,  7,  8,
   0,  5,  8,  1,
   5,  8,  9,  6,
   6,  7,  8,  9,

// case 5, 1=2,3>4
   8,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1,  7,
  10,  7,  1,  8,
  10,  8,  5,  9,
  10,  6,  7,  9,
  10,  7,  8,  9,
  10,  5,  6,  9,

// case 5, 1=2,3>4, a
   8,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1,  7,
  10,  7,  1,  8,
   7,  8,  5,  9,
  10,  6,  7,  5,
  10,  7,  8,  5,
   5,  9,  6,  7,

// case 5, 1=2,3>4, b
   8,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1,  7,
  10,  7,  1,  8,
   6,  8,  5,  9,
  10,  6,  7,  8,
  10,  6,  8,  5,
   8,  9,  6,  7,

// case 5, 1=2,3=4
  10,
  10,  6,  0,  7,
  10,  1,  5,  8,
  10,  0,  1, 11,
  10, 11,  1,  8,
  10,  0, 11,  7,
  10,  7, 11,  8,
  10,  8,  5,  9,
  10,  6,  7,  9,
  10,  7,  8,  9,
  10,  5,  6,  9,

// case 6_0
   4,
   7,  8,  9,  3,
   6,  5,  2,  9,
   4,  1,  5,  8,
   0,  4,  6,  7,

// case 6_1
   4,
   6,  4,  5,  8,
   6,  5,  9,  8,
   6,  9,  7,  8,
   6,  7,  4,  8,

// case 6_1, a
   4,
   5,  8,  9,  7,
   5,  9,  6,  7,
   5,  6,  4,  7,
   5,  4,  8,  7,

// case 6_1, b
   4,
   4,  5,  6,  9,
   4,  6,  7,  9,
   4,  7,  8,  9,
   4,  8,  5,  9,

};

