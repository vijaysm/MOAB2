/*
 * Copyright 2003 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

#include "moab/Refiner.h"
#include "moab/EdgeEvaluator.h"

#include "verdict.h"

namespace moab {

// how's this for avoiding namespace conflicts?! 8-)
static vtkTetra* argyle = 0;
static int argyleRef = 0;
static vtkPoints* goCallTheCops;


#undef UGLY_ASPECT_RATIO_HACK
#undef DBG_MIDPTS

#include <stack>
#include <algorithm>

#ifdef PARAVIEW_DEBUG_TESSELLATOR
#  define VTK_TESSELLATOR_INCR_CASE_COUNT(cs) this->CaseCounts[cs]++
#  define VTK_TESSELLATOR_INCR_SUBCASE_COUNT(cs,sc) this->SubcaseCounts[cs][sc]++
#else // PARAVIEW_DEBUG_TESSELLATOR
#  define VTK_TESSELLATOR_INCR_CASE_COUNT(cs)
#  define VTK_TESSELLATOR_INCR_SUBCASE_COUNT(cs,sc)
#endif // PARAVIEW_DEBUG_TESSELLATOR

Refiner::Refiner()
{
  this->PrivateData = 0;
  this->ConstPrivateData = 0;
  this->Algorithm = 0;
  this->Callback1 = 0;
  this->Callback2 = 0;
  this->Callback3 = 0;
  this->MaximumNumberOfSubdivisions = 3;
  for ( int i=0; i<4; ++i )
    {
    this->EmbeddingDimension[i] = i;
    this->PointDimension[i] = i+3; // By default, FieldSize = 0
    }
  if ( ! argyle )
    {
    argyle = vtkTetra::New();
    argyleRef = 1;
    goCallTheCops = argyle->GetPoints();
    }
  else
    {
    ++argyleRef;
    }
}

Refiner::~Refiner()
{
  if ( this->Algorithm )
    this->Algorithm->UnRegister( this );
  if ( ! (--argyleRef) )
    {
    argyle->Delete();
    argyle = 0;
    }
}

void Refiner::SetEmbeddingDimension( int k, int d )
{
  if ( d > 8 )
    {
    ErrorMacro( "Embedding dimension may not be > 8. (You asked for " << d << "). Truncating to 8." );
    d = 8;
    }

  if ( k == 0 || k < -1 || k >= 4 )
    {
    vtkWarningMacro( "Invalid argument k=" << k );
    return;
    }

  if ( k < 0 )
    {
    for ( k=0; k<4; k++ )
      if ( this->EmbeddingDimension[k] != d )
        {
        this->PointDimension[k] += d - this->EmbeddingDimension[k] ;
        this->EmbeddingDimension[k] = d;
        this->Modified();
        }
    return;
    }
  if ( this->EmbeddingDimension[k] != d )
    {
    this->PointDimension[k] += d - this->EmbeddingDimension[k] ;
    this->EmbeddingDimension[k] = d;
    this->Modified();
    }
}

void Refiner::SetFieldSize( int k, int s )
{
  if ( s > Refiner::MaxFieldSize )
    {
    ErrorMacro( "Embedding dimension may not be > " << MaxFieldSize << ". (You asked for " << s << "). Truncating to " << MaxFieldSize );
    s = Refiner::MaxFieldSize;
    }

  if ( k == 0 || k < -1 || k >= 4 )
    {
    vtkWarningMacro( "Invalid argument k=" << k );
    return;
    }

  if ( k < 0 )
    {
    // Use field size for all facet types (point, line, triangle, tet, ...)
    for ( k=0; k<4; k++ )
      if ( this->PointDimension[k] != s + this->EmbeddingDimension[k] + 3 )
        {
        this->PointDimension[k] = s + this->EmbeddingDimension[k] + 3;
        this->Modified();
        }
    return;
    }

  if ( this->PointDimension[k] != s + this->EmbeddingDimension[k] + 3 )
    {
    this->PointDimension[k] = s + this->EmbeddingDimension[k] + 3;
    this->Modified();
    }
}

void Refiner::SetMaximumNumberOfSubdivisions( int num_subdiv_in )
{
  if ( this->MaximumNumberOfSubdivisions == num_subdiv_in )
    return;

  if ( num_subdiv_in < 0 )
    {
    ErrorMacro( "MaximumNumberOfSubdivisions must be 0 or greater (you requested " << num_subdiv_in << ")" );
    return;
    }

  this->MaximumNumberOfSubdivisions = num_subdiv_in;
  this->Modified();
}

void Refiner::SetTriangleCallback( TriangleProcessorFunction f )
{
  this->Callback2 = f;
}

Refiner::TriangleProcessorFunction Refiner::GetTriangleCallback() const
{
  return this->Callback2;
}

void Refiner::SetTetrahedronCallback( TetrahedronProcessorFunction f )
{
  this->Callback3 = f;
}

Refiner::TetrahedronProcessorFunction Refiner::GetTetrahedronCallback() const
{
  return this->Callback3;
}

void Refiner::SetEdgeCallback( EdgeProcessorFunction f )
{
  this->Callback1 = f;
}

Refiner::EdgeProcessorFunction Refiner::GetEdgeCallback() const
{
  return this->Callback1;
}

void Refiner::SetPrivateData( void* Private )
{
  this->PrivateData = Private;
}

void* Refiner::GetPrivateData() const
{
  return this->PrivateData;
}

void Refiner::SetConstPrivateData( const void* ConstPrivate )
{
  this->ConstPrivateData = ConstPrivate;
}

const void* Refiner::GetConstPrivateData() const
{
  return this->ConstPrivateData;
}

void Refiner::SetSubdivisionAlgorithm( vtkSubdivisionAlgorithm* a )
{
  if ( a != this->Algorithm )
    {
    if ( this->Algorithm )
      this->Algorithm->UnRegister( this );

    this->Algorithm = a;
    this->Modified();

    if ( this->Algorithm )
      this->Algorithm->Register( this );
    }
}

// Returns true if || a0a1 || < || b0b1 ||
// We use this to test which triangulation has the best
// aspect ratio when there are 2 to choose from.
bool compare_Hopf_cross_string_dist( const double* a0, const double* a1, const double* b0, const double* b1 )
{
  double SqMagA = 0.;
  double SqMagB = 0.;
  for (int i=0; i<3; i++)
    {
    double tmp;
    tmp = a0[i] - a1[i];
    SqMagA += tmp*tmp;
    tmp = b0[i] - b1[i];
    SqMagB += tmp*tmp;
    }
  return SqMagA < SqMagB;
}


int Refiner::BestTets( int* connOffsets, double** verts, int permOffset, int sgn ) const
{
  int bestOffset = -1;
  double bestQuality = 0.;
  double currQuality;

  while ( *connOffsets >= 0 )
    {
    int nTets = TetrahedralDecompositions[*connOffsets];
    vtkIdType* conn = &TetrahedralDecompositions[*connOffsets +1];
    int v;
    currQuality = 0.;
    for (v=0; v<nTets; ++v)
      {
      goCallTheCops->SetPoint( 0, verts[ Refiner::PermutationsFromIndex[ permOffset ][ conn[sgn < 0 ? 1:0]] ] );
      goCallTheCops->SetPoint( 1, verts[ Refiner::PermutationsFromIndex[ permOffset ][ conn[sgn < 0 ? 0:1]] ] );
      goCallTheCops->SetPoint( 2, verts[ Refiner::PermutationsFromIndex[ permOffset ][ conn[2]] ] );
      goCallTheCops->SetPoint( 3, verts[ Refiner::PermutationsFromIndex[ permOffset ][ conn[3]] ] );
      currQuality += vtkMeshQuality::TetAspectFrobenius( argyle );
      conn += 4;
      }
    currQuality /= nTets;
    //std::cout << currQuality << " " << *connOffsets << " ";
    if ( bestQuality > currQuality || bestOffset < 0 )
      {
      bestQuality = currQuality;
      bestOffset = *connOffsets;
      }
      ++connOffsets;
    }
    //std::cout << "Choose " << bestOffset << "\n";
  return bestOffset;
}


void Refiner::AdaptivelySample1Facet( double* v0, double* v1, int maxDepth ) const
{
  int edgeCode = 0;

  double midpt0[11+Refiner::MaxFieldSize];
  // make valgrind happy
  std::fill(midpt0,midpt0+this->PointDimension[1],0.);

  if ( maxDepth-- > 0 )
    {
      for ( int i=0; i<this->PointDimension[1]; i++ )
        midpt0[i] = (v0[i] + v1[i])/2.;

      if ( this->Algorithm->EvaluateEdge( v0, midpt0, v1, 3+this->EmbeddingDimension[1] ) )
        edgeCode += 1;
    }

  switch (edgeCode) {
    // No edges to subdivide
    case 0:
      Callback1( v0, v1, this->Algorithm, this->PrivateData, this->ConstPrivateData );
      break ;

      // One edge to subdivide
    case 1:
      this->AdaptivelySample1Facet( v0, midpt0, maxDepth );
      this->AdaptivelySample1Facet( midpt0, v1, maxDepth );
      break;
  }
}

void Refiner::AdaptivelySample2Facet( double* v0, double* v1, double* v2, int maxDepth, int move ) const
{
  int edgeCode = 0;

  double midpt0[11+Refiner::MaxFieldSize];
  double midpt1[11+Refiner::MaxFieldSize];
  double midpt2[11+Refiner::MaxFieldSize];

  // Make valgrind happy
  std::fill(midpt0,midpt0+this->PointDimension[2],0.);
  std::fill(midpt1,midpt1+this->PointDimension[2],0.);
  std::fill(midpt2,midpt2+this->PointDimension[2],0.);

  if ( maxDepth-- > 0 )
    {

    for ( int i=0; i<this->PointDimension[2]; i++ )
      {
      midpt0[i] = (v0[i] + v1[i])/2.;
      midpt1[i] = (v1[i] + v2[i])/2.;
      midpt2[i] = (v2[i] + v0[i])/2.;
      }

    if ( (move & 1) && Algorithm->EvaluateEdge( v0, midpt0, v1, 3+this->EmbeddingDimension[2] ) )
      edgeCode += 1;
    if ( (move & 2) && Algorithm->EvaluateEdge( v1, midpt1, v2, 3+this->EmbeddingDimension[2] ) )
      edgeCode += 2;
    if ( (move & 4) && Algorithm->EvaluateEdge( v2, midpt2, v0, 3+this->EmbeddingDimension[2] ) )
      edgeCode += 4;
#ifdef UGLY_ASPECT_RATIO_HACK
    double dist0=0.;
    double dist1=0.;
    double dist2=0.;
    double tmp;
    for ( int j=0; j<3; ++j )
      {
      tmp = v0[j] - v1[j];
      dist0 += tmp*tmp;
      tmp = v1[j] - v2[j];
      dist1 += tmp*tmp;
      tmp = v2[j] - v0[j];
      dist2 += tmp*tmp;
      }

    if ( edgeCode & 1 ) dist0 /= 2.;
    if ( edgeCode & 2 ) dist1 /= 2.;
    if ( edgeCode & 4 ) dist2 /= 2.;

#define MAR2 2.25
    if ( (!(edgeCode & 1)) && (move&1) && ((dist0/dist1 > MAR2) || (dist0/dist2 > MAR2)) )
      {
      edgeCode += 1;
      move &= 6;
      }
    if ( (!(edgeCode & 2)) && (move&2) && ((dist1/dist0 > MAR2) || (dist1/dist2 > MAR2)) )
      {
      edgeCode += 2;
      move &= 5;
      }
    if ( (!(edgeCode & 4)) && (move&4) && ((dist2/dist1 > MAR2) || (dist2/dist0 > MAR2)) )
      {
      edgeCode += 4;
      move &= 3;
      }
#endif // UGLY_ASPECT_RATIO_HACK
    }

#ifdef DBG_MIDPTS
  if ( maxDepth == 0 )
    {
    fprintf( stderr, "midpoint of v%d (%g %g %g/%g %g %g)-v%d (%g %g %g/%g %g %g) = (%g %g %g/%g %g %g)\n",
      0, v0[0], v0[1], v0[2], v0[3], v0[4], v0[5],
      1, v1[0], v1[1], v1[2], v1[3], v1[4], v1[5],
         midpt0[0], midpt0[1], midpt0[2], midpt0[3], midpt0[4], midpt0[5]
    );

    fprintf( stderr, "midpoint of v%d (%g %g %g/%g %g %g)-v%d (%g %g %g/%g %g %g) = (%g %g %g/%g %g %g)\n",
      1, v1[0], v1[1], v1[2], v1[3], v1[4], v1[5],
      2, v2[0], v2[1], v2[2], v2[3], v2[4], v2[5],
         midpt1[0], midpt1[1], midpt1[2], midpt1[3], midpt1[4], midpt1[5]
    );

    fprintf( stderr, "midpoint of v%d (%g %g %g/%g %g %g)-v%d (%g %g %g/%g %g %g) = (%g %g %g/%g %g %g)\n\n",
      2, v2[0], v2[1], v2[2], v2[3], v2[4], v2[5],
      0, v0[0], v0[1], v0[2], v0[3], v0[4], v0[5],
         midpt2[0], midpt2[1], midpt2[2], midpt2[3], midpt2[4], midpt2[5]
    );
    }
#endif // DBG_MIDPTS

  switch (edgeCode)
    {
    // No edges to subdivide
  case 0:
    Callback2( v0, v1, v2, this->Algorithm, this->PrivateData, this->ConstPrivateData );
    break ;

    // One edge to subdivide
  case 1:
    this->AdaptivelySample2Facet( v0, midpt0, v2, maxDepth, move | 2 );
    this->AdaptivelySample2Facet( midpt0, v1, v2, maxDepth, move | 4 );
    break;
  case 2:
    this->AdaptivelySample2Facet( v0, v1, midpt1, maxDepth, move | 4 );
    this->AdaptivelySample2Facet( v0, midpt1, v2, maxDepth, move | 1 );
    break;
  case 4:
    this->AdaptivelySample2Facet( v0, v1, midpt2, maxDepth, move | 2 );
    this->AdaptivelySample2Facet( midpt2, v1, v2, maxDepth, move | 1 );
    break;

    // Two edges to subdivide
  case 3:
    this->AdaptivelySample2Facet( midpt0, v1, midpt1, maxDepth, move | 4 );
    if ( compare_Hopf_cross_string_dist( v2, midpt0, v0, midpt1 ) )
      {
      this->AdaptivelySample2Facet( midpt0, midpt1,   v2  , maxDepth, move | 5 );
      this->AdaptivelySample2Facet(   v0,   midpt0,   v2  , maxDepth, move | 2 );
      }
    else
      {
      this->AdaptivelySample2Facet(   v0  , midpt0, midpt1, maxDepth, move | 6 );
      this->AdaptivelySample2Facet(   v0,   midpt1,   v2  , maxDepth, move | 1 );
      }
    break;
  case 5:
    this->AdaptivelySample2Facet( v0, midpt0, midpt2, maxDepth, move | 2 );
    if ( compare_Hopf_cross_string_dist( v2, midpt0, v1, midpt2 ) )
      {
      this->AdaptivelySample2Facet( midpt0,   v1,     v2  , maxDepth, move | 4 );
      this->AdaptivelySample2Facet( midpt2, midpt0,   v2  , maxDepth, move | 3 );
      }
    else
      {
      this->AdaptivelySample2Facet( midpt0,   v1,   midpt2, maxDepth, move | 6 );
      this->AdaptivelySample2Facet( midpt2,   v1,     v2,   maxDepth, move | 1 );
      }
    break;
  case 6:
    this->AdaptivelySample2Facet( midpt2, midpt1, v2, maxDepth, move | 1 );
    if ( compare_Hopf_cross_string_dist( v0, midpt1, v1, midpt2 ) )
      {
      this->AdaptivelySample2Facet(   v0,   midpt1, midpt2, maxDepth, move | 3 );
      this->AdaptivelySample2Facet(   v0,     v1,   midpt1, maxDepth, move | 4 );
      }
    else
      {
      this->AdaptivelySample2Facet(   v0,     v1,   midpt2, maxDepth, move | 2 );
      this->AdaptivelySample2Facet( midpt2,   v1,   midpt1, maxDepth, move | 5 );
      }
    break;

    // Three edges to subdivide
  case 7:
    this->AdaptivelySample2Facet( midpt0, midpt1, midpt2, maxDepth, 7 );
    this->AdaptivelySample2Facet(   v0  , midpt0, midpt2, maxDepth, move | 2 );
    this->AdaptivelySample2Facet( midpt0,   v1  , midpt1, maxDepth, move | 4 );
    this->AdaptivelySample2Facet( midpt2, midpt1,   v2  , maxDepth, move | 1 );
    break;
    }
}

void Refiner::AdaptivelySample3Facet( double* v0, double* v1, double* v2, double* v3, int maxDepth ) const
{
  int edgeCode = 0;

  double midpt0[11+Refiner::MaxFieldSize];
  double midpt1[11+Refiner::MaxFieldSize];
  double midpt2[11+Refiner::MaxFieldSize];
  double midpt3[11+Refiner::MaxFieldSize];
  double midpt4[11+Refiner::MaxFieldSize];
  double midpt5[11+Refiner::MaxFieldSize];
#ifdef ALLOW_TET_INTERIOR_PT
  double midpt6[11+Refiner::MaxFieldSize];
#endif // ALLOW_TET_INTERIOR_PT
  double facept0[11+Refiner::MaxFieldSize];
  double facept1[11+Refiner::MaxFieldSize];
  double facept2[11+Refiner::MaxFieldSize];
  double facept3[11+Refiner::MaxFieldSize];

  // Make valgrind happy
  std::fill(midpt0,midpt0+this->PointDimension[3],0.);
  std::fill(midpt1,midpt1+this->PointDimension[3],0.);
  std::fill(midpt2,midpt2+this->PointDimension[3],0.);
  std::fill(midpt3,midpt3+this->PointDimension[3],0.);
  std::fill(midpt4,midpt4+this->PointDimension[3],0.);
  std::fill(midpt5,midpt5+this->PointDimension[3],0.);
#ifdef ALLOW_TET_INTERIOR_PT
  std::fill(midpt6,midpt6+this->PointDimension[3],0.);
#endif // ALLOW_TET_INTERIOR_PT

  double edgeLength2[6];
  if ( maxDepth-- > 0 )
    {
    for ( int i=0; i<this->PointDimension[3]; i++ )
      {
      midpt0[i] = (v0[i] + v1[i])/2.;
      midpt1[i] = (v1[i] + v2[i])/2.;
      midpt2[i] = (v2[i] + v0[i])/2.;
      midpt3[i] = (v0[i] + v3[i])/2.;
      midpt4[i] = (v1[i] + v3[i])/2.;
      midpt5[i] = (v2[i] + v3[i])/2.;
#ifdef ALLOW_TET_INTERIOR_PT
      midpt6[i] = (v1[i]+v2[i]-2*v0[i])/4. + (v3[i]-v0[i])/3. + v0[i];
#endif // ALLOW_TET_INTERIOR_PT
      }

    if ( Algorithm->EvaluateEdge( v0, midpt0, v1, 3+this->EmbeddingDimension[3] ) )
      edgeCode |=  1;
    if ( Algorithm->EvaluateEdge( v1, midpt1, v2, 3+this->EmbeddingDimension[3] ) )
      edgeCode |=  2;
    if ( Algorithm->EvaluateEdge( v2, midpt2, v0, 3+this->EmbeddingDimension[3] ) )
      edgeCode |=  4;

    if ( Algorithm->EvaluateEdge( v0, midpt3, v3, 3+this->EmbeddingDimension[3] ) )
      edgeCode |=  8;
    if ( Algorithm->EvaluateEdge( v1, midpt4, v3, 3+this->EmbeddingDimension[3] ) )
      edgeCode |= 16;
    if ( Algorithm->EvaluateEdge( v2, midpt5, v3, 3+this->EmbeddingDimension[3] ) )
      edgeCode |= 32;

    edgeLength2[0] = edgeLength2[1] = edgeLength2[2] = edgeLength2[3]
      = edgeLength2[4] = edgeLength2[5] = 0;
    for ( int c=0; c<3; ++c )
      {
      double tmp;
      tmp = v1[c] - v0[c];
      edgeLength2[0] += tmp*tmp;
      tmp = v2[c] - v1[c];
      edgeLength2[1] += tmp*tmp;
      tmp = v2[c] - v0[c];
      edgeLength2[2] += tmp*tmp;
      tmp = v3[c] - v0[c];
      edgeLength2[3] += tmp*tmp;
      tmp = v3[c] - v1[c];
      edgeLength2[4] += tmp*tmp;
      tmp = v3[c] - v2[c];
      edgeLength2[5] += tmp*tmp;
      }

#ifdef ALLOW_TET_INTERIOR_PT
    // Find the longest/shortest edges
    double shortest = edgeLength2[0];
    double longest = edgeLength2[0];
    for ( int e=1; e<6; e++ )
      {
      if ( edgeLength2[e] < shortest )
        shortest = edgeLength2[e];
      if ( edgeLength2[e] > longest )
        longest = edgeLength2[e];
      }
    // Divide at center if aspect ratio is > 4:1 (remember these are squares of edge lengths):
    if ( shortest / longest < 1./16. )
      edgeCode |= 64;
#if 0
    // Dunno if we need to have the subdivision algorithm evaluate it, but I suspect
    // it wouldn't be a bad idea.
    if ( Algorithm->EvaluatePoint( midpt6 ) )
      edgeCode |= 64;
#endif // 0
#endif // ALLOW_TET_INTERIOR_PT
    }

  if ( edgeCode == 0 )
    {
    // No edges to subdivide
    Callback3( v0, v1, v2, v3, this->Algorithm, this->PrivateData, this->ConstPrivateData );
    }
  else
    {
    // Do the subdivision
#ifdef ALLOW_TET_INTERIOR_PT
    double* vertices[11] =
    {
      v0, v1, v2, v3,     midpt0, midpt1, midpt2, midpt3, midpt4, midpt5, midpt6
    };
#else // ALLOW_TET_INTERIOR_PT
    double* vertices[10] =
    {
      v0, v1, v2, v3,     midpt0, midpt1, midpt2, midpt3, midpt4, midpt5
    };
#endif // ALLOW_TET_INTERIOR_PT

    // Generate tetrahedra that are compatible except when edge
    // lengths are equal on indeterminately subdivided faces.
    double* permuted[14];
    double permlen[6]; // permuted edge lengths
    int C = Refiner::EdgeCodesToCaseCodesPlusPermutation[ edgeCode ][0];
    int P = Refiner::EdgeCodesToCaseCodesPlusPermutation[ edgeCode ][1];
    int i;

    // 1. Permute the tetrahedron into our canonical configuration
    for ( i=0; i<4; ++i )
      {
      permuted[i] = vertices[ Refiner::PermutationsFromIndex[P][i] ];
      }
    for ( i=4; i<10; ++i )
      {
      // permute the edge lengths, too
      permuted[i] = vertices[ Refiner::PermutationsFromIndex[P][i] ];
      permlen[i-4]  = edgeLength2[ Refiner::PermutationsFromIndex[P][i] - 4 ];
      }
    // Add our local (heap) storage for face points to the list.
    permuted[10] = facept0;
    permuted[11] = facept1;
    permuted[12] = facept2;
    permuted[13] = facept3;

    int comparisonBits;
    std::stack<vtkIdType*> outputTets;
    std::stack<vtkIdType*> outputPerm;
    std::stack<int>        outputSign;

    // cout << "Case " << C << "  Permutation " << P << endl;
    // 2. Generate tetrahedra based on the configuration.
    //    Note that case 0 is handled above (edgeCode == 0).
    switch (C)
      {

    case 1: // Ruprecht-Müller Case 1
      VTK_TESSELLATOR_INCR_CASE_COUNT(0);
      outputTets.push( Refiner::TetrahedralDecompositions + 0 );
      outputPerm.push( Refiner::PermutationsFromIndex[0] );
      outputSign.push( 1 );
      VTK_TESSELLATOR_INCR_SUBCASE_COUNT(0,0);
      break;
    case 2: // Ruprecht-Müller Case 2a
      comparisonBits = 
        (permlen[0] <= permlen[1] ? 1 : 0) | (permlen[0] >= permlen[1] ? 2 : 0) |
        0;
      if ( (comparisonBits & 3) == 3 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[10][i] = (permuted[0][i] + permuted[2][i])*0.375 + permuted[1][i]/4.;
          }
        }
      VTK_TESSELLATOR_INCR_CASE_COUNT(1);
      outputTets.push( Refiner::TetrahedralDecompositions + 9 );
      outputPerm.push( Refiner::PermutationsFromIndex[0] );
      outputSign.push( 1 );
      VTK_TESSELLATOR_INCR_SUBCASE_COUNT(1,0);
      switch (comparisonBits)
        {
        case 2: // 0>1
          outputTets.push( Refiner::TetrahedralDecompositions + 14 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(1,1);
          break;
        case 1: // 1>0
          outputTets.push( Refiner::TetrahedralDecompositions + 14 );
          outputPerm.push( Refiner::PermutationsFromIndex[13] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(1,2);
          break;
        case 3: // 0=1
          outputTets.push( Refiner::TetrahedralDecompositions + 23 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(1,3);
          break;
        }
      break;
    case 3: // Ruprecht-Müller Case 2b
      VTK_TESSELLATOR_INCR_CASE_COUNT(2);
      outputTets.push( Refiner::TetrahedralDecompositions + 40 );
      outputPerm.push( Refiner::PermutationsFromIndex[0] );
      outputSign.push( 1 );
      VTK_TESSELLATOR_INCR_SUBCASE_COUNT(2,0);
      break;
    case 4: // Ruprecht-Müller Case 3a
      comparisonBits = 
        (permlen[0] <= permlen[3] ? 1 : 0) | (permlen[0] >= permlen[3] ? 2 : 0) |
        (permlen[2] <= permlen[3] ? 4 : 0) | (permlen[2] >= permlen[3] ? 8 : 0) |
        (permlen[0] <= permlen[2] ? 16 : 0) | (permlen[0] >= permlen[2] ? 32 : 0) |
        0;
      if ( (comparisonBits & 3) == 3 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[11][i] = (permuted[1][i] + permuted[3][i])*0.375 + permuted[0][i]/4.;
          }
        }
      if ( (comparisonBits & 12) == 12 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[13][i] = (permuted[2][i] + permuted[3][i])*0.375 + permuted[0][i]/4.;
          }
        }
      if ( (comparisonBits & 48) == 48 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[10][i] = (permuted[1][i] + permuted[2][i])*0.375 + permuted[0][i]/4.;
          }
        }
      VTK_TESSELLATOR_INCR_CASE_COUNT(3);
      outputTets.push( Refiner::TetrahedralDecompositions + 57 );
      outputPerm.push( Refiner::PermutationsFromIndex[0] );
      outputSign.push( 1 );
      VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,0);
      switch (comparisonBits)
        {
        case 42: // 0>2>3<0
          outputTets.push( Refiner::TetrahedralDecompositions + 62 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,1);
          break;
        case 25: // 2>3>0<2
          outputTets.push( Refiner::TetrahedralDecompositions + 62 );
          outputPerm.push( Refiner::PermutationsFromIndex[11] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,2);
          break;
        case 37: // 3>0>2<3
          outputTets.push( Refiner::TetrahedralDecompositions + 62 );
          outputPerm.push( Refiner::PermutationsFromIndex[3] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,3);
          break;
        case 21: // 3>2>0<3
          outputTets.push( Refiner::TetrahedralDecompositions + 62 );
          outputPerm.push( Refiner::PermutationsFromIndex[22] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,4);
          break;
        case 26: // 2>0>3<2
          outputTets.push( Refiner::TetrahedralDecompositions + 62 );
          outputPerm.push( Refiner::PermutationsFromIndex[12] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,5);
          break;
        case 38: // 0>3>2<0
          outputTets.push( Refiner::TetrahedralDecompositions + 62 );
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,6);
          break;
        case 58: // 0=2>3<0
          outputTets.push( Refiner::TetrahedralDecompositions + 75 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,7);
          break;
        case 29: // 2=3>0<2
          outputTets.push( Refiner::TetrahedralDecompositions + 75 );
          outputPerm.push( Refiner::PermutationsFromIndex[11] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,8);
          break;
        case 39: // 0=3>2<0
          outputTets.push( Refiner::TetrahedralDecompositions + 75 );
          outputPerm.push( Refiner::PermutationsFromIndex[3] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,9);
          break;
        case 53: // 3>0=2<3
          outputTets.push( Refiner::TetrahedralDecompositions + 96 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,10);
          break;
        case 46: // 0>2=3<0
          outputTets.push( Refiner::TetrahedralDecompositions + 96 );
          outputPerm.push( Refiner::PermutationsFromIndex[11] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,11);
          break;
        case 27: // 2>0=3<2
          outputTets.push( Refiner::TetrahedralDecompositions + 96 );
          outputPerm.push( Refiner::PermutationsFromIndex[3] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,12);
          break;
        case 63: // 0=2=3=0
          outputTets.push( Refiner::TetrahedralDecompositions + 117 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(3,13);
          break;
        }
      break;
    case 5: // Ruprecht-Müller Case 3b
      VTK_TESSELLATOR_INCR_CASE_COUNT(4);
      outputTets.push( Refiner::TetrahedralDecompositions + 162 );
      outputPerm.push( Refiner::PermutationsFromIndex[0] );
      outputSign.push( 1 );
      VTK_TESSELLATOR_INCR_SUBCASE_COUNT(4,0);
      break;
    case 6: // Ruprecht-Müller Case 3c
      comparisonBits = 
        (permlen[0] <= permlen[1] ? 1 : 0) | (permlen[0] >= permlen[1] ? 2 : 0) |
        (permlen[0] <= permlen[3] ? 4 : 0) | (permlen[0] >= permlen[3] ? 8 : 0) |
        0;
      if ( (comparisonBits & 3) == 3 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[10][i] = (permuted[0][i] + permuted[2][i])*0.375 + permuted[1][i]/4.;
          }
        }
      if ( (comparisonBits & 12) == 12 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[11][i] = (permuted[1][i] + permuted[3][i])*0.375 + permuted[0][i]/4.;
          }
        }
      VTK_TESSELLATOR_INCR_CASE_COUNT(5);
      switch (comparisonBits)
        {
        case 10: // 0>1,0>3
          outputTets.push( Refiner::TetrahedralDecompositions + 179 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(5,0);
          break;
        case 5: // 1>0,3>0
          outputTets.push( Refiner::TetrahedralDecompositions + 200 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(5,1);
          break;
        case 6: // 0>1,3>0
          outputTets.push( Refiner::TetrahedralDecompositions + 221 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(5,2);
          break;
        case 9: // 1>0,0>3
          outputTets.push( Refiner::TetrahedralDecompositions + 242 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(5,3);
          break;
        case 11: // 0=1,0>3
          outputTets.push( Refiner::TetrahedralDecompositions + 263 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(5,4);
          break;
        case 14: // 0=3,0>1
          outputTets.push( Refiner::TetrahedralDecompositions + 263 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(5,5);
          break;
        case 7: // 3>0,0=1
          outputTets.push( Refiner::TetrahedralDecompositions + 292 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(5,6);
          break;
        case 13: // 1>0,0=3
          outputTets.push( Refiner::TetrahedralDecompositions + 292 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(5,7);
          break;
        case 15: // 0=1,0=3
          outputTets.push( Refiner::TetrahedralDecompositions + 321 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(5,8);
          break;
        }
      break;
    case 7: // Ruprecht-Müller Case 3d
      comparisonBits = 
        (permlen[0] <= permlen[2] ? 1 : 0) | (permlen[0] >= permlen[2] ? 2 : 0) |
        (permlen[0] <= permlen[4] ? 4 : 0) | (permlen[0] >= permlen[4] ? 8 : 0) |
        0;
      if ( (comparisonBits & 3) == 3 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[10][i] = (permuted[1][i] + permuted[2][i])*0.375 + permuted[0][i]/4.;
          }
        }
      if ( (comparisonBits & 12) == 12 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[11][i] = (permuted[0][i] + permuted[3][i])*0.375 + permuted[1][i]/4.;
          }
        }
      VTK_TESSELLATOR_INCR_CASE_COUNT(6);
      switch (comparisonBits)
        {
        case 10: // 0>4,0>2
          outputTets.push( Refiner::TetrahedralDecompositions + 362 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(6,0);
          break;
        case 5: // 4>0,2>0
          outputTets.push( Refiner::TetrahedralDecompositions + 383 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(6,1);
          break;
        case 9: // 0>4,2>0
          outputTets.push( Refiner::TetrahedralDecompositions + 404 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(6,2);
          break;
        case 6: // 4>0,0>2
          outputTets.push( Refiner::TetrahedralDecompositions + 425 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(6,3);
          break;
        case 14: // 0=4,0>2
          outputTets.push( Refiner::TetrahedralDecompositions + 446 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(6,4);
          break;
        case 11: // 0=2,0>4
          outputTets.push( Refiner::TetrahedralDecompositions + 446 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(6,5);
          break;
        case 13: // 2>0,0=4
          outputTets.push( Refiner::TetrahedralDecompositions + 475 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(6,6);
          break;
        case 7: // 4>0,0=2
          outputTets.push( Refiner::TetrahedralDecompositions + 475 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(6,7);
          break;
        case 15: // 0=4,0=2
          outputTets.push( Refiner::TetrahedralDecompositions + 504 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(6,8);
          break;
        }
      break;
    case 8: // Ruprecht-Müller Case 4a
      comparisonBits = 
        (permlen[4] <= permlen[5] ? 1 : 0) | (permlen[4] >= permlen[5] ? 2 : 0) |
        (permlen[3] <= permlen[4] ? 4 : 0) | (permlen[3] >= permlen[4] ? 8 : 0) |
        0;
      if ( (comparisonBits & 3) == 3 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[12][i] = (permuted[1][i] + permuted[2][i])*0.375 + permuted[3][i]/4.;
          }
        }
      if ( (comparisonBits & 12) == 12 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[11][i] = (permuted[0][i] + permuted[1][i])*0.375 + permuted[3][i]/4.;
          }
        }
      VTK_TESSELLATOR_INCR_CASE_COUNT(7);
      outputTets.push( Refiner::TetrahedralDecompositions + 545 );
      outputPerm.push( Refiner::PermutationsFromIndex[0] );
      outputSign.push( 1 );
      VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,0);
      switch (comparisonBits)
        {
        case 5: // 5>4>3
          outputTets.push( Refiner::TetrahedralDecompositions + 554 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,1);
          break;
        case 10: // 3>4>5
          outputTets.push( Refiner::TetrahedralDecompositions + 554 );
          outputPerm.push( Refiner::PermutationsFromIndex[13] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,2);
          break;
        case 6: // 3<4>5
          outputTets.push( Refiner::TetrahedralDecompositions + 571 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,3);
          break;
        case 9: // 3>4<5
          outputTets.push( Refiner::TetrahedralDecompositions + 588 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,4);
          break;
        case 14: // 3=4>5
          outputTets.push( Refiner::TetrahedralDecompositions + 605 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,5);
          break;
        case 7: // 4=5,4>3
          outputTets.push( Refiner::TetrahedralDecompositions + 605 );
          outputPerm.push( Refiner::PermutationsFromIndex[13] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,6);
          break;
        case 13: // 5>4,3=4
          outputTets.push( Refiner::TetrahedralDecompositions + 630 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,7);
          break;
        case 11: // 3>4=5
          outputTets.push( Refiner::TetrahedralDecompositions + 630 );
          outputPerm.push( Refiner::PermutationsFromIndex[13] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,8);
          break;
        case 15: // 3=4=5
          outputTets.push( Refiner::TetrahedralDecompositions + 655 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(7,9);
          break;
        }
      break;
    case 9: // Ruprecht-Müller Case 4b
      comparisonBits = 
        (permlen[1] <= permlen[2] ? 1 : 0) | (permlen[1] >= permlen[2] ? 2 : 0) |
        (permlen[2] <= permlen[3] ? 4 : 0) | (permlen[2] >= permlen[3] ? 8 : 0) |
        (permlen[3] <= permlen[4] ? 16 : 0) | (permlen[3] >= permlen[4] ? 32 : 0) |
        (permlen[1] <= permlen[4] ? 64 : 0) | (permlen[1] >= permlen[4] ? 128 : 0) |
        0;
      if ( (comparisonBits & 3) == 3 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[10][i] = (permuted[1][i] + permuted[0][i])*0.375 + permuted[2][i]/4.;
          }
        }
      if ( (comparisonBits & 12) == 12 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[13][i] = (permuted[2][i] + permuted[3][i])*0.375 + permuted[0][i]/4.;
          }
        }
      if ( (comparisonBits & 48) == 48 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[11][i] = (permuted[0][i] + permuted[1][i])*0.375 + permuted[3][i]/4.;
          }
        }
      if ( (comparisonBits & 192) == 192 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[12][i] = (permuted[2][i] + permuted[3][i])*0.375 + permuted[1][i]/4.;
          }
        }
      VTK_TESSELLATOR_INCR_CASE_COUNT(8);
      switch (comparisonBits)
        {
        case 85: // 2>1,3>2,4>3,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 688 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,0);
          break;
        case 102: // 1>2,3>2,3>4,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 688 );
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,1);
          break;
        case 170: // 1>2,2>3,3>4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 688 );
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,2);
          break;
        case 153: // 2>1,2>3,4>3,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 688 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,3);
          break;
        case 90: // 1>2,2>3,4>3,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 688 );
          outputPerm.push( Refiner::PermutationsFromIndex[9] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,4);
          break;
        case 105: // 2>1,2>3,3>4,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 688 );
          outputPerm.push( Refiner::PermutationsFromIndex[7] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,5);
          break;
        case 165: // 2>1,3>2,3>4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 688 );
          outputPerm.push( Refiner::PermutationsFromIndex[19] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,6);
          break;
        case 150: // 1>2,3>2,4>3,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 688 );
          outputPerm.push( Refiner::PermutationsFromIndex[23] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,7);
          break;
        case 101: // 2>1,3>2,3>4,4>1
          {
            int alternates[] = { 713, 738, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 0, 1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,8);
          break;
        case 86: // 1>2,3>2,4>3,4>1
          {
            int alternates[] = {713, 738, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 14, -1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,9);
          break;
        case 154: // 1>2,2>3,4>3,1>4
          {
            int alternates[] = {713, 738, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 5, 1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,10);
          break;
        case 169: // 2>1,2>3,3>4,1>4
          {
            int alternates[] = {713, 738, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 15, -1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,11);
          break;
        case 89: // 2>1,2>3,4>3,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 763 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,12);
          break;
        case 166: // 1>2,3>2,3>4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 763 );
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,13);
          break;
        case 103: // 1=2,3>2,3>4,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 788 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,14);
          break;
        case 87: // 1=2,3>2,4>3,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 788 );
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,15);
          break;
        case 185: // 2>1,2>3,3=4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 788 );
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,16);
          break;
        case 186: // 1>2,2>3,3=4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 788 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,17);
          break;
        case 158: // 1>2,2=3,4>3,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 788 );
          outputPerm.push( Refiner::PermutationsFromIndex[9] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,18);
          break;
        case 229: // 2>1,3>2,3>4,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 788 );
          outputPerm.push( Refiner::PermutationsFromIndex[7] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,19);
          break;
        case 233: // 2>1,2>3,3>4,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 788 );
          outputPerm.push( Refiner::PermutationsFromIndex[19] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,20);
          break;
        case 94: // 1>2,2=3,4>3,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 788 );
          outputPerm.push( Refiner::PermutationsFromIndex[23] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,21);
          break;
        case 155: // 1=2,2>3,4>3,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 825 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,22);
          break;
        case 171: // 1=2,2>3,3>4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 825 );
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,23);
          break;
        case 118: // 1>2,3>2,3=4,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 825 );
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,24);
          break;
        case 117: // 2>1,3>2,3=4,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 825 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,25);
          break;
        case 109: // 2>1,2=3,3>4,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 825 );
          outputPerm.push( Refiner::PermutationsFromIndex[9] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,26);
          break;
        case 218: // 1>2,2>3,4>3,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 825 );
          outputPerm.push( Refiner::PermutationsFromIndex[7] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,27);
          break;
        case 214: // 1>2,3>2,4>3,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 825 );
          outputPerm.push( Refiner::PermutationsFromIndex[19] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,28);
          break;
        case 173: // 2>1,2=3,3>4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 825 );
          outputPerm.push( Refiner::PermutationsFromIndex[23] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,29);
          break;
        case 91: // 1=2,2>3,4>3,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 862 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,30);
          break;
        case 167: // 1=2,3>2,3>4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 862 );
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,31);
          break;
        case 182: // 1>2,3>2,3=4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 862 );
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,32);
          break;
        case 121: // 2>1,2>3,3=4,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 862 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,33);
          break;
        case 93: // 2>1,2=3,4>3,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 862 );
          outputPerm.push( Refiner::PermutationsFromIndex[9] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,34);
          break;
        case 217: // 2>1,2>3,4>3,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 862 );
          outputPerm.push( Refiner::PermutationsFromIndex[7] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,35);
          break;
        case 230: // 1>2,3>2,3>4,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 862 );
          outputPerm.push( Refiner::PermutationsFromIndex[19] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,36);
          break;
        case 174: // 1>2,2=3,3>4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 862 );
          outputPerm.push( Refiner::PermutationsFromIndex[23] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,37);
          break;
        case 119: // 1=2,3>2,3=4,4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 899 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,38);
          break;
        case 187: // 1=2>3=4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 899 );
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,39);
          break;
        case 222: // 1>2,2=3,4>3,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 899 );
          outputPerm.push( Refiner::PermutationsFromIndex[9] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,40);
          break;
        case 237: // 2>1,2=3,3>4,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 899 );
          outputPerm.push( Refiner::PermutationsFromIndex[7] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,41);
          break;
        case 95: // 4>1=2=3,4>3
          outputTets.push( Refiner::TetrahedralDecompositions + 944 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,42);
          break;
        case 231: // 1=2,3>2,3>4,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 944 );
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,43);
          break;
        case 190: // 1>2=3=4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 944 );
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,44);
          break;
        case 249: // 2>1,2>3,3=4,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 944 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,45);
          break;
        case 175: // 1=2=3>4,1>4
          outputTets.push( Refiner::TetrahedralDecompositions + 993 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,46);
          break;
        case 219: // 1=2>3,4>3,1=4
          outputTets.push( Refiner::TetrahedralDecompositions + 993 );
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,47);
          break;
        case 125: // 2>1,2=3=4>1
          outputTets.push( Refiner::TetrahedralDecompositions + 993 );
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,48);
          break;
        case 246: // 1>2,3>2,3=4=1
          outputTets.push( Refiner::TetrahedralDecompositions + 993 );
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,49);
          break;
        case 255: // 1=2=3=4=1
          outputTets.push( Refiner::TetrahedralDecompositions + 1042 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(8,50);
          break;
        }
      break;
    case 10: // Ruprecht-Müller Case 5
      comparisonBits = 
        (permlen[1] <= permlen[2] ? 1 : 0) | (permlen[1] >= permlen[2] ? 2 : 0) |
        (permlen[3] <= permlen[4] ? 4 : 0) | (permlen[3] >= permlen[4] ? 8 : 0) |
        0;
      if ( (comparisonBits & 3) == 3 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[10][i] = (permuted[1][i] + permuted[0][i])*0.375 + permuted[2][i]/4.;
          }
        }
      if ( (comparisonBits & 12) == 12 )
        {
        // Compute face point
        for ( i=0; i<this->PointDimension[3]; i++ )
          {
          permuted[11][i] = (permuted[0][i] + permuted[1][i])*0.375 + permuted[3][i]/4.;
          }
        }
      VTK_TESSELLATOR_INCR_CASE_COUNT(9);
      outputTets.push( Refiner::TetrahedralDecompositions + 1107 );
      outputPerm.push( Refiner::PermutationsFromIndex[0] );
      outputSign.push( 1 );
      VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,0);
      switch (comparisonBits)
        {
        case 10: // 1>2,3>4
          outputTets.push( Refiner::TetrahedralDecompositions + 1116 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,1);
          break;
        case 5: // 2>1,4>3
          outputTets.push( Refiner::TetrahedralDecompositions + 1116 );
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,2);
          break;
        case 6: // 1>2,4>3
          {
            int alternates[] = { 1137, 1158, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 0, 1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,3);
          break;
        case 9: // 2>1,3>4
          {
            int alternates[] = {1137, 1158, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 14, -1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,4);
          break;
        case 11: // 1=2,3>4
          {
            int alternates[] = { 1179, 1212, 1245, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 0, 1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,5);
          break;
        case 7: // 1=2,4>3
          {
            int alternates[] = {1179, 1212, 1245, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 14, -1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[14] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,6);
          break;
        case 14: // 3=4,1>2
          {
            int alternates[] = {1179, 1212, 1245, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 5, 1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[5] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,7);
          break;
        case 13: // 3=4,2>1
          {
            int alternates[] = {1179, 1212, 1245, -1 };
            outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 15, -1 ) );
          }
          outputPerm.push( Refiner::PermutationsFromIndex[15] );
          outputSign.push( -1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,8);
          break;
        case 15: // 1=2,3=4
          outputTets.push( Refiner::TetrahedralDecompositions + 1278 );
          outputPerm.push( Refiner::PermutationsFromIndex[0] );
          outputSign.push( 1 );
          VTK_TESSELLATOR_INCR_SUBCASE_COUNT(9,9);
          break;
        }
      break;
    case 11: // Ruprecht-Müller Case 6
      VTK_TESSELLATOR_INCR_CASE_COUNT(10);
      outputTets.push( Refiner::TetrahedralDecompositions + 1319 );
      outputPerm.push( Refiner::PermutationsFromIndex[0] );
      outputSign.push( 1 );
      VTK_TESSELLATOR_INCR_SUBCASE_COUNT(10,0);
        {
        int alternates[] = { 1336, 1353, 1370, -1 };
        outputTets.push( Refiner::TetrahedralDecompositions + this->BestTets( alternates, permuted, 0, 1 ) );
        }
      outputPerm.push( Refiner::PermutationsFromIndex[0] );
      outputSign.push( 1 );
      VTK_TESSELLATOR_INCR_SUBCASE_COUNT(10,1);
      break;
      }

    vtkIdType* tets;
    vtkIdType  ntets;
    vtkIdType* perm;
    int        sgn;
#ifdef PARAVIEW_DEBUG_TESSELLATOR
    if ( outputTets.empty() )
    {
    cout << "Argh! Case " << C << " Perm " << P << " has no output!" << endl;
    }
#endif // PARAVIEW_DEBUG_TESSELLATOR
    while ( ! outputTets.empty() )
      {
      tets = outputTets.top();
      ntets = *tets;
      tets++;
      perm = outputPerm.top();
      sgn = outputSign.top();

      outputTets.pop();
      outputPerm.pop();
      outputSign.pop();

      int t;
      if ( sgn > 0 )
        {
        for ( t = 0; t < ntets; ++t )
          {
          this->AdaptivelySample3Facet(
            permuted[ perm[ tets[0] ] ], permuted[ perm[ tets[1] ] ],
            permuted[ perm[ tets[2] ] ], permuted[ perm[ tets[3] ] ],
            maxDepth );
          tets += 4;
          }
        }
      else
        {
        // we have an inverted tet... reverse the first 2 vertices
        // so the orientation is positive.
        for ( t = 0; t < ntets; ++t )
          {
          this->AdaptivelySample3Facet(
            permuted[ perm[ tets[1] ] ], permuted[ perm[ tets[0] ] ],
            permuted[ perm[ tets[2] ] ], permuted[ perm[ tets[3] ] ],
            maxDepth );
          tets += 4;
          }
        }
      }
    }
}

} // namespace moab
