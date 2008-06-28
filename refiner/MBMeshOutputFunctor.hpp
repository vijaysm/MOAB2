#ifndef MB_MESHOUTPUTFUNCTOR_HPP
#define MB_MESHOUTPUTFUNCTOR_HPP

#include "MBTypes.h"
#include "MBEntityRefiner.hpp"

#include <iostream>
#include <map>

template< int _n >
class MBSplitVertexIndex
{
public:
  MBSplitVertexIndex() { }
  MBSplitVertexIndex( const MBEntityHandle* src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src[i]; std::sort( this->handles, this->handles + _n ); }
  MBSplitVertexIndex( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; }
  MBSplitVertexIndex& operator = ( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; return *this; }

  bool operator < ( const MBSplitVertexIndex<_n>& other ) const
    {
    for ( int i = 0; i < _n; ++ i )
      if ( this->handles[i] < other.handles[i] )
        return true;
      else if ( this->handles[i] > other.handles[i] )
        return false;
    return true;
    }

  MBEntityHandle handles[_n];
};

class MBSplitVerticesBase
{
public:
  MBSplitVerticesBase( MBInterface* m )
    {
    this->mesh = m;
    }
  virtual ~MBSplitVerticesBase() { }
  virtual bool find_or_create( const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle ) = 0;
  MBInterface* mesh;
};

template< int _n >
class MBSplitVertices : public std::map<MBSplitVertexIndex<_n>,MBEntityHandle>, public MBSplitVerticesBase
{
public:
  typedef std::map<MBSplitVertexIndex<_n>,MBEntityHandle> MapType;
  typedef typename std::map<MBSplitVertexIndex<_n>,MBEntityHandle>::iterator MapIteratorType;

  MBSplitVertices( MBInterface* m )
    : MBSplitVerticesBase( m )
    {
    }
  virtual ~MBSplitVertices() { }
  virtual bool find_or_create( const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle )
    {
    MapIteratorType it = this->find( MBSplitVertexIndex<_n>( split_src ) );
    if ( it == this->end() )
      {
      if ( this->mesh->create_vertex( coords, vert_handle ) != MB_SUCCESS )
        {
        return false;
        }
      return true;
      }
    vert_handle = it->second;
    return false;
    }
};

class MBMeshOutputFunctor : public MBEntityRefinerOutputFunctor
{
public:
  MBInterface* mesh;
  bool input_is_output;
  std::vector<MBSplitVerticesBase*> split_vertices;
  std::vector<MBEntityHandle> elem_vert;
  MBRefinerTagManager* tag_manager;
  MBEntityHandle destination_set;

  MBMeshOutputFunctor( MBRefinerTagManager* tag_mgr )
    {
    this->mesh = tag_mgr->get_output_mesh();
    this->input_is_output = ( tag_mgr->get_input_mesh() == this->mesh );
    this->tag_manager = tag_mgr;
    this->destination_set = 0; // don't place output entities in a set by default.

    this->split_vertices.resize( 4 );
    this->split_vertices[0] = 0; // Vertices (0-faces) cannot be split
    this->split_vertices[1] = new MBSplitVertices<1>( this->mesh );
    this->split_vertices[2] = new MBSplitVertices<2>( this->mesh );
    this->split_vertices[3] = new MBSplitVertices<3>( this->mesh );
    }

  ~MBMeshOutputFunctor()
    {
    for ( int i = 0; i < 4; ++ i )
      delete this->split_vertices[i];
    }

  void print_vert_crud( MBEntityHandle vout, int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags )
    {
    std::cout << "+ {";
    for ( int i = 0; i < nvhash; ++ i )
      std::cout << " " << vhash[i];
    std::cout << " } -> " << vout << " ";

    std::cout << "[ " << vcoords[0];
    for ( int i = 1; i < 6; ++i )
      std::cout << ", " << vcoords[i];
    std::cout << " ] ";

#if 0
    double* x = (double*)vtags;
    int* m = (int*)( (char*)vtags + 2 * sizeof( double ) );
    std::cout << "< " << x[0]
              << ", " << x[1];
    for ( int i = 0; i < 4; ++i )
      std::cout << ", " << m[i];
#endif // 0
    std::cout << " >\n";
    }

  void assign_tags( MBEntityHandle vhandle, const void* vtags )
    {
    if ( ! vhandle )
      return; // Ignore bad vertices

    int num_tags = this->tag_manager->get_number_of_vertex_tags();
    MBTag tag_handle;
    int tag_offset;
    for ( int i = 0; i < num_tags; ++i )
      {
      this->tag_manager->get_output_vertex_tag( i, tag_handle, tag_offset );
      this->mesh->tag_set_data( tag_handle, &vhandle, 1, vtags );
      }
    }

  virtual MBEntityHandle operator () ( MBEntityHandle vhash, const double* vcoords, const void* vtags )
    {
    if ( this->input_is_output )
      { // Don't copy the original vertex!
      this->print_vert_crud( vhash, 1, &vhash, vcoords, vtags );
      return vhash;
      }
    MBEntityHandle vertex_handle;
    if ( this->mesh->create_vertex( vcoords + 3, vertex_handle ) != MB_SUCCESS )
      {
      std::cerr << "Could not insert mid-edge vertex!\n";
      }
    this->assign_tags( vertex_handle, vtags );
    this->print_vert_crud( vertex_handle, 1, &vhash, vcoords, vtags );
    return vertex_handle;
    }

  virtual MBEntityHandle operator () ( int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags )
    {
    MBEntityHandle vertex_handle;
    if ( nvhash == 1 )
      {
      vertex_handle = (*this)( *vhash, vcoords, vtags );
      }
    else if ( nvhash < 4 )
      {
      bool newly_created = this->split_vertices[nvhash]->find_or_create( vhash, vcoords, vertex_handle );
      if ( newly_created )
        {
        this->assign_tags( vertex_handle, vtags );
        }
      if ( ! vertex_handle )
        {
        std::cerr << "Could not insert mid-edge vertex!\n";
        }
      this->print_vert_crud( vertex_handle, nvhash, vhash, vcoords, vtags );
      }
    else
      {
      vertex_handle = 0;
      std::cerr << "Not handling splits on faces with " << nvhash << " corners yet.\n";
      }
    return vertex_handle;
    }

  virtual void operator () ( MBEntityHandle h )
    {
    std::cout << h << " ";
    this->elem_vert.push_back( h );
    }

  virtual void operator () ( MBEntityType etyp )
    {
    MBEntityHandle elem_handle;
    if ( this->mesh->create_element( etyp, &this->elem_vert[0], this->elem_vert.size(), elem_handle ) == MB_FAILURE )
      {
      std::cerr << " *** ";
      }
    this->elem_vert.clear();
    std::cout << "---------> " << elem_handle << " ( " << etyp << " )\n\n";
    }
};

#endif // MB_MESHOUTPUTFUNCTOR_HPP
