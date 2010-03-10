#ifndef VERTEX_SEQUENCE_HPP
#define VERTEX_SEQUENCE_HPP

#include "EntitySequence.hpp"
#include "SequenceData.hpp"

class VertexSequence : public EntitySequence
{
public:

  VertexSequence( MBEntityHandle start,
                  MBEntityID count,
                  SequenceData* data )
    : EntitySequence( start, count, data )
    {}
  
  VertexSequence( MBEntityHandle start,
                  MBEntityID count,
                  MBEntityID data_size )
    : EntitySequence( start, count, new SequenceData( 3, start, start+data_size-1 ) )
    {
      data()->create_sequence_data( X, sizeof(double) );
      data()->create_sequence_data( Y, sizeof(double) );
      data()->create_sequence_data( Z, sizeof(double) );
    } 
  
  virtual ~VertexSequence();
  
  inline MBErrorCode get_coordinates( MBEntityHandle handle,
                                      double& x,
                                      double& y,
                                      double& z ) const;

  inline MBErrorCode get_coordinates( MBEntityHandle handle,
                                      double coords[3] ) const;

  inline MBErrorCode get_coordinates_ref( MBEntityHandle handle,
                                          const double*& x,
                                          const double*& y,
                                          const double*& z ) const;

  inline MBErrorCode set_coordinates( MBEntityHandle entity,
                                      double x, 
                                      double y,
                                      double z );

  inline MBErrorCode set_coordinates( MBEntityHandle entity,
                                      const double xyz[3] );

  inline MBErrorCode get_coordinate_arrays( double*& x, 
                                            double*& y, 
                                            double*& z );

  inline MBErrorCode get_coordinate_arrays( const double*& x,
                                            const double*& y,
                                            const double*& z ) const;
 
  EntitySequence* split( MBEntityHandle here );
  
  SequenceData* create_data_subset( MBEntityHandle start, MBEntityHandle end ) const;
  
  MBErrorCode push_front( MBEntityID count );
  MBErrorCode push_back( MBEntityID count );
  
  void get_const_memory_use( unsigned long& bytes_per_entity,
                             unsigned long& size_of_sequence ) const;
                             
private:

  enum Coord{ X = 0, Y = 1, Z = 2 };

  inline double* array( Coord coord )
  { 
    return reinterpret_cast<double*>(data()->get_sequence_data( coord ));
  }

  inline const double* array( Coord coord ) const
  { 
    return reinterpret_cast<const double*>(data()->get_sequence_data( coord ));
  }
  
  inline double* x_array() { return array(X); }
  inline double* y_array() { return array(Y); }
  inline double* z_array() { return array(Z); }
  
  inline const double* x_array() const { return array(X); }
  inline const double* y_array() const { return array(Y); }
  inline const double* z_array() const { return array(Z); }
  
  VertexSequence( VertexSequence& split_from, MBEntityHandle here )
    : EntitySequence( split_from, here )
    {}
};

  
MBErrorCode VertexSequence::get_coordinates( MBEntityHandle handle,
                                             double& x,
                                             double& y,
                                             double& z ) const
{
  MBEntityID offset = handle - data()->start_handle();
  x = x_array()[offset];
  y = y_array()[offset];
  z = z_array()[offset];
  return MB_SUCCESS;
}

MBErrorCode VertexSequence::get_coordinates( MBEntityHandle handle,
                                             double coords[3] ) const
{
  MBEntityID offset = handle - data()->start_handle();
  coords[X] = x_array()[offset];
  coords[Y] = y_array()[offset];
  coords[Z] = z_array()[offset];
  return MB_SUCCESS;
}
  

MBErrorCode VertexSequence::get_coordinates_ref( MBEntityHandle handle,
                                                 const double*& x,
                                                 const double*& y,
                                                 const double*& z ) const
{
  MBEntityID offset = handle - data()->start_handle();
  x = x_array()+offset;
  y = y_array()+offset;
  z = z_array()+offset;
  return MB_SUCCESS;
}

MBErrorCode VertexSequence::set_coordinates( MBEntityHandle entity,
                                             double x, 
                                             double y,
                                             double z )
{
  MBEntityID offset = entity - data()->start_handle();
  x_array()[offset] = x;
  y_array()[offset] = y;
  z_array()[offset] = z;
  return MB_SUCCESS;
}

MBErrorCode VertexSequence::set_coordinates( MBEntityHandle entity,
                                             const double* xyz )
{
  MBEntityID offset = entity - data()->start_handle();
  x_array()[offset] = xyz[0];
  y_array()[offset] = xyz[1];
  z_array()[offset] = xyz[2];
  return MB_SUCCESS;
}

MBErrorCode VertexSequence::get_coordinate_arrays( double*& x, 
                                                   double*& y, 
                                                   double*& z )
{
  MBEntityID offset = start_handle() - data()->start_handle();
  x = x_array()+offset;
  y = y_array()+offset;
  z = z_array()+offset;
  return MB_SUCCESS;
}
  
MBErrorCode VertexSequence::get_coordinate_arrays( const double*& x,
                                                   const double*& y,
                                                   const double*& z ) const
{
  return get_coordinates_ref( start_handle(), x, y, z ); 
}

#endif
