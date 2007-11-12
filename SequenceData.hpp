#ifndef SEQUENCE_DATA_HPP
#define SEQUENCE_DATA_HPP


#include "TypeSequenceManager.hpp"

#include <vector>

class TagServer;

class SequenceData
{
public:

  typedef std::vector<MBEntityHandle>* AdjacencyDataType;

  inline SequenceData( int num_sequence_arrays, 
                       MBEntityHandle start,
                       MBEntityHandle end );
  
  virtual ~SequenceData();
  
  MBEntityHandle start_handle() const 
    { return startHandle; }
  
  MBEntityHandle end_handle() const
    { return endHandle; }
    
  MBEntityID size() const
    { return endHandle + 1 - startHandle; }
    
  void*       get_sequence_data( int array_num )       
                { return arraySet[-1-array_num]; }
  void const* get_sequence_data( int array_num ) const 
                { return arraySet[-1-array_num]; }
  
  AdjacencyDataType*       get_adjacency_data( )       
                { return reinterpret_cast<AdjacencyDataType*>(arraySet[0]); }
  AdjacencyDataType const* get_adjacency_data( ) const 
                { return reinterpret_cast<AdjacencyDataType const*>(arraySet[0]); }
  
  void*       get_tag_data( int tag_num )              
                { return tag_num < numTagData  ? arraySet[tag_num+1] : 0; }
  void const* get_tag_data( int tag_num ) const        
                { return tag_num < numTagData  ? arraySet[tag_num+1] : 0; }
  
  void* create_sequence_data( int array_num, 
                              int bytes_per_ent,
                              void* initial_val = 0 );
                             
  void* create_custom_data( int array_num, size_t total_bytes );
  
  AdjacencyDataType* allocate_adjacency_data();
  
  void* create_tag_data( int tag_num, int bytes_per_ent, void* initial_val = 0 );
  
  SequenceData* subset( MBEntityHandle start, 
                        MBEntityHandle end,
                        const int* sequence_data_sizes,
                        const int* tag_data_sizes ) const;
  
  TypeSequenceManager::SequenceDataPtr seqManData;
  
  void move_tag_data( SequenceData* destination, TagServer* tag_server ) {}
  
protected:

  SequenceData( const SequenceData* subset_from,
                MBEntityHandle start, 
                MBEntityHandle end,
                const int* sequence_data_sizes,
                const int* tag_data_sizes );

private:

  void* create_data( int index, int bytes_per_ent, void* initial_val = 0 );
  void copy_data_subset( int index, 
                         int size_per_ent, 
                         const void* source, 
                         size_t offset, 
                         size_t count );

  const int numSequenceData;
  int numTagData;
  void** arraySet;
  MBEntityHandle startHandle, endHandle;
};

inline SequenceData::SequenceData( int num_sequence_arrays, 
                                   MBEntityHandle start,
                                   MBEntityHandle end )
  : numSequenceData(num_sequence_arrays),
    numTagData(0),
    startHandle(start),
    endHandle(end)
{
  const size_t size = sizeof(void*) * (num_sequence_arrays + 1);
  void** data = (void**)malloc( size );
  memset( data, 0, size );
  arraySet = data + num_sequence_arrays;
}


#endif
