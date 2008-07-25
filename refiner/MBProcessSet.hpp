#ifndef MB_PROCESSSET_HPP
#define MB_PROCESSSET_HPP

#include "MBTypes.h"
#include "MBParallelComm.hpp"

#include <iostream>
#include <vector>

/**\brief Represent a set of processes using a bit vector.
  *
  * This is used by the mesh refiner when determining where to record
  * split vertices so that labeling can be inferred across process
  * boundaries without communicating anything other than the number of
  * entities in a given partition.
  */
class MBProcessSet
{
public:
  enum
    {
    SHARED_PROC_BYTES = (MAX_SHARING_PROCS / 8 + (MAX_SHARING_PROCS % 8 ? 1 : 0))
    };

  MBProcessSet();
  MBProcessSet( const unsigned char* psetbits );
  ~MBProcessSet();

  void unite( const MBProcessSet& other );
  void intersect( const MBProcessSet& other );

  void clear();

  void set_process_member( int i );
  void set_process_members( const std::vector<int>& procs );

  bool get_process_members( int rank, std::vector<int>& procs );
  bool is_process_member( int i ) const;

  const unsigned char* data() const;

  bool operator < ( const MBProcessSet& other ) const;

  friend std::ostream& operator << ( std::ostream& os, const MBProcessSet& pset );

protected:
  unsigned char processes[SHARED_PROC_BYTES];
};

#endif /* MB_PROCESSSET_HPP */
