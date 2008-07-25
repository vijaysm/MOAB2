#include "MBProcessSet.hpp"

#include <assert.h>

MBProcessSet::MBProcessSet()
{
  this->clear();
}

MBProcessSet::MBProcessSet( const unsigned char* psetbits )
{
  for ( int i = 0; i < SHARED_PROC_BYTES; ++ i )
    this->processes[i] = psetbits[i];
}

MBProcessSet::~MBProcessSet()
{
}

void MBProcessSet::unite( const MBProcessSet& other )
{
  for ( int i = 0; i < SHARED_PROC_BYTES; ++ i )
    {
    this->processes[i] |= other.processes[i];
    }
}

void MBProcessSet::intersect( const MBProcessSet& other )
{
  for ( int i = 0; i < SHARED_PROC_BYTES; ++ i )
    {
    this->processes[i] &= other.processes[i];
    }
}

void MBProcessSet::clear()
{
  memset( this->processes, 0, SHARED_PROC_BYTES );
}

void MBProcessSet::set_process_member( int i )
{
  int byte = i / 8;
  int bitmask = 1 << ( i % 8 );
  this->processes[byte] |= bitmask;
}

void MBProcessSet::set_process_members( const std::vector<int>& procs )
{
  for ( std::vector<int>::const_iterator it = procs.begin(); it != procs.end() && *it != -1; ++ it )
    {
    this->set_process_member( *it );
    }
}

/**\brief Retrieve a vector containing processes in this set.
  *
  * @param [in] rank The rank of the local process. This integer will not be included in the output list.
  * @param [out] procs The vector in which the list of sharing processes is listed.
  * @return   True when \a rank is the owning process and false otherwise.
  */
bool MBProcessSet::get_process_members( int rank, std::vector<int>& procs )
{
  int i = 0;
  assert( rank >= 0 );
  procs.clear();
  bool rank_owner = false;
  for ( int byte = 0; byte < SHARED_PROC_BYTES; ++ byte )
    {
    char val = this->processes[byte];
    for ( int bit = 0; val && ( bit < 8 ); ++ bit, ++ i )
      {
      if ( ! val )
        {
        i += 8 - bit;
        break;
        }
      else if ( val & 0x1 )
        {
        if ( i != rank )
          procs.push_back( i );
        else if ( ! procs.size() )
          rank_owner = true;
        }
      }
    }
  return rank_owner;
}

bool MBProcessSet::is_process_member( int i ) const
{
  int byte = i / 8;
  int bitmask = 1 << ( i % 8 );
  return ( this->processes[byte] & bitmask ) ? true : false;
}

const unsigned char* MBProcessSet::data() const
{
  return this->processes;
}

bool MBProcessSet::operator < ( const MBProcessSet& other ) const
{
  for ( int i = 0; i < SHARED_PROC_BYTES; ++ i )
    {
    if ( this->processes[i] < other.processes[i] )
      return true;
    else if ( this->processes[i] > other.processes[i] )
      return false;
    }
  return false; // equality
}

std::ostream& operator << ( std::ostream& os, const MBProcessSet& pset )
{
  for ( int i = 0; i < MAX_SHARING_PROCS; ++ i )
    {
    os << ( pset.is_process_member( i ) ? "1" : "0" );
    }
  return os;
}


