#include "MBProcessSet.hpp"

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


