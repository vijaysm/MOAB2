/*  Filename   :     MBUnkonwnInterface.h
 *  Creator    :     Clinton Stimpson
 *
 *  Date       :     10 Jan 2002
 *
 *  Owner      :     Clinton Stimpson
 *
 *  Description:     Contains declarations for MBuuid which keeps
 *                   track of different interfaces.
 *                   Also contains the declaration for the base class
 *                   MBUknownInterface from which all interfaces are
 *                   derived from
 */

#ifndef MBUNKNOWNINTERFACE_HPP
#define MBUNKNOWNINTERFACE_HPP

#include <memory.h>

typedef unsigned char   MBuchar;
typedef unsigned short  MBushort;
typedef unsigned        MBuint;

//!  struct that handles universally unique id's for the Mesh Database

// note: this MBuuid is compliant with the windows GUID.  
// It is possible to do a memcpy() to copy the data from a MBuuid to a GUID
// if we want to support dll registration
struct MBuuid
{
   //! default constructor that initializes to zero
   MBuuid()
   {
      memset( this, 0, sizeof(MBuuid) );
   }
   //! constructor that takes initialization arguments
   MBuuid( MBuint l, MBushort w1, MBushort w2, 
         MBuchar b1, MBuchar b2, MBuchar b3, 
         MBuchar b4, MBuchar b5, MBuchar b6, 
         MBuchar b7, MBuchar b8 )
   {
      data1 = l;
      data2 = w1;
      data3 = w2;
      data4[0] = b1;
      data4[1] = b2;
      data4[2] = b3;
      data4[3] = b4;
      data4[4] = b5;
      data4[5] = b6;
      data4[6] = b7;
      data4[7] = b8;
   }
   //! copy constructor
   MBuuid( const MBuuid& mdbuuid )
   {
      memcpy( this, &mdbuuid, sizeof(MBuuid));
   }
   //! sets this uuid equal to another one
   MBuuid &operator=(const MBuuid& orig)
   {
      memcpy( this, &orig, sizeof(MBuuid));
      return *this;
   }
   //! returns whether two uuid's are equal
   bool operator==(const MBuuid& orig) const
   {
      return !memcmp(this, &orig, sizeof(MBuuid));
   }
   //! returns whether two uuid's are not equal
   bool operator!=(const MBuuid& orig) const
   {
      return!(*this == orig);
   }

   //! uuid data storage
   MBuint   data1;
   MBushort data2;
   MBushort data3;
   MBuchar  data4[8];
};
  
//! MBuuid for an unknown interface
//! this can be used to either return a default interface
//! or a NULL interface
static const MBuuid IDD_MBUnknown = MBuuid( 0xf4f6605e, 0x2a7e, 0x4760, 
   0xbb, 0x06, 0xb9, 0xed, 0x27, 0xe9, 0x4a, 0xec );


//! base class for all interface classes
class MBUnknownInterface
{
public:
   virtual int QueryInterface
      ( const MBuuid&, MBUnknownInterface** ) = 0;
   virtual ~MBUnknownInterface() {};
};


#endif  // MBUNKNOWNINTERFACE_HPP

