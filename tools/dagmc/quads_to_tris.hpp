// Takes a DagMC-style meshset of quads and converts it to triangles.
// It is assumed that quads are only in surface meshsets. Meshset
// membership for tris is only preserved for surfaces meshsets of their
// parent quads.
//

#include <iostream>
#include <assert.h>
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBRange.hpp"

MBErrorCode make_tris_from_quad( MBInterface *MBI,
                                 MBEntityHandle quad,  /* input  */
                                 MBEntityHandle &tri0, /* output */
				 MBEntityHandle &tri1  /* output */);

MBErrorCode quads_to_tris( MBInterface *MBI, MBEntityHandle input_meshset );
