#ifndef READWRITE_DEFINES_HPP
#define READWRITE_DEFINES_HPP


#ifdef WIN32
#pragma warning(disable : 4786)
#endif

/* CJS  -- we need to figure out how this fits in the component framework */
#error "don't include this file"


static const char* element_type_names[] =
{
  "BAR",
  "BAR2",
  "BAR3",
  "BEAM",
  "BEAM2",
  "BEAM3",
  "TRUSS",
  "TRUSS2",
  "TRUSS3",
  "QUAD",
  "QUAD4",
  "QUAD5",
  "QUAD8",
  "QUAD9",
  "SHELL",
  "SHELL4",
  "SHELL8",
  "SHELL9",
  "TRI",
  "TRI3",
  "TRI6",
  "TRI7",
  "HEX",
  "HEX8",
  "HEX9",
  "HEX20",
  "HEX27",
  "PYRAMID",
  "TETRA",
  "TETRA4",
  "TETRA8",
  "TETRA10",
  "TETRA14",
  "KNIFE",
  "WEDGE"
};


enum MBElementType 
{
  MB_BAR = 0,
  MB_BAR2,
  MB_BAR3,
  MB_BEAM,
  MB_BEAM2,
  MB_BEAM3,
  MB_TRUSS,
  MB_TRUSS2,
  MB_TRUSS3,
  MB_QUAD,
  MB_QUAD4,
  MB_QUAD5,
  MB_QUAD8,
  MB_QUAD9,
  MB_SHELL,
  MB_SHELL4,
  MB_SHELL8,
  MB_SHELL9,
  MB_TRI,
  MB_TRI3,
  MB_TRI6,
  MB_TRI7,
  MB_HEX,
  MB_HEX8,
  MB_HEX9,
  MB_HEX20,
  MB_HEX27,
  MB_PYRAMID,
  MB_TETRA,
  MB_TETRA4,
  MB_TETRA8,
  MB_TETRA10,
  MB_TETRA14,
  MB_KNIFE,
  MB_WEDGE,
  MB_MAX_ELEM_TYPE
};


#endif
