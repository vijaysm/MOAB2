/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#ifndef MOAB_TYPES_HPP
#define MOAB_TYPES_HPP

#ifdef __cplusplus
#include "moab/EntityType.hpp"
#include "moab/EntityHandle.hpp"
#endif

/**\name Types and names
 * Types used in the MOAB interface
 */
/*@{*/

#ifdef __cplusplus
namespace moab {
#endif

#ifdef WIN32
#ifdef MB_EXPORTS
#define MB_DLL_EXPORT __declspec(dllexport)
#else
#define MB_DLL_EXPORT
#endif
#else
#define MB_DLL_EXPORT
#endif

/** MOAB error codes */
enum ErrorCode { MB_SUCCESS = 0,
                   MB_INDEX_OUT_OF_RANGE,
                   MB_TYPE_OUT_OF_RANGE,
                   MB_MEMORY_ALLOCATION_FAILED,
                   MB_ENTITY_NOT_FOUND,
                   MB_MULTIPLE_ENTITIES_FOUND,
                   MB_TAG_NOT_FOUND,
                   MB_FILE_DOES_NOT_EXIST,
                   MB_FILE_WRITE_ERROR,
                   MB_NOT_IMPLEMENTED,
                   MB_ALREADY_ALLOCATED,
                   MB_VARIABLE_DATA_LENGTH,
                   MB_INVALID_SIZE,
                   MB_UNSUPPORTED_OPERATION,
                   MB_UNHANDLED_OPTION,
                   MB_FAILURE};

/** Misc. integer constants, declared in enum for portability */
enum Constants {
  MB_VARIABLE_LENGTH = -1 /**< Length value for variable-length tags */ 
};

/** Specify storage type for tags.  See MOAB users guide for more information. */
enum TagType {
  MB_TAG_BIT = 0, /**< size measured in bits instead of bytes, otherwise identical to sparse */
  MB_TAG_SPARSE,  /**< tags stored in (entity handle, tag value) pairs */
  MB_TAG_DENSE,   /**< tags stored in vectors directly on entity sequences, cheaper for tags which go on lots of entities */ 
  MB_TAG_MESH, 
  MB_TAG_LAST=MB_TAG_MESH};

/** Specify data type for tags. */
enum DataType {
  MB_TYPE_OPAQUE  = 0, /**< byte array */
  MB_TYPE_INTEGER = 1, /**< native 'int' type */
  MB_TYPE_DOUBLE  = 2, /**< native 'double' type */
  MB_TYPE_BIT     = 3, /**< mandatory type for tags with MB_TAG_BIT storage */
  MB_TYPE_HANDLE  = 4, /**< EntityHandle */
  MB_MAX_DATA_TYPE = MB_TYPE_HANDLE
};

/** Used to reference tags; since they're so different from entities, we
 *  use void** instead of a uint to prevent them from being confused as 
 *  entity handles.
 */
typedef void** Tag;

/** Meshset options: properties for meshset creation.
 *  Values are bit flags that may be combined with a bitwise OR (|)
 */
enum EntitySetProperty {
  MESHSET_TRACK_OWNER = 0x1, /**< create entity to meshset adjacencies */
  MESHSET_SET         = 0x2, /**< set contents are unique */
  MESHSET_ORDERED     = 0x4  /**< order of set contents is preserved */
};

#ifdef __cplusplus
} /* namespace moab */
#endif

/*@}*/

#endif
