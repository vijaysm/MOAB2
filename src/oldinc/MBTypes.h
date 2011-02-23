#ifndef MBTypes_HEADER
#define MBTypes_HEADER

#include "moab/Types.hpp"
#include "MBEntityType.h"
#include "MBEntityHandle.h"

typedef moab::Tag MBTag;

typedef moab::ErrorCode MBErrorCode;
using moab::MB_SUCCESS;
using moab::MB_INDEX_OUT_OF_RANGE;
using moab::MB_MEMORY_ALLOCATION_FAILED;
using moab::MB_ENTITY_NOT_FOUND;
using moab::MB_MULTIPLE_ENTITIES_FOUND;
using moab::MB_TAG_NOT_FOUND;
using moab::MB_FILE_DOES_NOT_EXIST;
using moab::MB_FILE_WRITE_ERROR;
using moab::MB_NOT_IMPLEMENTED;
using moab::MB_ALREADY_ALLOCATED;
using moab::MB_VARIABLE_DATA_LENGTH;
using moab::MB_INVALID_SIZE;
using moab::MB_UNSUPPORTED_OPERATION;
using moab::MB_UNHANDLED_OPTION;
using moab::MB_FAILURE;

typedef moab::Constants MBConstants;
using moab::MB_VARIABLE_LENGTH;

typedef moab::TagType MBTagType;
using moab::MB_TAG_BIT;
using moab::MB_TAG_SPARSE;
using moab::MB_TAG_DENSE;
using moab::MB_TAG_MESH;

typedef moab::DataType MBDataType;
using moab::MB_TYPE_OPAQUE;
using moab::MB_TYPE_INTEGER;
using moab::MB_TYPE_DOUBLE;
using moab::MB_TYPE_BIT;
using moab::MB_TYPE_HANDLE;
using moab::MB_MAX_DATA_TYPE;

typedef moab::EntitySetProperty MBEntitySetProperty;
using moab::MESHSET_TRACK_OWNER;
using moab::MESHSET_SET;
using moab::MESHSET_ORDERED;

#endif
