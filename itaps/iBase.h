#ifndef _ITAPS_iBase
#define _ITAPS_iBase

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compile time version number digits
 *
 * iBase maintains a major, minor and patch digit in its version number.
 * Technically speaking, there is not much practical value in patch digit
 * for an interface specification. A patch release is typically only used
 * for bug fix releases. Although it is rare, sometimes a bug fix
 * necessitates an API change. So, we define a patch digit for iMesh.
 *
 * Although each interface in ITAPS has been designed to support its own
 * uniqe version numbers, apart from other ITAPS interfaces, as currently
 * used, we require all ITAPS interfaces to use the same ITAPS-wide version
 * number derived from the version number defined by these three digits.
 ******************************************************************************/
#define IBASE_VERSION_MAJOR 1
#define IBASE_VERSION_MINOR 4
#define IBASE_VERSION_PATCH 0

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Version Comparison
 *
 * Evaluates to true at CPP time if the version of iBase currently being
 * compiled is greater than or equal to the version specified.
 ******************************************************************************/
#define IBASE_VERSION_GE(Maj,Min,Pat) \
    (((IBASE_VERSION_MAJOR==(Maj)) && (IBASE_VERSION_MINOR==(Min)) && (IBASE_VERSION_PATCH>=(Pat))) || \
     ((IBASE_VERSION_MAJOR==(Maj)) && (IBASE_VERSION_MINOR>(Min))) || \
      (IBASE_VERSION_MAJOR>(Maj)))

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compose compile-time string represention of the version number
 ******************************************************************************/
#define IBASE_VERSION_STRING___(I,X,Y,Z) #I "_Version_" #X "." #Y "." #Z
#define IBASE_VERSION_STRING__(I,X,Y,Z) IBASE_VERSION_STRING___(I,X,Y,Z)
#define IBASE_VERSION_STRING_(I) IBASE_VERSION_STRING__(I,IBASE_VERSION_MAJOR,IBASE_VERSION_MINOR,IBASE_VERSION_PATCH)
#define IBASE_VERSION_STRING IBASE_VERSION_STRING_(iBase)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compose compile-time symbol name derived from the version number.
 ******************************************************************************/
#define IBASE_VERSION_TAG__(I,X,Y,Z) I##_Version_##X##_##Y##_##Z
#define IBASE_VERSION_TAG_(I,X,Y,Z) IBASE_VERSION_TAG__(I,X,Y,Z)
#define IBASE_VERSION_TAG(I) IBASE_VERSION_TAG_(I,IBASE_VERSION_MAJOR,IBASE_VERSION_MINOR,IBASE_VERSION_PATCH)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief ITAPS-wide (across all ITAPS APIs) version handling
 ******************************************************************************/
#define ITAPS_VERSION_MAJOR IBASE_VERSION_MAJOR
#define ITAPS_VERSION_MINOR IBASE_VERSION_MINOR
#define ITAPS_VERSION_PATCH IBASE_VERSION_PATCH
#define ITAPS_VERSION_GE(Maj,Min,Pat) IBASE_VERSION_GE(Maj,Min,Pat)
#define ITAPS_VERSION_STRING_(I) IBASE_VERSION_STRING_(I)
#define ITAPS_VERSION_STRING ITAPS_VERSION_STRING_(ITAPS)
#define ITAPS_VERSION_TAG_(I) IBASE_VERSION_TAG(I)
#define ITAPS_VERSION_TAG ITAPS_VERSION_TAG_(I)

/***************************************************************************//**
 * \defgroup EnumIterators Enum-Iterators
 * \ingroup iBase
 * \brief Convenience macros for iterating over all possible values in an enum
 *
 * These convenience macros are provided to facilitate iterating over all
 * possible values in an enumerated type. To use these macros, for example...
 * \code 
 * for (iBase_EntityType i  = IBASE_MINENUM(iBase_EntityType);
 *                       i <= IBASE_MAXENUM(iBase_EntityType);
 *                            IBASE_INCENUM(i,iBase_EntityType))
 * {
 * }
 * \endcode
 * Be aware that some enumerated types include a <em>wild card</em> often used
 * in queries to represent all possible values and you may or may not want to
 * include such a value in your iteration.
 ******************************************************************************/

/***************************************************************************//**
 * \ingroup EnumIterators
 * @{
 ******************************************************************************/
#define IBASE_MINENUM(enumName) enumName ## _MIN
#define IBASE_MAXENUM(enumName) enumName ## _MAX
#define IBASE_NUMENUM(enumName) ((int)IBASE_MAXENUM(enumName) - (int)IBASE_MINENUM(enumName) + 1)
#define IBASE_INCENUM(enumName,I) (I = (enum enumName)((int)I+1))
/** @} */

#ifdef __cplusplus
extern "C" {
#endif

typedef void* iBase_Instance;
typedef struct iBase_EntityHandle_Private* iBase_EntityHandle;
typedef struct iBase_EntitySetHandle_Private* iBase_EntitySetHandle;
typedef struct iBase_TagHandle_Private* iBase_TagHandle;
typedef struct iBase_EntityIterator_Private* iBase_EntityIterator;
typedef struct iBase_EntityArrIterator_Private* iBase_EntityArrIterator;

enum iBase_EntityType {
    iBase_EntityType_MIN = 0,
        /**< facilitates iteration over all values */
    iBase_VERTEX = iBase_EntityType_MIN,
        /**< A topological dimension 0 entity */
    iBase_EDGE,
        /**< A topological dimension 1 entity */
    iBase_FACE,
        /**< A topological dimension 2 entity */
    iBase_REGION,
        /**< A topological dimension 3 entity */
    iBase_ALL_TYPES,
        /**< used only in queires to request information about all types */
    iBase_EntityType_MAX = iBase_ALL_TYPES
        /**< facilitates iteration over all values */
};

enum iBase_AdjacencyCost {
    iBase_AdjacencyCost_MIN = 0,
        /**< facilitates iteration over all values */
    iBase_UNAVAILABLE = iBase_AdjacencyCost_MIN,
        /**< Adjacency information not supported */
    iBase_ALL_ORDER_1,
        /**< No more than local mesh traversal required (i!=j) */
    iBase_ALL_ORDER_LOGN,
        /**< Global tree search (i!=j) */
    iBase_ALL_ORDER_N,
        /**< Global exhaustive search (i!=j) */
    iBase_SOME_ORDER_1,
        /**< Only some adjacency info, local (i!=j) */
    iBase_SOME_ORDER_LOGN,
        /**< Only some adjacency info, tree (i!=j) */
    iBase_SOME_ORDER_N,
        /**< Only some adjacency info, exhaustive (i!=j) */
    iBase_AVAILABLE,
        /**< ALL (intermediate) entities available. (i==j) */
    iBase_AdjacencyCost_MAX = iBase_AVAILABLE
        /**< facilitates iteration over all values */
};

enum iBase_CreationStatus {
    iBase_CreationStatus_MIN = 0,
        /**< facilitates iteration over all values */
    iBase_NEW = iBase_CreationStatus_MIN,
        /**< The entity was newly created */
    iBase_ALREADY_EXISTED,
        /**< The entity already existed and the handle for that
             already existing handle was returned */
    iBase_CREATED_DUPLICATE,
        /**< The entity already existed but a new, duplicate entity was
             nevertheless created */
    iBase_CREATION_FAILED,
        /**< Creation of the entity did not succeed */
    iBase_CreationStatus_MAX = iBase_CREATION_FAILED
        /**< facilitates iteration over all values */
};

enum iBase_ErrorType {
    iBase_ErrorType_MIN = 0,
        /**< facilitates iteration over all values */
    iBase_SUCCESS = iBase_ErrorType_MIN,
    iBase_MESH_ALREADY_LOADED,
    iBase_FILE_NOT_FOUND,
    iBase_FILE_WRITE_ERROR,
    iBase_NIL_ARRAY,
    iBase_BAD_ARRAY_SIZE,
    iBase_BAD_ARRAY_DIMENSION,
    iBase_INVALID_ENTITY_HANDLE,
    iBase_INVALID_ENTITY_COUNT,
    iBase_INVALID_ENTITY_TYPE,
    iBase_INVALID_ENTITY_TOPOLOGY,
    iBase_BAD_TYPE_AND_TOPO,
    iBase_ENTITY_CREATION_ERROR,
    iBase_INVALID_TAG_HANDLE,
    iBase_TAG_NOT_FOUND,
    iBase_TAG_ALREADY_EXISTS,
    iBase_TAG_IN_USE,
    iBase_INVALID_ENTITYSET_HANDLE,
    iBase_INVALID_ITERATOR_HANDLE,
    iBase_INVALID_ARGUMENT,
    iBase_MEMORY_ALLOCATION_FAILED,
    iBase_NOT_SUPPORTED,
    iBase_FAILURE,
    iBase_ErrorType_MAX = iBase_FAILURE
        /**< facilitates iteration over all values */
  };

/***************************************************************************//**
 * \details
 * Many of the functions in iMesh can return arrays of tuples; that is, arrays
 * of multi-valued type. For example, the function iMesh_getVtxArrCoords,
 * returns an array of xyz coordinate 3-tuples (or, perhaps for geometrically
 * 2D meshes, xy 2-tuples). In these situations, there are multiple ways the
 * data can be organized in memory. For example, it could be stored xyz,xyz,xyz
 * or xxx...,yyy...,zzz.... These two different storage orders are referred
 * to as INTERLEAVED and BLOCKED, respectively. For some functions in iMesh,
 * the storage order is explicitly specified as an argument to the function.
 * For other functions, the storage order is not explicitly specified. And,
 * in these cases, it shall always be implicitly assumed to be INTERLEAVED.
 * This fact will be mentioned in the documentation for each specific function
 * where it applies. For example, in case of iMesh_getEntArrAdj, the returned
 * array of adjacent entities is multi-valued in that it stores for each
 * entity queried, all its adjacent entities. Such an array will be stored
 * INTERLEAVED with all adjacent entities for the first entity in the query
 * followed by all adjacent entities for the second entity in the query and
 * so forth.
 ******************************************************************************/
enum iBase_StorageOrder {
    iBase_StorageOrder_MIN = 0,
        /**< facilitates iteration over all values */
    iBase_BLOCKED = iBase_StorageOrder_MIN,
        /**< xxx...yyy...zzz... */
    iBase_INTERLEAVED,
        /**< xyzxyzxyz... */
    iBase_StorageOrder_MAX = iBase_INTERLEAVED
        /**< facilitates iteration over all values */
};

enum iBase_TagValueType {
    iBase_TagValueType_MIN = 0,
        /**< facilitates iteration over all values */
    iBase_BYTES = iBase_TagValueType_MIN,
        /**< An opaque sequence of bytes, size always measured in bytes */
    iBase_INTEGER,
        /**< A value of type \c int */
    iBase_DOUBLE,
        /**< A value of type \c double */
    iBase_ENTITY_HANDLE,
        /**< A value of type \c iBase_EntityHandle */
    iBase_ENTITY_SET_HANDLE,
        /**< A value of type \c iBase_EntitySetHandle */
    iBase_TagValueType_MAX = iBase_ENTITY_SET_HANDLE
        /**< facilitates iteration over all values */
};

/***************************************************************************//**
 * \page ibase iBase: ITAPS Base Interface
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iBase iBase
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup VersionNumbers Version Numbers
 * \ingroup iBase
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup ErrorHandling Error Handling
 * \ingroup iBase
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup Datatypes Datatypes
 * \ingroup iBase
 ******************************************************************************/

/***************************************************************************//**
 * \mainpage The ITAPS Interfaces
 *
 * \subpage ibase
 *
 * \subpage imesh
 *
 * \subpage imeshp
 * 
 * \subpage igeom
 *
 * \subpage error
 *
 * \subpage trio
 *
 * \subpage strlen
 *
 * \subpage options
 *
 * \subpage numhops
 *
 * \page error Error Handling
 *
 * With few exceptions, every iMesh function includes an output argument,
 * 'int *err', which returns an error code indicating if the function call
 * may have failed. If the value returned for the 'err' argument is NOT
 * iBase_SUCCESS, the caller should NOT attempt to interpret (read the
 * values in) any of the other return arguments of the call. While some
 * implementations may actually return valid/useful results in other
 * return arguments of a call that has failed, there is no guarentee that
 * ALL implementations will do similarly and so depending on such behavior
 * is neither portable nor safe. This is true even if the returned values
 * are different from the values of the arguments before the call was
 * made.
 *
 * \page trio Array pointer, allocated and occupied sizes argument trio
 *
 * Many of the functions in iMesh have arguments corresponding to lists of 
 * objects.  In-type arguments for lists consist of a pointer to an array and
 * a list size.  Lists returned from functions are passed in three arguments,
 * a pointer to the array representing the list, and pointers to the
 * allocated and occupied lengths of the array.  These three arguments are 
 * inout-type arguments, because they can be allocated by the application and
 * passed into the interface to hold the results of the function.  Lists
 * which are pre-allocated must be large enough to hold the results of the
 * function; if this is not the case, an error is generated.  Otherwise, the
 * occupied size is changed to the size output from the function.  If a list
 * argument is unallocated (the list pointer points to a NULL value) or if
 * the incoming value of the allocated size is zero, the list storage will be
 * allocated by the implementation.
 *
 * IN ALL CASES, MEMORY ALLOCATED BY ITAPS INTERFACE IMPLEMENTATIONS IS DONE
 * USING THE C MALLOC FUNCTION, AND MUST BE DE-ALLOCATED USING THE C FREE
 * FUNCTION.
 *
 * \page strlen String Length Arguments
 *
 * Many of the functions in iMesh involve passing a string and also the length
 * of that string. How is the null character is handled?
 * For users of the iMesh interface calling iMesh functions, it is optional
 * as to whether or not to include the null character in computing the length
 * of the string. So, for example, calling iMesh from a C program, users could
 * pass strlen(my_string) or strlen(my_string)+1 as the length of the string.
 *
 * <em>Note to implementors</em>: However, it should be noted that the situation
 * is different for implementers of the iMesh interface. In implementing an
 * iMesh interface function, there can be no assumption that the string is
 * indeed null terminated. The length argument the caller passes in may or may
 * NOT include the null character and implementations must be coded to
 * accommodate this. This requirement is primarily due to differences in how
 * Fortran and C/C++ handle passing of strings as function arguments.
 * Furthermore, because of the way Fortran clients pass strings (Fortran always
 * passes the length of the string as declared in the source code), there
 * may be trailing spaces in the string that need to be truncated.
 *
 * \page numhops Indirection in Set-Inclusion and Parent-Child structures
 *
 * Various functions to query entities, entity sets and parent or child sets 
 * as well as the numbers of these involve a num_hops argument. If the set
 * upon which the query is originated is the root set, the num_hops argument
 * is irrelevant and is ignored by the implementation. Otherwise, the num_hops
 * argument represents the maximum number of levels of indirection employed in
 * satisfying the query not including the originating set. For example, using
 * value for num_hops of 0 (zero) in iMesh_getEntSets will return all the 
 * entity sets that are immediately contained in a given set. Likewise, a
 * value for num_hops of 1 (one) will return all entity sets that are 
 * immediately contained in the given set plus all entity sets that
 * are contained in those immediately contained sets (e.g. one level of
 * indirection). Using a value of -1 for num_hops will return results for
 * all possible levels of indirection. In other words, using a value of
 * -1 for num_hops is equivalent to setting the maximum number of levels
 * of indirection to infinity.
 *
 * \page options Option Strings
 *
 * A few of the functions in iMesh support arbitrary options passed as a
 * character string, called an 'Option String'. The format of and handling
 * of an Option String is as follows...
 *
 * 1. Option Strings are INsensitive to case.
 *
 * 2. Each option in an Option String is pre-pended with the implementation
 * name followed by a special character called the separator character.
 *
 * 3. The separator is a colon, ':'.
 *
 * 4. Multiple options existing in a single Option String are separated by a
 * special character called the delimiter character.
 *
 * 5. The delimiter character is a space, ' '.
 *
 * 6. The effect of multiple options in a single Option String is 
 * INsensitive to order of occurrence in the string.
 *
 * 7. By default, implementations silently ignore any options that
 * do not match on the implementation name part (everything before
 * the separator character). This way, a caller may included options
 * in a single string intended for multiple different implementations.
 *
 * 8. Implementations may (or may not) warn or error for option strings
 * that match on implementation name part but are found to be in error
 * for other reasons the implementation decides.
 *
 * 9. Whenever either the separator character, ':', or delimiter character,
 * ' ', need to appear in an option, they must be escaped with the
 * backslash character, '\'.
 *
 * For example, consider the Options String
 *
 *     "grummp:silant FMDB:TwoPhaseIO moab:mpiio_hints\ foo\:bar"
 *
 * In the above example, the space serves as the delimiter character
 * between multiple options in the string. The colon serves as the
 * implementation-name/option separator character. Because options are
 * required to be insensitive to case, the caller is free to use case as a
 * word separator as in 'TwoPhaseIO' and even in the implementation name,
 * as in 'FMDB:', although 'fmdb:twophaseio' and 'fmdb:TWOPHASEIO' would
 * all have the same effect. In the moab option, both the separator
 * character and delimiter character appear in the option and so are
 * pre-pended (e.g. escaped) with the backslash character.
 
 * GRUMMP will silently ignore the FMDB: and moab: options because they do
 * NOT match on the implementation name part. However, GRUMMP may
 * optionally error out, or warn or silently ignore 'grummp:silant' (it was
 * supposed to be spelled 'silent') as an invalid option.
 *
 * Note that iMesh itself currently does not define any options. In order
 * to discover options a given implementation defines, users are directed
 * to the developers of the respective implementations.
 ******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif /* #ifndef _ITAPS_iBase */
