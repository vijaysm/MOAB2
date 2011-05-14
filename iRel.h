#ifndef _ITAPS_iRel
#define _ITAPS_iRel

#define IREL_MAJOR_VERSION 1
#define IREL_MINOR_VERSION 1

  /** \mainpage The ITAPS Relations Interface iRel
   *
   * Each ITAPS interface encapsulates functionality that "belongs"
   * together, for example mesh or geometric model functionality.  In
   * some cases, however, data in several of these interfaces need to
   * be related together.  For example, a collection of mesh faces
   * should be related to the geometric model face which they
   * discretize.  The ITAPS Relations interface accomplishes this in a
   * way which allows the lower-level interfaces to remain
   * independent.
   *
   * iRel defines relations as pairwise relations between entities
   * or entity sets.  Related entities can be in the same or different
   * interfaces.  A given relation is created for a given pair of
   * interfaces and returned in the form of a \em Relation \em Handle.
   * After a specific relation pair has been created, concrete
   * relations for that pair can be assigned and retrieved for
   * specific entities using set and get functions on the iRel
   * interface.  A given interface instance can appear in one or many
   * relation pairs, each identified by the relation pair handle.
   *
   * \section Types Relation Types
   *
   * Relations are also distinguished by a pair of relation types.
   * For each interface in a relation pair, a corresponding type
   * indicates whether the relation applies to entities, entity sets,
   * or both entities and sets in the corresponding interface in the
   * pair.  If only one of the interfaces in a given pair has a
   * 'both'-type, entities and entity sets in that
   * interface are each related to either entities or sets in the other
   * interface in the pair.  If both of the sides of a relation are of 
   * 'both'-type, entities and sets on one side of a relation point to 
   * sets on the other side.
   *
   * \section Status Relation Status
   *
   * Relations are also distinguished by a pair of relation statuses.
   * For each interface in a relation pair, a corresponding status indicates
   * whether the relation on that side is kept up to date, or stored at all.
   * Allowable values for status are iRel_ACTIVE, iRel_INACTIVE, and iRel_NOTEXIST,
   * defined in the iRel_RelationStatus enumeration.  Status for a given side
   * can be changed from iRel_ACTIVE to iRel_INACTIVE and vice versa, or from
   * either of those to iRel_NOTEXIST.  However, once changed to iRel_NOTEXIST
   * (or created that way), a side cannot be changed back to the other two.
   * Changing a side to be iRel_INACTIVE can be used when frequent changes to
   * the underlying entities are being made, e.g. during adaptive mesh refinement.
   * Changing from iRel_INACTIVE to iRel_ACTIVE implies a traversal of all entities
   * on the iRel_ACTIVE side to recover which entities on the iRel_INACTIVE side
   * must have their relations updated.
   *
   * \section ArgOrder Argument Order
   *
   * Many functions in the iRel interface take as input two entities,
   * or two lists of entities, along with a relation pair handle.  For
   * these functions, the entities or lists are assumed to be in the
   * same order as the interfaces used to create that relation pair.
   * For example, if a relation pair is created by calling:
   * \code 
   * iRel_createRelation(instance, iface1, ent_or_set1, type1, 
   *                     iface2, ent_or_set2, type2,
   *                     &relation_handle, &ierr)
   * \endcode
   * and relations set by calling
   * \code
   * iRel_setEntEntRelation(instance, relation_handle,
   *                        ent1, is_set1, ent2, is_set2, &ierr)
   * \endcode
   * it is assumed that ent1 is contained in iface1 and ent2 in
   * iface2.
   *
   * For functions taking only one entity or list as input, and
   * returning an entity or list, an additional argument indicates
   * whether the input entity or list belongs to the first or second
   * interface in that relation pair.
   * 
   */

#include "iBase.h"
#include "iRel_protos.h"

#ifdef __cplusplus

extern "C" 
{
#endif

    /**\brief  Type used to store iRel interface handle
     *
     * Type used to store iRel interface handle
     */
  typedef void* iRel_Instance;

    /**\brief  Type used to store references to relation pairs
     *
     * Type used to store references to relation pairs
     */
  typedef struct iRel_PairHandle_Private* iRel_PairHandle;

    /**\brief  \enum IfaceType Enumerator specifying interface types
     *
     * Enumerator specifying interface types.  This enumeration is
     * necessary because functions to get entities of a given dimension
     * are part of the higher-level interfaces (e.g. iGeom, iMesh) instead
     * of iBase.
     */
  enum iRel_IfaceType 
  {iRel_IfaceType_MIN = 0,
   iRel_IGEOM_IFACE = iRel_IfaceType_MIN, 
   iRel_IMESH_IFACE, 
   iRel_IFIELD_IFACE, 
   iRel_IREL_IFACE,
   iRel_IfaceType_MAX = iRel_IREL_IFACE};

    /**\brief  \enum RelationType Enumerator specifying relation types
     *
     * Enumerator specifying relation types.  A relation has two types, one
     * for each side of the relation.
     */
  enum iRel_RelationType 
  {iRel_RelationType_MIN = 0,
   iRel_ENTITY = iRel_RelationType_MIN, 
   iRel_SET, 
   iRel_BOTH,
   iRel_RelationType_MAX = iRel_BOTH};

    /**\brief  \enum RelationStatus Enumerator specifying relation status
     *
     * Enumerator specifying relation status.  A relation has two statuses, one
     * for each side of the relation.  Allowed values of this enumeration are:
     *    iRel_ACTIVE: the relation on this side is active and up to date
     *    iRel_INACTIVE: the relation on this side is inactive, and may be out of date
     *    iRel_NOTEXIST: the relation on this side is not stored
     * It is an error to request relations from a side that does not have iRel_ACTIVE status
     */
  enum iRel_RelationStatus
  {iRel_RelationStatus_MIN = 0,
   iRel_ACTIVE = iRel_RelationType_MIN, 
   iRel_INACTIVE, 
   iRel_NOTEXIST,
   iRel_RelationStatus_MAX = iRel_NOTEXIST};

    /**\brief  Get the error type returned from the last iRel function
     *
     * Get the error type returned from the last iRel function.  Value
     * returned is a member of the iBase_ErrorType enumeration.
     * \param instance iRel instance handle
     * \param *error_type Error type returned from last iRel function
     */
  void iRel_getErrorType(
    iRel_Instance instance,
    /*out*/ int *error_type);

    /**\brief  Get a description of the error returned from the last iRel
     *         function
     *
     * Get a description of the error returned from the last iRel function
     * \param instance iRel instance handle
     * \param descr Pointer to a character string to be filled with a
     *        description of the error from the last iRel function
     * \param descr_len Length of the character string pointed to by descr
     */
  void iRel_getDescription(
    iRel_Instance instance,
    /*inout*/ char *descr, 
    /*in*/ int descr_len);

    /**\brief  Create a new iRel instance
     *
     * Create a new iRel instance.  Currently no options are implemented.
     * \param options Options for the implementation
     * \param *instance Interface instance
     * \param *err Pointer to error value, returned from function
     * \param options_len Length of options string
     */
  void iRel_create(
    const char *options,
    /*out*/ iRel_Instance *instance,
    /*out*/ int *err,
    /*in*/ const int options_len);

    /**\brief  Destroy the interface object
     *
     * Calls destructor on interface object
     * \param instance Interface object handle to destroy
     * \param *err Pointer to error value, returned from function
     */
  void iRel_destroy(
    iRel_Instance instance,
    /*out*/ int *err);

    /**\brief  Create a relation pair between two interfaces
     *
     * Creates a relation pair between two interfaces, passing
     * back a handle to the pair.  It is an error to create a relation pair
     * having both sides iRel_NOTEXIST.  If a relation pair has a side with status
     * iRel_NOTEXIST, the relation for that side is never stored, and the status
     * cannot change over the life of the relation pair.  
     * \param instance Interface instance
     * \param iface1 1st interface object in the relation pair
     * \param ent_or_set1 This relation relates entities, sets, or both from
     *        1st interface object
     * \param iface_type1 Type of 1st interface
     * \param irel_status1 The status of 1st side
     * \param iface2 2nd interface object in the relation pair
     * \param ent_or_set2 This relation relates entities, sets, or both from
     *        2nd interface object
     * \param iface_type2 Type of 2nd interface
     * \param irel_status2 The status of 2nd side
     * \param *pair Pointer to relation pair handle, returned from function
     * \param *err Pointer to error value, returned from function
     */
  void iRel_createPair(
    iRel_Instance instance,
    /*in*/ iBase_Instance iface1,
    /*in*/ const int ent_or_set1,
    /*in*/ const int iface_type1,
    /*in*/ const int irel_status1,
    /*in*/ iBase_Instance iface2,
    /*in*/ const int ent_or_set2,
    /*in*/ const int iface_type2,
    /*in*/ const int irel_status2,
    /*out*/ iRel_PairHandle *pair,
    /*out*/ int *err);

    /**\brief  Get information for this relation handle
     *
     * Get information about the interfaces and relation type for this
     * relation.  Relation type for each side is passed back as integers,
     * but values will be from RelationType enumeration.
     * \param instance Interface instance
     * \param pair Handle of relation pair being queried
     * \param *iface1 Side 1 instance for this relation
     * \param *ent_or_set1 Relation type for side 1 of this relation
     * \param *iface_type1 Interface type for side 1 of this relation
     * \param *irel_status1 The status of 1st side of this relation
     * \param *iface2 Side 2 instance for this relation
     * \param *ent_or_set2 Relation type for side 2 of this relation
     * \param *iface_type2 Interface type for side 2 of this relation
     * \param *irel_status2 The status of 2nd side of this relation
     * \param *err Pointer to error value, returned from function
     */
  void iRel_getPairInfo(
      iRel_Instance instance,
      /*in*/ iRel_PairHandle pair,
      /*out*/ iBase_Instance *iface1,
      /*out*/ int *ent_or_set1,
      /*out*/ int *iface_type1,
      /*out*/ int *irel_status1,
      /*out*/ iBase_Instance *iface2,
      /*out*/ int *ent_or_set2,
      /*out*/ int *iface_type2,
      /*out*/ int *irel_status2,
      /*out*/ int *err);

    /**\brief  Change the relation type
     *
     * Change the type of one or both sides of a relation.  Only changes that
     * result in no lost information are allowed, e.g. changing a type from SET
     * to BOTH or vice versa.
     *  \param instance Interface instance
     *  \param pair Relation pair handle being changed
     *  \param ent_or_set1 The new type of side 1 of this relation pair
     *  \param ent_or_set2 The new type of side 2 of this relation pair
     *  \param *err Pointer to error value, returned from function
     */
  void iRel_changePairType(
      /*in*/ iRel_Instance instance,
      /*in*/ iRel_PairHandle pair,
      /*in*/ int ent_or_set1,
      /*in*/ int ent_or_set2,
      /*out*/ int *err);

    /**\brief  Change the relation status
     *
     * Change the status of one or both sides of a relation.  It is an error to
     * change the status of both sides to iRel_NOTEXIST.  If a side is changed to
     * iRel_NOTEXIST, it will no longer be changeable back to iRel_ACTIVE or iRel_INACTIVE.
     * Changing a side from iRel_INACTIVE to iRel_ACTIVE implies a traversal of all
     * related entities on the other side, to recover the relations on the side being changed.
     * Changing both sides from iRel_ACTIVE to something else is an error, since in that 
     * case neither will be able to be updated to iRel_ACTIVE.
     *  \param instance Interface instance
     *  \param pair Relation pair handle being changed
     *  \param status1 The new status of side 1 of this relation pair
     *  \param status2 The new status of side 2 of this relation pair
     *  \param *err Pointer to error value, returned from function
     */
  void iRel_changePairStatus(
      /*in*/ iRel_Instance instance,
      /*in*/ iRel_PairHandle pair,
      /*in*/ int irel_status1,
      /*in*/ int irel_status2,
      /*out*/ int *err);

    /**\brief  Destroy a relation pair
     *
     * Destroy the relation pair corresponding to the handle input
     * \param instance Interface instance
     * \param pair Handle of relation pair to destroy
     * \param *err Pointer to error value, returned from function
     */
  void iRel_destroyPair(
    iRel_Instance instance, 
    /*in*/ iRel_PairHandle pair,
    /*out*/ int *err);

    /**\brief  Get relations containing specified interface
     *
     * Get relations containing the specified interface
     * \param instance Interface instance
     * \param iface Specified interface 
     * \param pairs Pointer to array holding returned relation pairs
     *        containing specified interface
     * \param pairs_allocated Pointer to allocated size of relation pairs list
     * \param pairs_size Pointer to occupied size of relation pairs list
     * \param *err Pointer to error value, returned from function
     */
  void iRel_findPairs(
    iRel_Instance instance,
    /*in*/ iBase_Instance iface,
    /*inout*/ iRel_PairHandle **pairs,
    /*inout*/ int *pairs_allocated,
    /*out*/ int *pairs_size,
    /*out*/ int *err);

    /**\brief  Set a relation between two entities
     *
     * Set a relation between an entity and several entities.  It is an error to set
     * a relation on a pair with both sides not iRel_ACTIVE.
     * \param instance Interface instance
     * \param pair Relation pair handle being queried
     * \param ent1 1st entity of relation being set
     * \param ent2 2nd entity of relation being set
     * \param *err Pointer to error value, returned from function
     */
  void iRel_setEntEntRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ iBase_EntityHandle ent2,
    /*out*/ int *err);
  void iRel_setEntSetRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ iBase_EntitySetHandle entset2,
    /*out*/ int *err);
  void iRel_setSetEntRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ iBase_EntityHandle ent2,
    /*out*/ int *err);
  void iRel_setSetSetRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ iBase_EntitySetHandle entset2,
    /*out*/ int *err);

    /**\brief Set relations between arrays of entities pairwise, 
     *        ent_array_1[i]<->ent_array_2[i]
     *
     * Set relations between arrays of entities pairwise, 
     * ent_array_1[i]<->ent_array_2[i].  If either array
     * contains sets and that side of the relation is 'both'-type, 
     * set relations for individual entities in those sets too.  It is an error to set
     * a relation on a pair with both sides not iRel_ACTIVE.
     * \param instance Interface instance
     * \param pair Relation pair handle being queried
     * \param ent_array_1 1st array of entities of relation being set
     * \param num_ent1 Number of entities in 1st array
     * \param ent_array_2 2nd array of entities of relation being set
     * \param num_ent2 Number of entities in 2nd array
     * \param *err Pointer to error value, returned from function
     */
  void iRel_setEntArrEntArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntityHandle *ent_array_2,
    /*in*/ int num_ent2,
    /*out*/ int *err);
  void iRel_setSetArrEntArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntityHandle *ent_array_2,
    /*in*/ int num_ent2,
    /*out*/ int *err);
  void iRel_setEntArrSetArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntitySetHandle *entset_array_2,
    /*in*/ int num_ent2,
    /*out*/ int *err);
  void iRel_setSetArrSetArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntitySetHandle *entset_array_2,
    /*in*/ int num_ent2,
    /*out*/ int *err);

    /**\brief  Get entity related to specified entity and relation handle
     *
     * Get entity related to specified entity and relation handle.  Also
     * returns whether the related entity is an entity or a set.  It is an error to get
     * a relation for a side with status iRel_NOTEXIST.
     * \param instance Interface instance
     * \param pair Relation pair handle being queried
     * \param ent1 1st entity of relation being queried
     * \param switch_order 1st entity is related to 1st interface (=0) or 2nd
     *        interface (=1) of relation pair
     * \param *ent2 Pointer to entity related to ent1
     * \param *err Pointer to error value, returned from function
    */
  void iRel_getEntEntRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ int switch_order,
    /*out*/ iBase_EntityHandle *ent2,
    /*out*/ int *err);
  void iRel_getEntSetRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ int switch_order,
    /*out*/ iBase_EntitySetHandle *entset2,
    /*out*/ int *err);
  void iRel_getSetEntRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ int switch_order,
    /*out*/ iBase_EntityHandle *ent2,
    /*out*/ int *err);
  void iRel_getSetSetRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ int switch_order,
    /*out*/ iBase_EntitySetHandle *entset2,
    /*out*/ int *err);
  void iRel_getEntSetIterRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ int switch_order,
    /*out*/ iBase_EntityIterator *entset2,
    /*out*/ int *err);

    /**\brief  Get entities related to those in specified array and relation,
     *         pairwise
     *
     * Get entities related to those in specified array and relation, pairwise.
     * Returns sets or entities, depending on relation type and entities in 
     * ent_array_1.  It is an error to get a relation for a side with status iRel_NOTEXIST.
     * \param instance Interface instance
     * \param pair Relation pair handle being queried
     * \param ent_array_1 Array of entities whose relations are being queried
     * \param ent_array_1_size Number of entities in ent_array_1
     * \param switch_order Entities in ent_array_1 are related with 1st (=0) 
     *        or 2nd (=1) interface of this relation pair
     * \param *ent_array_2 Pointer to array of entity handles returned from
     *        function
     * \param *ent_array_2_allocated Pointer to allocated size of ent_array_2
     * \param *ent_array_2_size Pointer to occupied size of ent_array_2
     * \param *err Pointer to error value, returned from function
     */
  void iRel_getEntArrEntArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int ent_array_1_size,
    /*in*/ int switch_order,
    /*inout*/ iBase_EntityHandle **ent_array_2,
    /*inout*/ int *ent_array_2_allocated,
    /*out*/ int *ent_array_2_size,
    /*out*/ int *err);
  void iRel_getEntArrSetArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int ent_array_1_size,
    /*in*/ int switch_order,
    /*inout*/ iBase_EntitySetHandle **entset_array_2,
    /*inout*/ int *entset_array_2_allocated,
    /*out*/ int *entset_array_2_size,
    /*out*/ int *err);
  void iRel_getSetArrEntArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int entset_array_1_size,
    /*in*/ int switch_order,
    /*inout*/ iBase_EntityHandle **ent_array_2,
    /*inout*/ int *ent_array_2_allocated,
    /*out*/ int *ent_array_2_size,
    /*out*/ int *err);
  void iRel_getSetArrSetArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int entset_array_1_size,
    /*in*/ int switch_order,
    /*inout*/ iBase_EntitySetHandle **entset_array_2,
    /*inout*/ int *entset_array_2_allocated,
    /*out*/ int *entset_array_2_size,
    /*out*/ int *err);
  void iRel_getEntArrSetIterArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int ent_array_1_size,
    /*in*/ int switch_order,
    /*inout*/ iBase_EntityIterator **entiter,
    /*inout*/ int *entiter_allocated,
    /*out*/ int *entiter_size,
    /*out*/ int *err);

    /**\brief  Remove a relation between two entities
     *
     * Remove a relation between an entity and several entities.
     * \param instance Interface instance
     * \param pair Relation pair handle being queried
     * \param ent1 1st entity of relation being removed
     * \param ent2 2nd entity of relation being removed
     * \param *err Pointer to error value, returned from function
     */
  void iRel_rmvEntEntRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ iBase_EntityHandle ent2,
    /*out*/ int *err);
  void iRel_rmvEntSetRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,
    /*in*/ iBase_EntityHandle ent1,
    /*in*/ iBase_EntitySetHandle entset2,
    /*out*/ int *err);
  void iRel_rmvSetEntRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ iBase_EntityHandle ent2,
    /*out*/ int *err);
  void iRel_rmvSetSetRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,
    /*in*/ iBase_EntitySetHandle entset1,
    /*in*/ iBase_EntitySetHandle entset2,
    /*out*/ int *err);

    /**\brief Remove relations between arrays of entities pairwise, 
     *        ent_array_1[i]<->ent_array_2[i]
     *
     * Remove relations between arrays of entities pairwise, 
     * ent_array_1[i]<->ent_array_2[i].  If either array
     * contains sets and that side of the relation is 'both'-type, 
     * remove relations for individual entities in those sets too.
     * \param instance Interface instance
     * \param pair Relation pair handle being queried
     * \param ent_array_1 1st array of entities of relation being removed
     * \param num_ent1 Number of entities in 1st array
     * \param ent_array_2 2nd array of entities of relation being removed
     * \param num_ent2 Number of entities in 2nd array
     * \param *err Pointer to error value, returned from function
     */
  void iRel_rmvEntArrEntArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntityHandle *ent_array_2,
    /*in*/ int num_ent2,
    /*out*/ int *err);
  void iRel_rmvSetArrEntArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntityHandle *ent_array_2,
    /*in*/ int num_ent2,
    /*out*/ int *err);
  void iRel_rmvEntArrSetArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle *ent_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntitySetHandle *entset_array_2,
    /*in*/ int num_ent2,
    /*out*/ int *err);
  void iRel_rmvSetArrSetArrRelation(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle *entset_array_1,
    /*in*/ int num_ent1,
    /*in*/ iBase_EntitySetHandle *entset_array_2,
    /*in*/ int num_ent2,
    /*out*/ int *err);

    /**\brief  Infer relations between entities in specified pair of interfaces
     *
     * Infer relations between entities in specified pair of interfaces.  The
     * criteria used to infer these relations depends on the interfaces in
     * the pair, the iRel implementation, and the source of the data in those
     * interfaces.
     * \param instance Interface instance
     * \param pair Relation pair handle being queried
     * \param *err Pointer to error value, returned from function
     */
  void iRel_inferAllRelations(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,
    /*out*/ int *err);

    /**\brief  Infer relations and relation type between entities in specified 
     *         pair of interfaces 
     *
     * Infer relations between entities in specified pair of interfaces, and the
     * relation type used by this iRel implementation.  The criteria used to
     * infer these relations depends on the interfaces in the pair, the iRel
     * implementation, and the source of the data in those interfaces.
     * \param instance Interface instance
     * \param pair Relation pair handle created by implementation
     * \param *err Pointer to error value, returned from function
     */
  void iRel_inferAllRelationsAndType(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle *pair,
    /*out*/ int *err);

    /**\brief  Infer relations corresponding to specified entity and relation
     *         pair
     *
     * Infer relations corresponding to specified entity and relation pair.  The
     * criteria used to infer these relations depends on the interfaces in
     * the pair, the iRel implementation, and the source of the data in those
     * interfaces.
     * \param instance Interface instance
     * \param pair Relation pair handle being queried
     * \param entity Entity whose relations are being inferred
     * \param iface_no Entity corresponds to 1st (=0) or 2nd (=1) interface
     *        in relation pair
     * \param *err Pointer to error value, returned from function
     */
  void iRel_inferEntRelations(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle entity,
    /*in*/ int iface_no,
    /*out*/ int *err);
  void iRel_inferSetRelations(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle entity_set,
    /*in*/ int iface_no,
    /*out*/ int *err);

    /**\brief  Infer relations corresponding to specified entities and relation
     *         pair
     *
     * Infer relations corresponding to specified entities and relation pair.
     * The criteria used to infer these relations depends on the interfaces in
     * the pair, the iRel implementation, and the source of the data in those
     * interfaces.
     * \param instance Interface instance
     * \param pair Relation pair handle being queried
     * \param entities Array of entities whose relation are being inferred
     * \param entities_size Number of entities in array
     * \param iface_no Entities correspond to 1st (=0) or 2nd (=1) interface
     *        in relation pair
     * \param *err Pointer to error value, returned from function
     */
  void iRel_inferEntArrRelations(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntityHandle *entities,
    /*in*/ int entities_size,
    /*in*/ int iface_no,
    /*out*/ int *err);
  void iRel_inferSetArrRelations(
    iRel_Instance instance,
    /*in*/ iRel_PairHandle pair,    
    /*in*/ iBase_EntitySetHandle *entity_sets,
    /*in*/ int entities_size,
    /*in*/ int iface_no,
    /*out*/ int *err);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* #ifndef _ITAPS_iRel */

