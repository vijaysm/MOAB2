//! Conventional tag names used for some often-used sets

/* MATERIAL_SET_TAG_NAME tag:
 * Represents sets of elements having a common material (corresponds to
 * element blocks in ExodusII)
 * size = sizeof(int)
 * type = int
 * value = integer id for this set (block id from ExodusII)
 */
#define MATERIAL_SET_TAG_NAME  "MATERIAL_SET"

/* DIRICHLET_SET_TAG_NAME tag:
 * Represents dirichlet-type boundary condition, usually contains only mesh vertices
 * (corresponds to nodesets in ExodusII)
 * size = sizeof(int)
 * type = int
 * value = integer id for this set (nodeset id from ExodusII)
 */
#define DIRICHLET_SET_TAG_NAME "DIRICHLET_SET"

/* NEUMANN_SET_TAG_NAME  tag:
 * Represents neumann-type boundary condition, usually contains elements with dimension
 * one lower than those found in material sets (i.e. edges in FE quad/tri models, quads/tris
 * in FE hex/tet models) (corresponds to sidesets in ExodusII)
 * size = sizeof(int)
 * type = int
 * value = integer id for this set (sideset id from ExodusII)
 */
#define NEUMANN_SET_TAG_NAME   "NEUMANN_SET"

/* HAS_MID_NODES_TAG_NAM tag:
 * Flags telling whether elements in a given set have mid-(edge, face, region) vertices/nodes;
 * index 0 is a place holder, so this datum can be indexed by dimension, e.g. has_mid_nodes[dim]
 * indicates whether mesh entities of dimension dim have mid nodes
 * size = 4*sizeof(int)
 * type = int[4]
 * value = 1 (has mid nodes), 0 (does not have mid nodes)
 */
#define HAS_MID_NODES_TAG_NAME "HAS_MID_NODES"

/* GEOM_DIMENSION tag: 
 * Represents entities "owned" by a given topological entity in a geometric model
 * size = sizeof(int)
 * type = int
 * value = dimension of geom entity 
 */
#define GEOM_DIMENSION_TAG_NAME "GEOM_DIMENSION"

/* MESH_TRANSFORM tag:
 * Represents homogeneous transform to be applied to mesh; used in ExodusII writer to apply
 * transform before writing nodal coordinates
 * size = 16*sizeof(double)
 * type = double[16]
 * value = 4x4 homogenous transform matrix
 */
#define MESH_TRANSFORM_TAG_NAME "MESH_TRANSFORM"

/* GLOBAL_ID tag:
 * Represents global id of entities (sets or mesh entities); this id is different than the id
 * embedded in the entity handle
 * size = sizeof(int)
 * type = int
 * value = global id
 */
#define GLOBAL_ID_TAG_NAME "GLOBAL_ID"

/* CATEGORY tag:
 * String name indicating generic "category" if the entity to which it is assigned (usually
 * sets); used e.g. to indicate a set represents geometric vertex/edge/face/region, 
 * dual surface/curve, etc.
 * size = CATEGORY_TAG_NAME_LENGTH (defined below)
 * type = char[CATEGORY_TAG_NAME_LENGTH]
 * value = NULL-terminated string denoting category name
 */
#define CATEGORY_TAG_NAME "CATEGORY"
#define CATEGORY_TAG_NAME_LENGTH 32

