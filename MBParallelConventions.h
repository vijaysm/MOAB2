#ifndef MB_PARALLEL_CONVENTIONS_H
#define MB_PARALLEL_CONVENTIONS_H

/** Tag conventions for naming parallel things.  Note this header
 * file belongs in the main MOAB directory because even serial
 * applications (e.g. partitioners) may write tags for use in
 * parallel applications.
 */

/** \brief Meshset tag name for interfaces between processors
 *
 * Meshset containing the interface between two processors.
 * The contents of the mesh set is a group of mesh sets
 * (typically corresponding to geometric topology) containing
 * mesh entities.
 *
 * The data of the tag is struct { int; int; }; and contains 
 * the processor numbers (MPI rank) for the pair of processors
 * sharing the interface.
 */
#define PARALLEL_INTERFACE_TAG_NAME "Interface"

/** \brief Global identifier for interface mesh
 *
 * An integer identifier common to the corresponding mesh entity
 * instances on each processor for a mesh entity on the interface.
 */
#define PARALLEL_GLOBAL_ID_TAG_NAME "GlobalHandle"

/** \brief Dimension of a meshset corresponding to geometric topology.
 *
 * A tag to identify a meshset that corresponds to some geometric
 * topology.  The value is an integer containing the dimension of
 * the corresponding geomteric entity (0=vertex, 1=curve, etc.)
 */
#define PARALLEL_GEOM_TOPO_TAG_NAME "GeometryTopology"

/** \brief Tag on a meshset representing a parallel partition.
 *
 * When the mesh is partitioned for use in a parallel environment,
 * the each CPUs partiiton of the mesh is stored in a meshset with
 * this tag.  The value of the tag is an integer containing the
 * processor ID (MPI rank).
 */
#define PARALLEL_PARTITION_TAG_NAME "PARALLEL_PARTITION"
 
/** \brief Tag storing which other processor a given entity is shared with
 *
 * This single-valued tag implies an entity is shared with one other proc
 */
#define PARALLEL_SHARED_PROC_TAG_NAME "PARALLEL_SHARED_PROC"
 
/** \brief Tag storing which other processorS a given entity is shared with
 *
 * This multiple-valued tag implies an entity is shared with multiple
 * other processors.  Length of tag is application-dependent, and depends on
 * what the maximum number of processors is which share an entity
 */
#define PARALLEL_SHARED_PROCS_TAG_NAME "PARALLEL_SHARED_PROCS"
 
#endif
