/*! \page metadata I/O and Meta-Data Storage Conventions in MOAB

	<Center> <H3>   Timothy J. Tautges </H3> </Center>

  	\subpage  md-contents

  	\subpage  md-tables
*/

/*!  \page md-contents Table of Contents

  \ref meta-introduction

  \ref meta-conventions

  \ref meta-options

  \ref meta-references

  \ref appendixA

  \ref appendixB

  \ref appendixC

  \ref appendixD

  \ref appendixE

  \section meta-introduction  Introduction

The Mesh-Oriented datABase (MOAB) is a library for representing finite element and other types of mesh data [1].  Various types of meta-data are often used in conjunction with a mesh.  Examples include boundary condition groupings, material types, and provenance information for the mesh.  Because the data model used in MOAB is so abstract, conventions are useful for describing how meta-data is stored into that data model.  This document describes those conventions for several types of data commonly found in meshes stored in MOAB.  Because the data models used by MOAB and iMesh, the ITAPS mesh interface [2], are so similar, the conventions described here apply almost unmodified to iMesh as well as to MOAB.

The meshes represented in MOAB originate in a variety of forms, including mesh read from files of various formats (e.g. CUBIT “.cub” file, VTK, etc.) as well as mesh written into MOAB directly by various software libraries (e.g. MeshKit).  Although there is no standard for naming or storing meta-data with a mesh, there is a great deal of commonality in the types of meta-data typically found with mesh data.  This document describes conventions that have been established for commonly encountered meta-data.  Various mesh readers implemented in MOAB attempt to read meta-data from a file and write it into the MOAB data model using these conventions.  Although there is no requirement to store a given type of meta-data in the form described here, a number of services have been written to handle meta-data using these conventions, no matter the source of the meta-data being processed.

Several specific tools are often used in concert with MOAB and bear special mention here.  The CUBIT toolkit generates finite element meshes, and saves them to a native save file (referred to as a “.cub” file) which MOAB is able to read.  Reading CUBIT meshes into MOAB through the .cub file format is preferred over other formats, since most other mesh formats written by CUBIT do not save most meta-data.  The MeshKit library also generates mesh using CGM and MOAB, and uses the same conventions for storing meshes into MOAB.  Finally, MOAB includes a CGM reader which can read a geometric model into a faceted representation in MOAB.  Meta-data from all these tools are stored in MOAB using the conventions described here.

The MOAB data model consists of the following basic types:
- <B>Entity</B>: The basic elements of topology, e.g. vertex, edge, triangle, tetrahedron, etc.  MOAB represents all types in the finite element zoo, plus polygons and polyhedra.
- <B>Entity %Set</B>: An arbitrary collection of entities and other sets.  Sets can have parent/child relations with other sets, and these relations are distinct from “contains” relations.
- <B>Interface</B>: The interface object through which other entities are accessed, in the sense of object-oriented-programming.  iMesh refers to the interface as the “root” set.
- <B>Tag</B>: A piece of data that can be assigned a distinct value to each entity and entity set, and to the interface itself.  Tags have a prescribed name, size in bytes, and data type; allowed data types are integer, double, entity handle, and byte or opaque.
.

The following section describes each meta-data tag convention in detail; these conventions are also summarized in Table 1.

\ref md-contents "Top"

  \section meta-conventions  Meta-Data Conventions

Meta-data is stored in MOAB and iMesh in the form of tags applied to either entities or entity sets.  For meta-data represented as entity sets, the contents of those sets are determined by the convention, with tags on those sets identifying them with the convention and adding any other semantic data.

Each meta-data convention is described in a subsection below.  Each convention begins with a short description of:

- Whether tags associated with the convention are assigned to entities or entity sets
- The tag(s) associated with the convention; information for each tag includes the name, the data type (I=integer, D=double, C=character, H=handle), and the tag length.  Tag lengths are specified after an asterisk (*); for example, C*32 implies a tag with character type and length 32.  Unspecified lengths correspond to length one.
.

<H3>Name</H3>
(Data: Entity sets, entities; Tag(s): NAME/C*32)

Character strings are used in many different contexts in applications.  MOAB uses the “NAME” tag to store character strings used to name entities.  This tag is of byte-type and is of length 32 bytes.  Note that the string stored in this tag may or may not be terminated with a NULL character.  It is always prudent account for missing NULL terminator, to avoid buffer overflow errors in the application.  Applications are free to define their own version of the NAME tag with a longer length, though this definition may conflict with other services attempting to use this tag with the conventional size.  Applications needing a string tag with a longer or variable length can also use MOAB’s variable-length tag type, though this will not be compatible with iMesh.

<H3>Title </H3>
(Data: Entity sets (file or instance); Tag(s): TITLE/C*strlen)

The title tag is meant to hold the overall identifier of a mesh, written at generation time or read from a file generated with a non-MOAB tool.  The tag length is variable, and is set by the application directly (by calling the tag_create function) or indirectly (by embedding the title in a file read by MOAB).

<H3> Global Identifier </H3>
(Data: Entity sets, entities; Tag(s): GLOBAL_ID/I)

Global identifiers are used in many different contexts in applications.  Geometric model entities are identified by dimension and id, e.g. “Volume 1”.  Mesh vertices and elements are identified similarly in mesh generation codes.  Boundary conditions and material types are identified similarly.  This tag is used to store such information.  This tag is currently stored in a 32-byte integer, though this may change in the future.

<H3> Geometric Model Information </H3>
(Data: Entity sets; Tag(s): GEOM_DIMENSION/I, GLOBAL_ID/I, NAME/C*32, CATEGORY/C*32, GEOM_SENSE_2(EH[2]), GEOM_SENSE_N_ENTS(EH*N), GEOM_SENSE_N_SENSES(I*N))

Mesh generation is often performed starting from a geometric model, represented in some form of CAD engine.  Many of the meshes used by MOAB are generated based on the CGM library.  Geometric models contain both topological information (the topological entities in the geometric model) and shape information (the geometric shape of those entities), as well as other meta-data written to the entities in a model.  When a mesh is read from a CUBIT .cub file, meta-data from the geometric model is read and represented in the MOAB data model, as described below. <B> Note that although MOAB reads and represents meta-data associated with the geometric model, it does not represent the geometric model itself.</B>  Therefore, shape-related information, e.g. the arc length of an edge or surface normal at a given point, can be retrieved only from the model represented in CGM or another geometric modeling engine.

The information contained in a geometric model, read into and represented in MOAB, consists of:
- Model entities (vertex, edge, face, volume)
- Topological relationships between model entities
- Groups of model entities
- Model entity/group ids
- Model entity/group names
.

The storage of this information into MOAB's data model is described for each type is described below.

- <B>Entities </B>

 in the geometric model (VERTEX, EDGE, FACE, VOLUME) are each represented by an entity set<sup>1</sup>.  These sets are tagged with the “GEOM_DIMENSION” tag, with integer value equal to the topological dimension of the entity (VERTEX = 0, EDGE = 1, etc.)  These sets contain the mesh owned by the corresponding entity in the geometric model.  Note this does not include mesh owned by bounding entities; thus, the set for a FACE will not contain the mesh vertices owned by bounding EDGEs in the geometric model.  These sets may or may not contain mesh entities of intermediate dimension, e.g. mesh edges owned by a FACE or faces owned by a VOLUME, depending on the application generating the mesh or the file from which the mesh was read.  These sets are all set-types, i.e. the order of entities in the sets is not significant, except in the case of EDGE sets, where order of the mesh vertices and edges corresponds to the relative order of vertices and edges at the time of mesh generation.  In MOAB, these sets are non-tracking by default, i.e. entities do not have knowledge of which geometry sets they are members of.

<sup>1</sup>Body-type entities from CUBIT are not explicitly represented in MOAB.

- <B> Topological Relationships </B>

In the geometric model, each FACE is bounded by zero or more EDGEs; other topological relationships between geometric entities exist in a similar manner.  These relationships are embedded in the data model using parent/child relations between entity sets.  For example, the entity set corresponding to a FACE will have child sets, each corresponding to a bounding EDGE, and parent sets, each corresponding to a VOLUME bounded by that FACE.  The relative order of sets in those parent/child lists is not significant, thus, “loops” bounding a FACE cannot reliably be inferred from this data.

- <B> Groups </B>

Geometric entities are sometimes assigned to application-specific groups.  These groups are represented using entity sets, tagged with a “GROUP” tag whose value equals the group id.  Group sets are “set”-type, and are not tracking sets.  These sets contain the sets corresponding to geometric entities contained in the groups in the geometric model, as well as any mesh entities assigned to the group.

- <B> Sense </B>

A geometric face has a natural orientation, indicated by the direction of the normal to the face; similarly, edges have a natural orientation determined by the direction of the tangent.  When faces bound regions, or edges bound faces, they do so with a sense; if a region includes a face with forward sense, that means the face's natural normal direction points out of the volume.  If a face includes an edge with forward sense, that means that if one moves along the edge in the direction of its tangent, the material of the face is on the left hand side.  The sense of a face (edge) with respect to a region (face) it bounds is stored using tags on the face (edge).

Most models allow a face to be part of only two regions.  Therefore, to store the sense of a face with respect to regions including it, a tag with two values is used.  This tag is named GEOM_SENSE_2, and has 2 EntityHandle values.  The first value corresponds to the entity set for the region for which that face has a forward sense, and the second to the region for which that face has a reverse sense.

Edges can bound more than two faces.  Therefore, two variable-length tags are used, one to store the EntityHandles of the faces the edge bounds, and the other to store the sense with which the edge bounds the corresponding face.  These tags are named GEOM_SENSE_N_ENTS and GEOM_SENSE_N_SENSES, respectively.  These are stored as variable-length tags; see the MOAB user's guide for information on how to work with tags of this type.

The following sense values are used:
- 0: forward
- 1: reverse
- -1: unnknown
.

<H3> Material Type </H3>
(Data: Entity sets; Tag(s): MATERIAL_SET/I)

Most finite element and other PDE-based analysis codes require a material type for each cell or element in the simulation.  MOAB uses entity sets to store this information, in the form of entity sets.  The MATERIAL_SET tag is used to identify these sets.  The value of this tag is conventionally an integer; in most cases this stores a user-assigned identifier associated with that material.

CUBIT assigns material types using what it calls “element blocks”, with each element block given a user-assigned id number and optionally a name.  The CUBIT and Exodus file readers in MOAB read element blocks into MATERIAL_SET sets.

In CUBIT, materials are typically assigned by assigning geometric volumes to element blocks.  Therefore, material sets often contain entity sets corresponding to those volumes.  Thus, a materrial set in MOAB is unlikely to contain mesh entities directly; rather, that set contains other sets which contain mesh entities.  In these cases, mesh entities can be retrieved by passing a “recursive” flag to the appropriate function (MOAB), or by calling the getEntitiesRec extension function (iMesh) provided by MOAB.

<H3> Boundary Conditions (Dirichlet, Neumann)</H3>
Data: Entity sets; Tag(s): DIRICHLET_SET/I, NEUMANN_SET/I)

Boundary conditions are often specified in terms of geometric model entities, similar to material types.  MOAB uses entity sets to store this information as well.  The DIRICHLET_SET and NEUMANN_SET tags are used to represent Dirichlet- and Neumann-type boundary condition sets, resp.  By convention, Neumann sets usually contain (indirectly) intermediate-dimension entities like edges in a 2D mesh or faces in a 3D mesh, while Dirichlet sets usually contain vertices.  In addition, Neumann sets are represented as sets of faces, rather than as sides of elements.  Faces can be ordered “forward” or “reverse” with respect to one of the bounding elements, depending on whether the right-hand normal points into or out of the element.  Forward-sense faces are added to the Neumann set.  Reverse-sense faces are put into a separate set; that set is tagged with the NEUSET_SENSE tag, with value = -1; and that reverse set is added to the Neummann set.

<H3> Parallel Mesh Constructs </H3>
(Data: Entity sets, entities; Tag(s): PARALLEL_PARTITION/I, PSTATUS/C*1, PARALLEL_SHARED_PROC/I, PARALLEL/SHARED_HANDLE/H, PARALLEL_SHARED_PROCS/I*NP, PARALLEL_SHARED_HANDLES/H*NP)

On a parallel computer, MOAB can represent the mesh on each processor as well as information about entities shared with neighboring processors.  Some of this information is also relevant even when the mesh is represented on a serial machine.  MOAB uses several tag and set conventions to describe the parallel nature of a mesh.  This information is summarized here; for a more complete description of MOAB’s parallel mesh representation and functionality, see [3].

- <B> Parallel partition, parts </B>

Most parallel mesh applications use a domain decomposition approach, where each processor solves for a subset of the domain.  The set of entities solved by a given processor is referred to as a part, and the collection of parts together is called the partition.  MOAB stores each part in an entity set, marked with the PARALLEL_PARTITION tag, whose value is the rank of the processor assigned that part; an entity set which contains all part sets is given the PARALLEL_PARTITIONING_TAG_NAME tag, whose value is currently meaningless.  The MBZoltan tool included as a tool in MOAB can partition a mesh for parallel solution, and writes the partition to the mesh in the form of parts and partitions.  Both these types of sets can be accessed in a serial mesh, e.g. for visualization.

- <B> Part interfaces </B>

When a partitioned mesh has been loaded on a parallel computer, the part on a given processor may share portions of its boundary with parts on other processors.  These shared regions are called part interfaces, and are also represented using entity sets.  These sets are marked with the PARALLEL_INTERFACE tag, whose value is currently meaningless.

- <B> Shared processor and handle </B>

For entities shared between processors, it is helpful to know locally which other processor shares an entity, and what the entity’s handle is on the remote processor.  There are two cases which are useful to distinguish, first where an entity is shared with only one other processor (referred to as shared), and second when a processor is shared by more than one other processor (referred to as multi-shared).   Shared entities are given the PARALLEL_SHARED_PROC and PARALLEL_SHARED_HANDLE tags, which store the rank of the sharing processor and the handle of the entity on that processor, respectively.  Multi-shared entities are marked with the PARALLEL_SHARED_PROCS and PARALLEL_SHARED_HANDLES tags; these tags have a length NP assigned at compile time in MOAB, with default values of -1 for processor rank and zero for handle (which are each invalid values for the corresponding data).  The processors/handles sharing a given entity are then written on the front of the arrays.  So, for example, an entity on processor rank 0, shared by processors 1 and 2, would have a PARALLEL_SHARED_PROCS tag whose values would be [1, 2, -1, -1, …], with PARALLEL_SHARED_HANDLES values of [m, n, 0, 0, …], where m and n would be the handles of that entity on processors 1 and 2.  The shared versions of these tags are “dense”, with default values which denote unshared entities.  The multi-shared tags are sparse tags in MOAB, with no default value.

- <B> Parallel status </B>

In addition to the tags above, MOAB also defines the PSTATUS tag, whose bits contain information about the parallel status of a given entity.  Starting with least significant bit, these bits represent whether an entity is 1) not owned, 2) shared, 3) multi-shared, 4) interface, 5) a ghost entity.  The first bit being set indicates “not owned” so that the default value for this tag, of zero, corresponds to an owned, unshared entity, which will be the state of most entities on a given processor.

<H3>Structured Mesh Parameters </H3>

MOAB has a structured mesh interface for creating structured mesh (see “ScdInterface.hpp” header file in MOAB source code).  Along with an internal representation that is more memory-efficient (since it does not need to store connectivity), MOAB also creates and tags entity sets with structured mesh parameters, which can be accessed through the normal tag and set interface.  The following tags are used:

- <B>BOX_DIMS</B>: This tag stores the ijk coordinates of the lower and upper corner of the structured mesh box(es).
- <B>GLOBAL_BOX_DIMS</B>: If specified when the structured mesh is created, a tag with this name stores the global box dimensions (which may be different than the local box dimensions).
- <B>BOX_PERIODIC</B>: Stores whether the box is periodic in the i (BOX_PERIODIC[0]) and j (BOX_PERIODIC[1]) directions.
- <B>__BOX_SET</B>: Pointer to the ScdBox instance corresponding to this entity set.<sup>2</sup>
.
Although the structured mesh is not saved as such in HDF5-format files, the entity sets and corresponding tags will be saved and restored.

<sup>2</sup>The double-underscore in the tag name implies that this tag will not be saved in a file, in this case because the ScdBox instances are not preserved in a file.

<H3>Spectral Mesh Constructs </H3>

The Spectral Element Method (SEM) is a high-order method, using a polynomial Legendre interpolation basis with Gauss-Lobatto quadrature points, in contrast to the Lagrange basis used in (linear) finite elements.  A spectral mesh with order O contains quadrilateral or hexahedral elements comprised of (O+1)d vertices.  Spectral meshes are usually represented in one of two ways, either as coarse elements which point to an array of higher-order vertices (and with corner vertices represented in the normal manner), or as linear quads/hexes formed from the higher-order vertices, with each original coarse quad/hex represented by Od fine quads/hexes.  Similarly, the spectral variables, which are normally computed at fine vertex positions, are stored either on those vertices, or in lexicographically-ordered arrays on elements (with tag values repeated on neighboring elements).  MOAB can read spectral meshes from a variety of formats (at this time, including CAM-SE, HOMME, and Nek5000).  Which of the above two representations are controlled by read options and are indicated by certain tags:

- SPECTRAL_MESH: read option indicating that spectral elements should be represented as coarse linear quads/hexes and each element containing an array of lexicographically-ordered vertex handles

- TAG_SPECTRAL_ELEMENTS: read option; if given, spectral variables are represented as lexicographically-ordered arrays on elements

- TAG_SPECTRAL_VERTICES: read option; if given, spectral variables are represented as tags on vertices

- CONN=<filename>: in CAM-SE, the connectivity of the spectral mesh is stored by default in a file named “HommeMapping.nc”; this option can be given to read the connectivity from a different file

- SPECTRAL_VERTICES: tag name for array of vertex handles

- SPECTRAL_ORDER: tag name for spectral order, written to file set or (if no file set given) to interface after a spectral mesh is read

.

\ref md-contents "Top"

  \section meta-options Reader/Writer Options

All mesh file readers and writers in MOAB take an option string as an argument.  By default, the semicolon (“;”) delimits individual options in the option string.  Options used in multiple readers are described in this section; the options enabled in specific readers/writers are described in the corresponding appendix at the end of this document.

<H3>variable=\<var_name\>[,...]</H3>

By default, all field data stored with the mesh is read with the mesh, and stored as tags on the associated mesh entities.  This option lists specific variables that should be read along with the mesh (note also the “nomesh” option, described elsewhere in this document).  The variable name listed will be read into a tag with the same name.  For time-dependent variables, the time step number will be appended to the variable name to form the tag name.  If no “timestep” or “timeval” option is given, all time steps will be read, resulting in several tags being created.  If the “nomesh” option is given, the application must pass the entity set resulting from the original mesh read in to the function, that this set must contain the mesh read only from that file.  The mesh in the file is checked against the mesh in the set to verify that the two correspond.  The special name “MOAB_ALL_VARIABLES” can be used to indicate that all variables should be read.  Multiple variable names can be specified, separated from each other by commas.

<H3>nomesh </H3>

Indicates that no mesh should be read from the file.  This option is used in conjunction with the “variable=” option, to read variables and assign them as tags to a previously-read mesh.  If this option is used, applications should pass an entity set to the read function, which should contain the mesh previously read from the file.

<H3>timestep=\<step_number\>[, ...] </H3>

Read the time step number whose time value is equal to or greater than the specified time value, for the specified variable(s).  Tag names for the variable(s) will be formed by appending the time step number to the variable name.  Multiple time step values can be specified, separated from each other by commas.

<H3>timeval=\<time_value\>[, ...]</H3>

Read the time step number whose time value is equal to or greater than the specified time value, for the
specified variable(s). Tag names for the variable(s) will be formed by appending the time step number
to the variable name. Multiple time step values can be specified, separated from each other by commas.

<H3>gather_set[=\<rank\>] </H3>

Create a gather set (associated with tag GATHER_SET) on one processor with the specified rank, to duplicate entities on other processors. If the rank is not specified, it will be rank 0 by default. If an invalid rank is passed, no gather set will be created. Gather set is specially used by HOMME, MPAS, and any other unstructured grid.

<H3>no_mixed_elements </H3>

Indicates that no mixed elements (e.g. pentagons and hexagons) should be created by the MPAS reader. If this option is used, a common parameter maxEdgesPerCell will be computed to be used across all processors (instead of the one reported in the MPAS file header, which is usually 10), and each cell is created with maxEdgesPerCell edges. Any cell that has less actual edges will be padded by duplicating the last vertex in the connectivity array. As a result, all created cells will be in one contiguous chunk.

<H3>no_edges </H3>

Indicates that no edges should be created and no edge variables will be read. This option can be used when there is no need to read variables on edges. For a huge MPAS file with 65M cells, it can save more than 3GB MOAB internal storage for edge connectivity.

\ref md-contents "Top"

  \section meta-references References

[1]     T.J. Tautges, R. Meyers, K. Merkley, C. Stimpson, and C. Ernst, MOAB: A Mesh-Oriented Database, Sandia National Laboratories, 2004.

[2]     L. Diachin, A. Bauer, B. Fix, J. Kraftcheck, K. Jansen, X. Luo, M. Miller, C. Ollivier-Gooch, M.S. Shephard, T. Tautges, and H. Trease, “Interoperable mesh and geometry tools for advanced petascale simulations,” Journal of Physics: Conference Series,  vol. 78, 2007, p. 012015.
[3]     T.J. Tautges, J.A. Kraftcheck, N. Bertram, V. Sachdeva, and J. Magerlein,  "Mesh Interface Resolution and Ghost Exchange in a Parallel Mesh Representation", In Proceedings of the 2012 IEEE 26th International Parallel and Distributed Processing Symposium Workshops & PhD Forum (IPDPSW '12), 2012.

\ref md-contents "Top"

  \section appendixA Appendix A: Summary

\subsection table1 Table 1: Summary of MOAB meta-data conventions.

<table border="1">
<tr>
<th>Convention</th>
<th>Applies to (E=ent, S=set)</th>
<th>Tag(s) (type/length)</th>
<th>Description</th>
</tr>
<tr>
<td>Name</td>
<td>E, S</td>
<td>NAME/C*32</td>
<td></td>
</tr>
<tr>
<td>Title</td>
<td>S</td>
<td>TITLE/C*strlen</td>
<td>Title of mesh</td>
</tr>
<tr>
<td>Global identifier</td>
<td>E, S</td>
<td>GLOBAL_ID/I</td>
<td></td>
</tr>
<tr>
<td>Geometric topology</td>
<td>S</td>
<td>GEOM_DIMENSION/I, GLOBAL_ID/I,^
NAME/C*32,
CATEGORY/C*32.
GEOM_SENSE_2/EH[2],
GEOM_SENSE_N_ENTS/EH*N,
GEOM_SENSE_N_SENSES/I*N</td>
<td>%Sets contain mesh owned by that entity; parent/child links to bounded/bounding entities in geometric model</td>
</tr>
<tr>
<td>Material type</td>
<td>S</td>
<td>MATERIAL_SET/I</td>
<td>%Set contains entities or sets assigned a common material type</td>
</tr>
<tr>
<td>Boundary condition</td>
<td>S</td>
<td>DIRICHLET_SET/I, NEUMANN_SET/I</td>
<td>%Set contains entities or sets assigned a particular boundary condition; neumann sets usually contain edges (2D) or faces (3D)</td>
</tr>
<tr>
<td>Parallel mesh constructs</td>
<td>E, S</td>
<td>PARALLEL_MESH_PARTITIONING/I, PARALLEL_PARTITION/I, PSTATUS/C*1, PARALLEL_SHARED_PROC/I, PARALLEL/SHARED_HANDLE/H, PARALLEL_SHARED_PROCS/I*NP, PARALLEL_SHARED_HANDLES/H*NP</td>
<td> Data which describes parallel mesh</td>
</tr>
<tr>
<td>Structured mesh constructs</td>
<td>S</td>
<td>BOX_DIMS/I*6, GLOBAL_BOX_DIMS/I*6, BOX_PERIODIC/2*I, __BOX_SET/O</td>
<td>Data describing structured mesh </td>
</tr>
<tr>
<td>Spectral mesh constructs </td>
<td>E, S</td>
<td>SPECTRAL_ORDER/I, SPECTRAL_VERTICES/I*(O+1)^2</td>
<td>Data marking spectral mesh constructs</td>
</tr>
</table>
 
  \ref meta-introduction "Back to Introduction"

  \subsection table2 Table 2: Summary of MOAB conventional tag names, types, and purposes.  Data types are I=integer, D=double, C=character, H=entity handle,O=opaque.  Data type with *x denote length of x elements of that data type.

<Table border="1">
<tr>
<th>Tag name</th>
<th>Data type</th>
<th>Applies to (E=entity, S=set)</th>
<th>Purpose</th>
</tr>
<tr>
<td>BOX_DIMS</td>
<td>I*6</td>
<td>S</td>
<td>Lower and upper ijk dimensions of box, ordered (ilo, jlo, klo, ihi, jhi, khi)</td>
</tr>
<tr>
<td>BOX_PERIODIC</td>
<td>I*2</td>
<td>S</td>
<td>Indicates whether box is periodic in i (BOX_PERIODIC[0]) or j (BOX_PERIODIC[1])</td>
</tr>
<tr>
<td>__BOX_SET</td>
<td>O</td>
<td>S</td>
<td>Pointer to corresponding ScdBox instance</td>
</tr>
<tr>
<td>CATEGORY</td>
<td>C*32</td>
<td>S</td>
<td>String describing purpose of set; examples include “group”, “vertex”, “edge”, “surface”, “volume”</td>
</tr>
<tr>
<td>DIRICHLET_SET </td>
<td>I</td>
<td>SO</td>
<td>Entities or sets with common boundary condition</td>
</tr>
<tr>
<td>GEOM_DIMENSION</td>
<td>I</td>
<td>S</td>
<td>Identifies mesh entities resolving a given geometric model entity</td>
</tr>
<tr>
<td>GEOM_SENSE_2</td>
<td>EH*2</td>
<td>S</td>
<td> Stored on face-type geometric topology sets, values store regions having forward and reverse sense</td>
</tr>
<tr>
<td>GEOM_SENSE_N_ENTS</td>
<td>EH*N</td>
<td>S</td>
<td>Stored on edge-type geometric topology sets, values store faces whose senses are stored in GEOM_SENSE_N_SENSES.</td>
</tr>
<tr>
<td>GEOM_SENSE_N_SENSES</td>
<td>I*N</td>
<td>S</td>
<td>Stored on edge-type geometric topology sets, values store senses of the edge with respect to faces stored in GEOM_SENSE_N_ENTS.</td>
</tr>
<tr>
<td>GLOBAL_ID</td>
<td>I</td>
<td>E,S</td>
<td>Application-specific entity id</td>
</tr>
<tr>
<td>MATERIAL_SET</td>
<td>I</td>
<td>S</td>
<td>Entities or sets grouped by material type</td>
</tr>
<tr>
<td>NAME</td>
<td>C*32</td>
<td>E, S</td>
<td>User-assigned entity name(s); multiple names delimited with ?</td>
</tr>
<tr>
<td>NEUMANN_SET</td>
<td>I</td>
<td>S</td>
<td>Entities or sets with common boundary condition </td>
</tr>
<tr>
<td>PARALLEL_PARTITION </td>
<td>I</td>
<td>S</td>
<td>Represent a part in a partition</td>
</tr>
<tr>
<td>PARALLEL_MESH_PARTITIONING</td>
<td>I</td>
<td>S</td>
<td>Represents a partition of the mesh for parallel solution, which is a collection of parts</td>
</tr>
<tr>
<td>__PARALLEL_SHARED_HANDLE</td>
<td>H</td>
<td>E, S</td>
<td> Handle of this entity/set on sharing processor</td>
</tr>
<tr>
<td>__PARALLEL_SHARED_PROC</td>
<td>I</td>
<td>E,S</td>
<td>Rank of other processor sharing this entity/set </td>
</tr>
<tr>
<td>__PARALLEL_SHARED_HANDLES</td>
<td>H*NP</td>
<td>E,S</td>
<td>Handles of this entity/set on sharing processors </td>
</tr>
<tr>
<td>__PARALLEL_SHARED_PROCS</td>
<td>I*NP</td>
<td>E,S</td>
<td>Ranks of other processors sharing this entity/set </td>
</tr>
<tr>
<td>__PARALLEL_STATUS</td>
<td>C*1</td>
<td>E,S</td>
<td>Bit-field indicating various parallel information </td>
</tr>
<tr>
<td>SPECTRAL_ORDER</td>
<td>I</td>
<td>S</td>
<td> Order of a spectral mesh </td>
</tr>
<tr>
<td>SPECTRAL_VERTICES</td>
<td>H*(O+1)^d</td>
<td>E</td>
<td> Vertices comprising a spectral element, ordered lexicographically; here, O=value of SPECTRAL_ORDER tag. </td>
</tr>
</table>

\ref md-contents "Top"

  \section appendixB Appendix B: CCMIO (Star-CD, Star-CCM+) Reader/Writer Conventions

  \subsection table3 Table 3: Translation between CCMIO options and MOAB tags.
<Table border="1">
<tr>
<th> %Set Type</th>
<th>CCMIO Construct</th>
<th>MOAB Tag Name, Type</th>
</tr>
<tr>
<td rowspan="2">File set / Interface</td>
<td>Title (option)</td>
<td>“Title” (C*32)</td>
</tr>
<tr>
<td>CreatingProgram</td>
<td>“CreatingProgram” (C*32)</td>
</tr>
<tr>
<td rowspan="13">Material sets</td>
<td>Index</td>
<td>MATERIAL_SET</td>
</tr>
<tr>
<td>Label<sup>1</sup></td>
<td>NAME</td>
</tr>
<tr>
<td>MaterialId</td>
<td>“MaterialId” (I)</td>
</tr>
<tr>
<td>Radiation</td>
<td>“Radiation” (I)</td>
</tr>
<tr>
<td>PorosityId</td>
<td>“PorosityId” (I)</td>
</tr>
<tr>
<td>SpinId</td>
<td>“SpinId” (I)</td>
</tr>
<tr>
<td>GroupId</td>
<td>“GroupId” (I)</td>
</tr>
<tr>
<td>ColorIdx</td>
<td>“ColorIdx” (I)</td>
</tr>
<tr>
<td>ProcessorId</td>
<td>“ProcessorId” (I)</td>
</tr>
<tr>
<td>LightMaterial</td>
<td>“LightMaterial” (I)</td>
</tr>
<tr>
<td>FreeSurfaceMaterial</td>
<td>“Thickness” (F)</td>
</tr>
<tr>
<td>Thickness</td>
<td>“Thickness” (F)</td>
</tr>
<tr>
<td>MaterialType</td>
<td>“MaterialType” (C*32)</td>
</tr>
<tr>
<td rowspan="5">Neumann sets</td>
<td>Index</td>
<td>NEUMANN_SET</td>
</tr>
<tr>
<td>Label</td>
<td>NEUMANN_SET</td>
</tr>
<tr>
<td>BoundaryName</td>
<td>NAME</td>
</tr>
<tr>
<td>BoundaryType</td>
<td>“BoundaryType” (C*32)</td>
</tr>
<tr>
<td>ProstarRegionNumber</td>
<td>“ProstarRegionNumber” (I)</td>
</tr>
</table>

Notes:
1. If no name is present, labels the material group with “MaterialX”, where X is the index of that group.

\ref md-contents "Top"

  \section appendixC Appendix C: ExodusII Reader/Writer Conventions 

  \subsection table4 Table 4: Translation between ExodusII constructs and MOAB tags.
<Table border="1">
<tr>
<th> Data Type</th>
<th>ExodusII Construct</th>
<th>MOAB Tag Name, Type</th>
</tr>
<tr>
<td></td>
<td>QA records</td>
<td>“qaRecord” (C*(v))<sup>2</sup></td>
</tr>
<tr>
<td rowspan="2">Material sets</td>
<td>Block number</td>
<td>MATERIAL_SET</td>
</tr>
<tr>
<td>Block element type</td>
<td>Entity type, # vertices per entity</td>
</tr>
<tr>
<td rowspan="2">Dirichlet sets<sup>3</sup></td>
<td>Nodeset number</td>
<td>DIRICHLET_SET</td>
</tr>
<tr>
<td>Distribution factors</td>
<td>“distFactor” (D*(v))<sup>1</sup></td>
</tr>
<tr>
<td>Neumann sets</td>
<td>Sideset number</td>
<td>NEUMANN_SET</td>
</tr>
<tr>
<td rowspan="2">Neumann sets, reverse faces3<sup>3</sup></td>
<td>Distribution factors</td>
<td>“distFactor” (D*(v))<sup>1</sup></td>
</tr>
<tr>
<td>Sides</td>
<td>SENSE</td>
</tr>
<tr>
<td>Nodes, elements</td>
<td>node_num_map, elem_map</td>
<td>GLOBAL_ID on nodes/elements</td>
</tr>
</table>

Notes:
-# Variable-length tag used for distribution factors; length for each set is the number of entities in
each set, such that there is one distribution factor for each entity in the set.
-# QA records are stored as variable-length tags on file set specified on read. Tag is a
concatenation of QA record strings into a single string, with '\0' used to delimit lines.
-# MOAB represents sidesets as sets of faces, rather than as sides of elements. Faces can be
ordered “forward” or “reverse” with respect to one of the bounding elements, depending on
whether the right-hand normal points into or out of the element. Forward-sense faces are added
to the Neumann set. Reverse-sense faces are put into a separate set; that set is tagged with the SENSE tag, with value = -1; and that reverse set is added to the Neummann set.
.

  \ref md-contents "Top"

  \section appendixD Appendix D: NC (Climate Data) Reader/Writer Conventions

The climate data reader in MOAB reads files with the '.nc' filename extension. By default, this reader
reads the whole mesh in the file and creates it as structured mesh in MOAB, with the mesh accessible
through MOAB's structured mesh interface. By default, all variables and timesteps are read from the
file, and written as tags on the mesh vertices from that file. This behavior is controlled by the
“variable”, “nomesh”, “timestep”, and “timeval” options described earlier in this document. If MOAB
is compiled for parallel execution and configured with a pnetcdf reader, the mesh is read in parallel,
with a 1D or 2D decomposition designed to balance read performance and communication interface
size (for details on the partitioning method used, see the src/io/ReadNC.cpp source file).

Mesh is put into the entity set provided to the load_file function. This entity set is also annotated with
various tags representing information read from the file. These tags are described in Table 5.

Reading unstructured NC files in the HOMME format is also supported. Currently a trivial
element-based partition is the only option for parallel reading. As the data is unstructured, it is necessary to have a connectivity file to define the vertex adjacencies. The default convention is to have a file called HommeMapping.nc in the same directory as the the variable data file. If this convention is not followed, the connectivity file can be specified with the option -O CONN=”/path/to/connectivity.nc”. An example of mbconvert using the parallel read capability is shown below:

<B>  mpiexec -np 2 tools/mbconvert -O TRIVIAL -O DEBUG_IO=1 -o DEBUG_IO=9 -o PARALLEL=WRITE_PART /nfs2/hayes6/meshlab/homme_data/camrun.cam2.h0.0000-01-01-16200.nc output.h5m </B>

Several other things to note about reading climate data files into MOAB:
- Time-dependent variables: MOAB currently has no mechanism for time-dependent tags. Therefore, time-dependent variables are represented using one tag per timestep, with the tag name set as the variable name plus the timestep index. Thus, the first few timesteps for the variable TEMPERATURE would be represented in tags named TEMPERATURE0, TEMPERATURE1, etc.
- Cell- and face-centered variables: The climate data reader currently does not do cell- and face-
centered variables correctly.
.
  \subsection table5 Table 5: Summary of MOAB conventional tag names, types, and purposes. Data types are I=integer, D=double, C=character, H=entity handle. Data type with *x denote length of x elements of that data type; data type with *var denote variable-length tag. Tag names with two underscores prepended (“__”) denote tags not written to a file by MOAB.

<Table border="1">
<tr>
<th> Tag name </th>
<th>Data type </th>
<th> Applies to (E=entity, S=set) </th>
<th>Purpose </th>
</tr>
<tr>
<td>__NUM_DIMS </td>
<td>I</td>
<td>S</td>
<td>The number of dimensions in the netcdf file.</td>
</tr>
<tr>
<td>__NUM_VARS</td> 
<td>I</td>
<td>S</td>
<td>The number of variables in the netcdf file.</td>
</tr>
<tr>
<td>__DIM_NAMES </td>
<td>C*var</td>
<td>S</td>
<td>The dimension names, concatenated into a
character string, with '\0' terminating each name.
</td>
</tr>
<tr>
<td>__DIM_LENS </td>
<td>I*var</td>
<td>S</td>
<td>A vector of integers, storing the length of
each dimension.
</td>
</tr>
<tr>
<td>__VAR_NAMES
</td>
<td>C*var</td>
<td>S</td>
<td>The variable names, concatenated into a
character string, with '\0' terminating each name.
</td>
</tr>
<tr>
<td><dim_name> 
</td>
<td>(I or D)*var</td>
<td>S</td>
<td>For each dimension, the values for the dimension.
The data type for this tag corresponds to that in the
netcdf file. The length of this tag is the number of
values stored for the dimension in the netcdf file.</td>
</tr>
<tr>
<td>__<dim_name>_LOC_MIN_MAX</td>
<td>(I or D)*2</td>
<td>S</td>
<td>The indices (0-based) of the local min and max
values of dimension stored locally. For spatial
dimensions like lon or lat, this will store the
minimum and maximum indices in the local partition
of the grid. For dimensions like time, where each
processor represents the entire dimension, this will
likely store 0 and the number of values for that
dimension. Only one of __<dim_name>_LOC_VALS and
__<dim_name>_LOC_MIN_MAX can be used for a given
dimension.</td>
</tr>
<tr>
<td>__<dim_name>_LOC_VAL </td>
<td>(I or D)*var</td>
<td>S</td>
<td>The indices (0-based) of the dimension stored
locally. This tag only makes sense for dimensions
that can be read in multiple pieces, such as time.
Only one of __<dim_name>_LOC_VALS and
__<dim_name>_LOC_MIN_MAX can be used for a given
dimension.</td>
</tr>
<tr>
<td>__<dim_name>_GLOBAL_MIN_MAX</td>
<td>(I or D)*2</td>
<td>S</td>
<td>The indices (0-based) of the global min and max
values of dimension.</td>
</tr>
<tr>
<td>__<var_name>_DIMS 
</td>
<td>H*n 
</td>
<td>S</td>
<td>For each variable, this tag stores the tag
handles for the n dimensions defining this variable,
in netcdf ordering (last dimension varying fastest).
The size of this tag is n * sizeof(TagHandle).
</td>
</tr>
<tr>
<td><var_name><timestep_ind> 
</td>
<td>(data type)</td>
<td>E</td>
<td>Values of the variable for timestep <timestep_ind>
for vertices. The data type of this tag corresponds
to that of the variable from the netcdf file.
Timestep index is 0-based.
</td>
</tr>
<tr>
<td>__GLOBAL_ATTRIBS 
</td>
<td>C*var
</td>
<td>S</td>
<td>The global attributes, concatenated into a character
string, with ‘\0’ terminating each attribute name, ‘;’
       separating the data type and value, and ‘;’
          separating one name/data type/value from the next.
</td>
</tr>
<tr>
<td>__GLOBAL_ATTRIBS_LEN 
</td>
<td>I*var
</td>
<td>S</td>
<td>A vector of integers, marking the end position of
each attribute (name/data type/value) in __GLOBAL_ATTRIBS tag.
</td>
</tr>
<tr>
<td>__<var_name>_ATTRIBS 
</td>
<td>C*var
</td>
<td>S</td>
<td>The variable attributes, concatenated into a
character string, with ‘\0’ terminating each attribute
   name, ‘;’ separating the data type and value, and ‘;’
          separating one name/data type/value from the next.
</td>
</tr>
<tr>
<td>__<var_name>_ATTRIBS_LEN 
</td>
<td>I*var
</td>
<td>S</td>
<td>A vector of integers, marking the end position of
each attribute (name/data type/value) in
__<var_name>_ATTRIBS tags
</td>
</tr>
</table>

  \ref md-contents "Top"

  \section appendixE Appendix E: Nek5000 Reader/Writer Conventions

Nek5000, or Nek, is a code that uses the spectral element method to model fluid, heat transfer,
electromagnetics, and other physics. Nek uses unstructured hexahedral meshes, with each hex element
resolved by a structured grid of “Gauss Lebato Legendre” (GLL) points. Nek can read meshes through
MOAB, and can output physics variables and GLL points through MOAB as well.

Since fluid is a single material in Nek, no material sets are needed. Boundary conditions are mapped to
Nek's cbc array using Neumann sets and a user-provided “usr_moab2nek” subroutine (for an example
of this subroutine, see examples/moab/pipe.usr in the Nek source code). GLL point locations and fluid
variables on those points are stored in tags on the hex elements. All hex elements have the same
number of GLL points. The number of GLL points in each direction is stored in a tag on the mesh
instance. These tags are described in Table 6.

GLL point locations and fluid variables are stored in lexicographic order, similar to their storage order
inside the Nek code.

  \subsection table6 Table 6: Summary of MOAB conventional tag names, types, and purposes for Nek. Data types are I=integer, D=double, C=character, H=entity handle. Data type with *x denote length of x elements of that data type; data type with *var denote variable-length tag. Tag names with two underscores prepended (“__”) denote tags not written to a file by MOAB.
<Table border="1">
<tr>
<th> Tag name </th>
<th> Data Type</th>
<th>Applies to (E=entity, S=set)</th>
<th>Purpose</th>
</tr>
<tr>
<td>SEM_DIMS</td>
<td>I*3</td>
<td>S</td>
<td>The dimensions of the GLL mesh in each hex
element.
</td>
</tr>
<tr>
<td>SEM_X</td>
<td>D*nx*ny*nz</td>
<td>E</td>
<td>X position of GLL points (having nx*ny*nz
values)
</td>
</tr>
<tr>
<td>SEM_Y</td>
<td>D*nx*ny*nz</td>
<td>E</td>
<td>Y position of GLL points (having nx*ny*nz values)</td>
</tr>
<tr>
<td>SEM_Z</td>
<td>D*nx*ny*nz</td>
<td>E</td>
<td>Z position of GLL points (having nx*ny*nz values)</td>
</tr>
<tr>
<td>VEL_X</td>
<td>D*nx*ny*nz</td>
<td>E</td>
<td>Fluid velocities in the x direction for GLL point
array (having nx*ny*nz values)</td>
</tr>
<tr>
<td>VEL_Y</td>
<td>D*nx*ny*nz</td>
<td>E</td>
<td>Fluid velocities in the y direction for GLL point
array (having nx*ny*nz values)</td>
</tr>
<tr>
<td>VEL_Z</td>
<td>D*nx*ny*nz</td>
<td>E</td>
<td>Fluid velocities in the z direction for GLL point
array (having nx*ny*nz values)</td>
</tr>
<tr>
<td>TEMP</td>
<td>D*nx*ny*nz</td>
<td>E</td>
<td>Fluid temperature for GLL point array (having
nx*ny*nz values)
</td>
</tr>
<tr>
<td>PRESS</td>
<td>D*nx*ny*nz</td>
<td>E</td>
<td>Fluid pressure for GLL point array (having
nx*ny*nz values)
</td>
</tr>
</table>
  \ref md-contents "Top"
        */

/*!  \page md-tables List of Tables
    \ref table1

    \ref table2

    \ref table3

    \ref table4

    \ref table5

    \ref table6

*/

