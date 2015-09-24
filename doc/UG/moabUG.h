/*! \page userguide User's Guide (MOAB 4.8.2)
 
  \subpage team 
 
  \subpage contents
 
  \subpage figures
 
  \subpage tables
 
  \subpage differences

  \subpage building

  \page team MOAB team members
 <h2>The MOAB Team, including: </h2>
 
 - Vijay S. Mahadevan (Argonne National Lab)
 - Timothy J. Tautges (CD-Adapco, Univ Wisconsin-Madison)
 - Iulian Grindeanu (Argonne National Lab) 
 - Rajeev Jain (Argonne National Lab)
 - Danqing Wu  (Argonne National Lab)
 - Navamita Ray (Argonne National Lab)
 - Jane Hu (Univ Wisconsin-Madison)
 - Paul Wilson (Univ Wisconsin-Madison)
 - Patrick Shriwise (Univ Wisconsin-Madison)
 - Anthony Scopatz (Univ Wisconsin-Madison)


 <h2>Emeritus members:</h2>
 
 - Jason A. Kraftcheck
 - Brandon M. Smith
 - Hong-Jun Kim
 - Jim Porter
 - Xiabing Xu
 
  \page contents Table of Contents
 
  \ref introduction  

  \ref interface     

	\ref twoone    

	\ref twotwo     

	\ref twothree       

	\ref twofour   

  \ref api     

  \ref services      

    \ref fourone    

    \ref fourtwo   

    \ref fourthree  

    \ref fourfour      

    \ref fourfive    

    \ref foursix

    \ref fourseven

    \ref foureight

  \ref parallel      

    \ref fiveone    

    \ref fivetwo     

    \ref fivethree    

    \ref fivefour      

  \ref applications   

  \ref implementation         

  \ref representation     

  \ref element    

    \ref nineone  

    \ref ninetwo        

    \ref ninethree      

  \ref performance   

  \ref error-handling

  \ref conclusions    

  \ref references 

  \section introduction 1.Introduction

In scientific computing, systems of partial differential equations (PDEs) are solved on computers.  One of the most widely used methods to solve PDEs numerically is to solve over discrete neighborhoods or “elements” of the domain.  Popular discretization methods include Finite Difference (FD), Finite Element (FE), and Finite Volume (FV).  These methods require the decomposition of the domain into a discretized representation, which is referred to as a “mesh”.  The mesh is one of the fundamental types of data linking the various tools in the analysis process (mesh generation, analysis, visualization, etc.).  Thus, the representation of mesh data and operations on those data play a very important role in PDE-based simulations.
 
MOAB is a component for representing and evaluating mesh data.  MOAB can store structured and unstructured mesh, consisting of elements in the finite element “zoo”, along with polygons and polyhedra.  The functional interface to MOAB is simple, consisting of only four fundamental data types.  This data is quite powerful, allowing the representation of most types of metadata commonly found on the mesh.  Internally MOAB uses array-based storage for fine-grained data, which in many cases provides more efficient access, especially for large portions of mesh and associated data.  MOAB is optimized for efficiency in space and time, based on access to mesh in chunks rather than through individual entities, while also versatile enough to support individual entity access.

The MOAB data model consists of the following four fundamental types: mesh interface instance, mesh entities (vertex, edge, tri, etc.), sets, and tags.  Entities are addressed through handles rather than pointers, to allow the underlying representation of an entity to change without changing the handle to that entity.  Sets are arbitrary groupings of mesh entities and other sets.  Sets also support parent/child relationships as a relation distinct from sets containing other sets.  The directed graph provided by set parent/child relationships is useful for embedding graphs whose nodes include collections of mesh entities; this approach has been used to represent a wide variety of application-specific data, including geometric model topology, processor partitions, and various types of search trees.  Tags are named data which can be assigned to the mesh as a whole, individual entities, or sets.  Tags are a mechanism for attaching data to individual entities, and sets are a mechanism for describing relations between entities; the combination of these two mechanisms is a powerful yet simple interface for representing metadata or application-specific data.

Various mesh-related tools are provided with MOAB or can be used directly with MOAB.  These tools can be used for mesh format translation (mbconvert), mesh skinning (Skinner class), solution transfer between meshes (MBCoupler tool), ray tracing and other geometric searches (OrientedBoxTreeTool, AdaptiveKDTree), visualization (vtkMOABReader tool), and relation between mesh and geometric models (the separately-packed Lasso tool).  These tools are described later in this document.

MOAB is written in the C++ programming language, with applications interacting with MOAB mostly through its moab::Interface class.  All of the MOAB functions and classes are isolated in and accessed through the moab namespace<sup>1</sup>. The remainder of this report gives class and function names without the “moab::” namespace qualification; unless otherwise noted, the namespace qualifier should be added to all class and function names referenced here.  MOAB also implements the iMesh interface, which is specified in C but can be called directly from other languages.  Almost all of the functionality in MOAB can be accessed through the iMesh interface.  MOAB is developed and supported on the Linux and MacOS operating systems, as well as various HPC operating systems.  MOAB can be used on parallel computing systems as well, including both clusters and high-end parallel systems like IBM BG/P and Cray systems.  MOAB is released under a standard LGPL open source software license.

MOAB is used in several ways in various applications.  MOAB serves as the underlying mesh data representation in several scientific computing applications [1].  MOAB can also be used as a mesh format translator, using readers and writers included in MOAB.  MOAB has also been used as a bridge to couple results in multi-physics analysis and to link these applications with other mesh services [2].

The remainder of this report is organized as follows.  Section 2, “Getting Started”, provides a few simple examples of using MOAB to perform simple tasks on a mesh.  Section 3 discusses the MOAB data model in more detail, including some aspects of the implementation.  Section 4 summarizes the MOAB function API.  Section 5 describes some of the tools included with MOAB, and the implementation of mesh readers/writers for MOAB.  Section 6 describes how to build MOAB-based applications.  Section 7 contains a brief description of MOAB’s relation to the iMesh mesh interface.  Sections 8 and 9 discuss MOAB's representations of structured and spectral element meshes, respectively.  Section 10 gives helpful hints for accessing MOAB in an efficient manner from applications.  Section 11 gives a conclusion and future plans for MOAB development.  Section 12 gives references cited in this report.

Several other sources of information about MOAB may also be of interest to readers.  Meta-data conventions define how sets and /or tags are used together to represent various commonly-used simulation constructs; conventions used by MOAB are described in Ref [4], which is also included in the MOAB source distribution.  This document is maintained separately from this document, since it is expected to change over time.  The MOAB project maintains a wiki [5], which links to most MOAB-related information.  MOAB also uses several mailing lists [6],[7] for MOAB-related discussions.  Potential users are encouraged to interact with the MOAB team using these mailing lists.

<sup>1</sup> Non-namespaced names are also provided for backward compatibility, with the “MB” prefix added to the class or variable name.

 \ref contents

 \section interface 2.MOAB Data Model
The MOAB data model describes the basic types used in MOAB and the language used to communicate that data to applications.  This chapter describes that data model, along with some of the reasons for some of the design choices in MOAB.

 \ref contents

 \subsection twoone 2.1. MOAB Interface 

MOAB is written in C++.  The primary interface with applications is through member functions of the abstract base class Interface.  The MOAB library is created by instantiating Core, which implements the Interface API.  Multiple instances of MOAB can exist concurrently in the same application; mesh entities are not shared between these instances<sup>2</sup>.  MOAB is most easily viewed as a database of mesh objects accessed through the instance.  No other assumptions explicitly made about the nature of the mesh stored there; for example, there is no fundamental requirement that elements fill space or do not overlap each other geometrically.
 
<sup>2</sup> One exception to this statement is when the parallel interface to MOAB is used; in this case, entity sharing between instances is handled explicitly using message passing.  This is described in more detail in Section 5 of this document.

 \ref contents

 \subsection twotwo 2.2. Mesh Entities
MOAB represents the following topological mesh entities: vertex, edge, triangle, quadrilateral, polygon, tetrahedron, pyramid, prism, knife, hexahedron, polyhedron.  MOAB uses the EntityType enumeration to refer to these entity types (see Table 1).  This enumeration has several special characteristics, chosen intentionally: the types begin with vertex, entity types are grouped by topological dimension, with lower-dimensional entities appearing before higher dimensions; the enumeration includes an entity type for sets (described in the next section); and MBMAXTYPE is included at the end of this enumeration, and can be used to terminate loops over type.  In addition to these defined values, the an increment operator (++) is defined such that variables of type EntityType can be used as iterators in loops.
MOAB refers to entities using “handles”.  Handles are implemented as long integer data types, with the four highest-order bits used to store the entity type (mesh vertex, edge, tri, etc.) and the remaining bits storing the entity id.  This scheme is convenient for applications because:
- Handles sort lexicographically by type and dimension; this can be useful for grouping and iterating over entities by type.
- The type of an entity is indicated by the handle itself, without needing to call a function.
- Entities allocated in sequence will typically have contiguous handles; this characteristic can be used to efficiently store and operate on large lists of handles.
.

This handle implementation is exposed to applications intentionally, because of optimizations that it enables, and is unlikely to change in future versions.

  \subsection tableone Table 1: Values defined for the EntityType enumerated type.
<table border="1">
<tr>
<td>MBVERTEX = 0</td>
<td>MBPRISM</td>
</tr>
<tr>
<td>MBEDGE</td>
<td>MBKNIFE</td>
</tr>
<tr>
<td>MBTRI</td>
<td>MBHEX</td>
</tr>
<tr>
<td>MBQUAD</td>
<td>MBPOLYHEDRON</td>
</tr>
<tr>
<td>MBPOLYGON</td>
<td>MBENTITYSET</td>
</tr>
<tr>
<td>MBTET</td>
<td>MBMAXTYPE</td>
</tr>
<tr>
<td>MBPYRAMID</td>
<td></td>
</tr>
</table>

MOAB defines a special class for storing lists of entity handles, named Range.  This class stores handles as a series of (start_handle, end_handle) subrange tuples.  If a list of handles has large contiguous ranges, it can be represented in almost constant size using Range.  Since entities are typically created in groups, e.g. during mesh generation or file import, a high degree of contiguity in handle space is typical.  Range provides an interface similar to C++ STL containers like std::vector, containing iterator data types and functions for initializing and iterating over entity handles stored in the range.  Range also provides functions for efficient Boolean operations like subtraction and intersection.  Most API functions in MOAB come in both range-based and vector-based variants.  By definition, a list of entities stored in an Range is always sorted, and can contain a given entity handle only once.  Range cannot store the handle 0 (zero).

Typical usage of an Range object would look like:

\code
using namespace moab;
   int my_function(Range &from_range) {
          int num_in_range = from_range.size();
          Range to_range;
          Range::iterator rit;
    for (rit = from_range.begin(); rit != from_range.end(); ++rit) {
            EntityHandle this_ent = *rit;
            to_range.insert(this_ent);
          }
        }
\endcode

Here, the range is iterated similar to how std::vector is iterated.

  \ref contents

 \subsection adjacencies 2.2.1. Adjacencies & AEntities 

The term adjacencies is used to refer to those entities topologically connected to a given entity, e.g. the faces bounded by a given edge or the vertices bounding a given region.  The same term is used for both higher-dimensional (or bounded) and lower-dimensional (or bounding) adjacent entities.  MOAB provides functions for querying adjacent entities by target dimension, using the same functions for higher- and lower-dimension adjacencies.  By default, MOAB stores the minimum data necessary to recover adjacencies between entities.  When a mesh is initially loaded into MOAB, only entity-vertex (i.e. “downward”) adjacencies are stored, in the form of entity connectivity.  When “upward” adjacencies are requested for the first time, e.g. from vertices to regions, MOAB stores all vertex-entity adjacencies explicitly, for all entities in the mesh.  Non-vertex entity to entity adjacencies are never stored, unless explicitly requested by the application.

In its most fundamental form, a mesh need only be represented by its vertices and the entities of maximal topological dimension.  For example, a hexahedral mesh can be represented as the connectivity of the hex elements and the vertices forming the hexes.  Edges and faces in a 3D mesh need not be explicitly represented.  We refer to such entities as “AEntities”, where ‘A’ refers to “Auxiliary”, “Ancillary”, and a number of other words mostly beginning with ‘A’.  Individual AEntities are created only when requested by applications, either using mesh modification functions or by requesting adjacencies with a special “create if missing” flag passed as “true”.  This reduces the overall memory usage when representing large meshes.  Note entities must be explicitly represented before they can be assigned tag values or added to entity sets (described in following Sections).

\ref contents

 \subsection twothree 2.3. Entity Sets
Entity sets are used to store arbitrary collections of entities and other sets.  Sets are used for a variety of things in mesh-based applications, from the set of entities discretizing a given geometric model entity to the entities partitioned to a specific processor in a parallel finite element application.  MOAB entity sets can also store parent/child relations with other entity sets, with these relations distinct from contains relations.  Parent/child relations are useful for building directed graphs with graph nodes representing collections of mesh entities; this construct can be used, for example, to represent an interface of mesh faces shared by two distinct collections of mesh regions.  MOAB also defines one special set, the “root set” or the interface itself; all entities are part of this set by definition.  Defining a root set allows the use of a single set of MOAB API functions to query entities in the overall mesh as well as its subsets.

MOAB entity sets can be one of two distinct types: list-type entity sets preserve the order in which entities are added to the set, and can store a given entity handle multiple times in the same set; set-type sets are always ordered by handle, regardless of the order of addition to the set, and can store a given entity handle only once.  This characteristic is assigned when the set is created, and cannot be changed during the set’s lifetime.

MOAB provides the option to track or not track entities in a set.  When entities (and sets) are deleted by other operations in MOAB, they will also be removed from containing sets for which tracking has been enabled.  This behavior is assigned when the set is created, and cannot be changed during the set’s lifetime.  The cost of turning tracking on for a given set is sizeof(EntityHandle) for each entity added to the set; MOAB stores containing sets in the same list which stores adjacencies to other entities.

Using an entity set looks like the following:
\code
using namespace moab;
// load a file using MOAB, putting the loaded mesh into a file set
EntityHandle file_set;
ErrorCode rval = moab->create_meshset(MESHSET_SET, file_set);
rval = moab->load_file(“fname.vtk”, &file_set);
Range set_ents;
// get all the 3D entities in the set
rval = moab->get_entities_by_dimension(file_set, 3, set_ents);
\endcode

Entity sets are often used in conjunction with tags (described in the next section), and provide a powerful mechanism to store a variety of meta-data with meshes.

\ref contents

 \subsection twofour 2.4. Tags 

Applications of a mesh database often need to attach data to mesh entities.  The types of attached data are often not known at compile time, and can vary across individual entities and entity types.  MOAB refers to this attached data as a “tag”.  Tags can be thought of loosely as a variable, which can be given a distinct value for individual entities, entity sets, or for the interface itself.  A tag is referenced using a handle, similarly to how entities are referenced in MOAB.  Each MOAB tag has the following characteristics, which can be queried through the MOAB interface:
- Name
- Size (in bytes)
- Storage type
- Data type (integer, double, opaque, entity handle)
- Handle
.

The storage type determines how tag values are stored on entities.  

- Dense: Dense tag values are stored in arrays which match arrays of contiguous entity handles.  Dense tags are more efficient in both storage and memory if large numbers of entities are assigned the same tag.  Storage for a given dense tag is not allocated until a tag value is set on an entity; memory for a given dense tag is allocated for all entities in a given sequence at the same time.
- Sparse: Sparse tags are stored as a list of (entity handle, tag value) tuples, one list per sparse tag, sorted by entity handle.
- Bit: Bit tags are stored similarly to dense tags, but with special handling to allow allocation in bit-size amounts per entity.
.

MOAB also supports variable-length tags, which can have a different length for each entity they are assigned to.  Variable length tags are stored similarly to sparse tags.

The data type of a tag can either be one understood at compile time (integer, double, entity handle), in which case the tag value can be saved and restored properly to/from files and between computers of different architecture (MOAB provides a native HDF5-based save/restore format for this purpose; see Section 4.6).  The opaque data type is used for character strings, or for allocating “raw memory” for use by applications (e.g. for storage application-defined structures or other abstract data types).  These tags are saved and restored as raw memory, with no special handling for endian or precision differences.

An application would use the following code to attach a double-precision tag to vertices in a mesh, e.g. to assign a temperature field to those vertices:

\code
using namespace moab;
// load a file using MOAB and get the vertices
ErrorCode rval = moab->load_file(“fname.vtk”);
Range verts;
rval = moab->get_entities_by_dimension(0, 0, verts);
// create a tag called “TEMPERATURE”
Tag temperature;
double def_val = -1.0d-300, new_val = 273.0;
rval = moab->tag_create(“TEMPERATURE”, sizeof(double), MB_TAG_DENSE, 
                        MB_TYPE_DOUBLE, temperature, &def_val);
// assign a value to vertices
for (Range::iterator vit = verts.begin(); 
     vit != verts.end(); ++vit)
  rval = moab->tag_set_data(temperature, &(*rit), 1, &new_val);

\endcode

The semantic meaning of a tag is determined by applications using it.  However, to promote interoperability between applications, there are a number of tag names reserved by MOAB which are intended to be used by convention.  Mesh readers and writers in MOAB use these tag conventions, and applications can use them as well to access the same data. Ref. [4] maintains an up-to-date list of conventions for meta-data usage in MOAB.

  \ref contents

  \section api 3.MOAB API Design Philosophy and Summary

This section describes the design philosophy behind MOAB, and summarizes the functions, data types and enumerated variables in the MOAB API.  A complete description of the MOAB API is available in online documentation in the MOAB distribution [8].

MOAB is designed to operate efficiently on collections of entities.  Entities are often created or referenced in groups (e.g. the mesh faces discretizing a given geometric face, the 3D elements read from a file), with those groups having some form of temporal or spatial locality.  The interface provides special mechanisms for reading data directly into the native storage used in MOAB, and for writing large collections of entities directly from that storage, to avoid data copies.  MOAB applications structured to take advantage of that locality will typically operate more efficiently.

MOAB has been designed to maximize the flexibility of mesh data which can be represented.  There is no explicit constraint on the geometric structure of meshes represented in MOAB, or on the connectivity between elements.  In particular, MOAB allows the representation of multiple entities with the same exact connectivity; however, in these cases, explicit adjacencies must be used to distinguish adjacencies with AEntities bounding such entities.

The number of vertices used to represent a given topological entity can vary, depending on analysis needs; this is often the case in FEA.  For example, applications often use “quadratic” or 10-vertex tetrahedral, with vertices at edge midpoints as well as corners.  MOAB does not distinguish these variants by entity type, referring to all variants as “tetrahedra”.  The number of vertices for a given entity is used to distinguish the variants, with canonical numbering conventions used to determine placement of the vertices [9].  This is similar to how such variations are represented in the Exodus [10] and Patran [11] file formats.  In practice, we find that this simplifies coding in applications, since in many cases the handling of entities depends only on the number of corner vertices in the element.  Some MOAB API functions provide a flag which determines whether corner or all vertices are requested.

The MOAB API is designed to balance complexity and ease of use.  This balance is evident in the following general design characteristics:

- Entity lists: Lists of entities are passed to and from MOAB in a variety of forms.  Lists output from MOAB are passed as either STL vector or Range data types.  Either of these constructs may be more efficient in both time and memory, depending on the semantics of the data being requested.  Input lists are passed as either Range’s, or as a pointer to EntityHandle and a size.  The latter allows the same function to be used when passing individual entities, without requiring construction of an otherwise unneeded STL vector.
- Entity sets: Most query functions accept an entity set as input.  Applications can pass zero to indicate a request for the whole interface.  Note that this convention applies only to query functions; attempts to add or subtract entities to/from the interface using set-based modification functions, or to add parents or children to the interface set, will fail.  Allowing specification of the interface set in this manner avoids the need for a separate set of API functions to query the database as a whole.
- Implicit Booleans in output lists: A number of query functions in MOAB allow specification of a Boolean operation (Interface::INTERSECT or Interface::UNION).  This operation is applied to the results of the query, often eliminating the need for code the application would need to otherwise implement.  For example, to find the set of vertices shared by a collection of quadrilaterals, the application would pass that list of quadrilaterals to a request for vertex adjacencies, with Interface::INTERSECT passed for the Boolean flag.  The list of vertices returned would be the same as if the application called that function for each individual entity, and computed the intersection of the results over all the quadrilaterals.  Applications may also input non-empty lists to store the results, in which case the intersection is also performed with entities already in the list.  In many cases, this allows optimizations in both time and memory inside the MOAB implementation. 
.

Since these objectives are at odds with each other, tradeoffs had to be made between them.  Some specific issues that came up are:

- Using ranges: Where possible, entities can be referenced using either ranges (which allow efficient storage of long lists) or STL vectors (which allow list order to be preserved), in both input and output arguments.
- Entities in sets: Accessing the entities in a set is done using the same functions which access entities in the entire mesh.  The whole mesh is referenced by specifying a set handle of zero<sup>3</sup>.
- Entity vectors on input: Functions which could normally take a single entity as input are specified to take a vector of handles instead.  Single entities are specified by taking the address of that entity handle and specifying a list length of one.  This minimizes the number of functions, while preserving the ability to input single entities.<sup>4</sup>
.

Table 2 lists basic data types and enumerated variables defined and used by MOAB.  Values of the ErrorCode enumeration are returned from most MOAB functions, and can be compared to those listed in Appendix [ref-appendix].

MOAB uses several pre-defined tag names to define data commonly found in various mesh-based analyses.  Ref. [4] describes these meta-data conventions in more detail.  These conventions will be added to as new conventions emerge for using sets and tags in MOAB applications.

  \subsection tabletwo Table 2: Basic data types and enums defined in MOAB.

<table border="1">
<tr>
<th>Enum / Type</th>
<th>Description</th>
</tr>
<tr>
<td>ErrorCode</td>
<td>Specific error codes returned from MOAB</td>
</tr>
<tr>
<td>EntityHandle</td>
<td>Type used to represent entity handles</td>
</tr>
<tr>
<td>Tag</td>
<td>Type used to represent tag handles</td>
</tr>
<tr>
<td>TagType</td>
<td>Type used to represent tag storage type</td>
</tr>
<tr>
<td>DataType</td>
<td>Type used to represent tag data type</td>
</tr>
</table>

Table 3 lists the various groups of functions that comprise the MOAB API.  This is listed here strictly as a reference to the various types of functionality supported by MOAB; for a more detailed description of the scope and syntax of the MOAB API, see the online documentation [7].

  \subsection tablethree Table 3: Groups of functions in MOAB API.  See Ref. [7] for more details.

<table border="1">
<tr>
<th>Function group</th>
<th>Examples</th>
<th>Description</th>
</tr>
<tr>
<td>Constructor, destructor, interface</td>
<td>Interface, ~Core, query_interface</td>
<td>Construct/destroy interface; get pointer to read/write interface</td>
</tr>
<tr>
<td>Entity query</td>
<td>get_entities_by_dimension, get_entities_by_handle</td>
<td>Get entities by dimension, type, etc.</td>
</tr>
<tr>
<td>Adjacencies</td>
<td>get_adjacencies, set_adjacencies, add_adjacencies</td>
<td>Get topologically adjacent entities; set or add explicit adjacencies</td>
</tr>
<tr>
<td>Vertex coordinates</td>
<td>get_coords, set_coords</td>
<td>Get/set vertex coordinates</td>
</tr>
<tr>
<td>Connectivity</td>
<td>get_connectivity, set_connectivity</td>
<td>Get/set connectivity of non-vertex entities</td>
</tr>
<tr>
<td>Sets</td>
<td>create_meshset, add_entities, add_parent_child</td>
<td>Create and work with entity sets</td>
</tr>
<tr>
<td>Tags</td>
<td>tag_get_data, tag_create</td>
<td>Create, read, write tag data</td>
</tr>
<tr>
<td>Handles</td>
<td>type_from_handle, id_from_handle</td>
<td>Go between handles and types/ids</td>
</tr>
<tr>
<td>File handling</td>
<td>load_mesh, save_mesh</td>
<td>Read/write mesh files</td>
</tr>
<tr>
<td>Geometric dimension</td>
<td>get_dimension, set_dimension</td>
<td>Get/set geometric dimension of mesh</td>
</tr>
<tr>
<td>Mesh modification</td>
<td>create_vertex, delete_entity</td>
<td>Create or delete mesh entities</td>
</tr>
<tr>
<td>Information</td>
<td>list_entities, get_last_error</td>
<td>Get or print certain information</td>
</tr>
<tr>
<td>High-order nodes</td>
<td>high_order_node</td>
<td>Get information on high-order nodes</td>
</tr>
<tr>
<td>Canonical numbering</td>
<td>side_number</td>
<td>Get canonical numbering information</td>
</tr>
</table>

<sup>3</sup>In iMesh, the whole mesh is specified by a special entity set handle, referred to as the “root set”.

<sup>4</sup>Note that STL vectors of entity handles can be input in this manner by using &vector[0] and vector.size() for the 1d vector address and size, respectively.

 \ref contents

 \section services 4.Related Mesh Services

A number of mesh-based services are often used in conjunction with a mesh library.  For example, parallel applications often need to visualize the mesh and associated data.  Other services, like spatial interpolation or finding the faces on the “skin” of a 3D mesh, can be implemented more efficiently using knowledge of specific data structures in MOAB.  Several of these services provided with MOAB are described in this chapter.

 \ref contents

  \subsection fourone 4.1. Visualization

Visualization is one of the most common needs associated with meshes.  The primary tool used to visualize MOAB meshes is VisIt [12].  Users can specify that VisIt read mesh directly out of the MOAB instance, by specifying the ITAPS_MOAB mesh format and a file readable by MOAB (see http://sigma.mcs.anl.gov/?p=429).

There are some initial capabilities in VisIt for limited viewing and manipulation of tag data and some types of entity sets.  Tag data is visualized using the same mechanisms used to view other field data in VisIt, e.g. using a pseudocolor plot; sets are viewed using VisIt’s SIL window, accessed by selecting the SIL icon in the data selection window.  xxx shows a vertex-based radiation temperature field computed by the Cooper rad-hydro code [1] for a subset of geometric volumes in a mesh.   

Reorganization of VisIt’s set handling is also underway, to increase versatility and flexibility of this important mechanism.

 \ref contents

  \subsection fourtwo 4.2. Parallel Decomposition

To support parallel simulation, applications often need to partition a mesh into parts, designed to balance the load and minimize communication between sets.  MOAB includes the MBZoltan tool for this purpose, constructed on the well-known Zoltan partitioning library [13].  After computing the partition using Zoltan, MBZoltan stores the partition as either tags on individual entities in the partition, or as tagged sets, one set per part.  Since a partition often exhibits locality similar to how the entities were created, storing it as sets (based on Range’s) is often more memory-efficient than an entity tag-based representation.  Figure below shows a couple of partitioned meshes computed with MBZoltan (and visualized in VisIt).


 \image html vis_part.png



 \ref contents

  \subsection fourthree 4.3. Skinner

An operation commonly applied to mesh is to compute the outermost “skin” bounding a contiguous block of elements.  This skin consists of elements of one fewer topological dimension, arranged in one or more topological balls on the boundary of the elements.  The Skinner tool computes the skin of a mesh in a memory-efficient manner.  Skinner uses knowledge about whether vertex-entity adjacencies and AEntities exist to minimize memory requirements and searching time required during the skinning process.  This skin can be provided as a single collection of entities, or as sets of entities distinguished by forward and reverse orientation with respect to higher-dimensional entities in the set being skinned.

The following code fragment shows how Skinner can be used to compute the skin of a range of hex elements:

 \code
using namespace moab;
Range hexes, faces;
ErrorCode rval = moab->get_entities_by_dimension(0, 3, hexes);
Skinner myskinner(moab);
bool verts_too = false;
ErrorCode rval = myskinner.find_skin(hexes, verts_too, faces);
\endcode

Skinner can also skin a mesh based on geometric topology groupings imported with the mesh.  The geometric topology groupings contain information about the mesh “owned” by each of the entities in the geometric model, e.g. the model vertices, edges, etc.  Links between the mesh sets corresponding to those entities can be inferred directly from the mesh.  Skinning a mesh this way will typically be much faster than doing so on the actual mesh elements, because there is no need to create and destroy interior faces on the mesh.

 \ref contents

  \subsection fourfour 4.4. Tree Decompositions

MOAB provides several mechanisms for spatial decomposition and searching in a mesh:

- AdaptiveKDTree: Adaptive KD tree, a space-filling decomposition with axis-aligned splitting planes, enabling fast searching.
- BSPTree: Binary Space Partition tree, with non-axis-aligned partitions, for fast spatial searches with slightly better memory efficiency than KD trees.
- OrientedBoxTreeTool: Oriented Bounding Box tree hierarchy, useful for fast ray-tracing on collections of mesh facets.
.

These trees have various space and time searching efficiencies.  All are implemented based on entity sets and parent/child relations between those sets, allowing storage of a tree decomposition using MOAB’s native file storage mechanism (see Section 4.6.1).  MOAB’s entity set implementation is specialized for memory efficiency when representing binary trees.  Tree decompositions in MOAB have been used to implement fast ray tracing to support radiation transport [14], solution coupling between meshes [2], and embedded boundary mesh generation [15].  MOAB also includes the DAGMC tool, supporting Monte Carlo radiation transport.

The following code fragment shows very basic use of AdaptiveKDTree.  A range of entities is put in the tree; the leaf containing a given point is found, and the entities in that leaf are returned.
\code
using namespace moab;
// create the adaptive kd tree from a range of tets
EntityHandle tree_root
AdaptiveKDTree myTree(moab);
ErrorCode rval = myTree.build_tree(tets, tree_root);

// get the overall bounding box corners
double boxmax[3], boxmin;
rval = myTree.get_tree_box(tree_root, boxmax, boxmin);

// get the tree leaf containing point xyz, and the tets in that leaf
AdaptiveKDTreeIter treeiter;
rval = myTree.leaf_containing_point(tree_root, xyz, treeiter);
Range leaf_tets;
rval = moab->get_entities_by_dimension(treeiter.handle(), 3, 
                                       leaf_tets, false);
\endcode

More detailed examples of using the various tree decompositions in MOAB can be found in [ref-treeexamples].

 \ref contents

  \subsection fourfive 4.5. File Reader/Writer Interfaces

Mesh readers and writers communicate mesh into/out of MOAB from/to disk files.  Reading a mesh often involves importing large sets of data, for example coordinates of all the nodes in the mesh.  Normally, this process would involve reading data from the file into a temporary data buffer, then copying data from there into its destination in MOAB.  To avoid the expense of copying data, MOAB has implemented a reader/writer interface that provides direct access to blocks of memory used to represent mesh.

The reader interface, declared in ReadUtilIface, is used to request blocks of memory for storing coordinate positions and element connectivity.  The pointers returned from these functions point to the actual memory used to represent those data in MOAB.  Once data is written to that memory, no further copying is done.  This not only saves time, but it also eliminates the need to allocate a large memory buffer for intermediate storage of these data. 

MOAB allocates memory for nodes and elements (and their corresponding dense tags) in chunks, to avoid frequent allocation/de-allocation of small chunks of memory.  The chunk size used depends on from where the mesh is being created, and can strongly affect the performance (and memory layout) of MOAB.  Since dense tags are allocated at the chunk size, this can also affect overall memory usage in cases where the mesh size is small but the number of dense tags or dense tag size is large.  When creating vertices and elements through the normal MOAB API, default chunk sizes defined in the SequenceManager class are used.  However, most of the file readers in MOAB allocate the exact amount of space necessary to represent the mesh being read.  There are also a few exceptions to this:

- When compiled in parallel, this space is increased by a factor of 1.5, to allow subsequent creation of ghost vertices/elements in the same chunk as the original mesh.
- The .cub file reader, which creates nodes and elements for individual geometric entities in separate calls, allocates using the default vertex/element sequence sizes, which are defined in the SequenceManager class in MOAB.
.

Applications calling the reader interface functions directly can specify the allocation chunk size as an optional parameter.

The reader interface consists of the following functions:

- get_node_coords: Given the number of vertices requested, the number of geometric dimensions, and a requested start id, allocates a block of vertex handles and returns pointers to coordinate arrays in memory, along with the actual start handle for that block of vertices.
- get_element_connect: Given the number of elements requested, the number of vertices per element, the element type and the requested start id, allocates the block of elements, and returns a pointer to the connectivity array for those elements and the actual start handle for that block.  The number of vertices per element is necessary because those elements may include higher-order nodes, and MOAB stores these as part of the normal connectivity array.
- update_adjacencies: This function takes the start handle for a block of elements and the connectivity of those elements, and updates adjacencies for those elements.  Which adjacencies are updated depends on the options set in AEntityFactory.
.

The following code fragment illustrates the use of ReadUtilIface to read a mesh directly into MOAB’s native representation.  This code assumes that connectivity is specified in terms of vertex indices, with vertex indices starting from 1.

\code
// get the read iface from moab
ReadUtilIface *iface;
ErrorCode rval = moab->get_interface("ReadUtilIface", &iface);

// allocate a block of vertex handles and read xyz’s into them
std::vector<double*> arrays;
EntityHandle startv, *starth;
rval = iface->get_node_coords(3, num_nodes, 0, startv, arrays);
for (int i = 0; i < num_nodes; i++)
  infile >> arrays[0][i] >> arrays[1][i] >> arrays[2][i];

// allocate block of hex handles and read connectivity into them
rval = iface->get_element_connect(num_hexes, 8, MBHEX, 0, starth);
for (int i = 0; i < 8*num_hexes; i++)
  infile >> starth[i];

// change connectivity indices to vertex handles
for (int i = 0; i < 8*num_hexes; i++)
  starth[i] += startv-1;
\endcode

The writer interface, declared in WriteUtilIface, provides functions that support writing vertex coordinates and element connectivity to storage locations input by the application.  Assembling these data is a common task for writing mesh, and can be non-trivial when exporting only subsets of a mesh.  The writer interface declares the following functions:

- get_node_coords: Given already-allocated memory and the number of vertices and dimensions, and a range of vertices, this function writes vertex coordinates to that memory.  If a tag is input, that tag is also written with integer vertex ids, starting with 1, corresponding to the order the vertices appear in that sequence (these ids are used to write the connectivity array in the form of vertex indices).
- get_element_connect: Given a range of elements and the tag holding vertex ids, and a pointer to memory, the connectivity of the specified elements are written to that memory, in terms of the indices referenced by the specified tag.  Again, the number of vertices per element is input, to allow the direct output of higher-order vertices.
- gather_nodes_from_elements: Given a range of elements, this function returns the range of vertices used by those elements.  If a bit-type tag is input, vertices returned are also marked with 0x1 using that tag.  If no tag is input, the implementation of this function uses its own bit tag for marking, to avoid using an n2 algorithm for gathering vertices.
- reorder: Given a permutation vector, this function reorders the connectivity for entities with specified type and number of vertices per entity to match that permutation.  This function is needed for writing connectivity into numbering systems other than that used internally in MOAB.
.

The following code fragment shows how to use WriteUtilIface to write the vertex coordinates and connectivity indices for a subset of entities.

\code
using namespace moab;
// get the write iface from moab
WriteUtilIface *iface;
ErrorCode rval = moab->get_interface("WriteUtilIface", &iface);

// get all hexes the model, and choose the first 10 of those
Range tmp_hexes, hexes, verts;
rval = moab->get_entities_by_type(0, MBHEX, tmp_hexes);
for (int i = 0; i < 10; i++) hexes.insert(tmp_hexes[i]);
rval = iface->gather_nodes_from_elements(hexes, 0, verts);

// assign vertex ids
iface->assign_ids(verts, 0, 1);

// allocate space for coordinates & write them
std::vector<double*> arrays(3);
for (int i = 0; i < 3; i++) arrays[i] = new double[verts.size()];
iface->get_node_coords(3, verts.size(), verts, 0, 1, arrays);

// put connect’y in array, in the form of indices into vertex array
std::vector<int> conn(8*hexes.size());
iface->get_element_connect(hexes.size(), 8, 0, hexes, 0, 1, &conn[0]);
\endcode

 \ref contents

  \subsection foursix 4.6. File Readers/Writers Packaged With MOAB

MOAB has been designed to efficiently represent data and metadata commonly found in finite element mesh files.  Readers and writers are included with MOAB which import/export specific types of metadata in terms of MOAB sets and tags, as described earlier in this document.  The number of readers and writers in MOAB will probably grow over time, and so they are not enumerated here.  See the src/io/README file in the MOAB source distribution for a current list of supported formats.

Because of its generic support for readers and writers, described in the previous section, MOAB is also a good environment for constructing new mesh readers and writers.  The ReadTemplate and WriteTemplate classes in src/io are useful starting points for constructing new file readers/writers; applications are encouraged to submit their own readers/writers for inclusion in MOAB’s contrib/io directory in the MOAB source. 

The usefulness of a file reader/writer is determined not only by its ability to read and write nodes and elements, but also in its ability to store the various types of meta-data included with the typical mesh.  MOAB readers and writers are distinguished by their ability to preserve meta-data in meshes that they read and write.  For example, MOAB’s CUB reader imports not only the mesh saved from CUBIT, but also the grouping of mesh entities into sets which reflect the geometric topology of the model used to generate the mesh.  See [4] for a more detailed description of meta-data conventions used in MOAB’s file readers and writers, and the individual file reader/writer header files in src/io for details about the specific readers and writers.

Three specific file readers in MOAB bear further discussion: MOAB’s native HDF5-based file reader/writer; the CUB reader, used to import mesh and meta-data represented in CUBIT; and the CGM reader, which imports geometric models.  These are described next.

 \ref contents

  \subsection native 4.6.1. Native HD5-Based Reader/Writer

A mesh database must be able to save and restore the data stored in its data model, at least to the extent to which it understands the semantics of that data.  MOAB defines an HDF5-based file format that can store data embedded in MOAB.  By convention, these files are given an “.h5m” file extension.  When reading or writing large amounts of data, it is recommended to use this file format, as it is the most complete and also the most efficient of the file readers/writers in MOAB. 

  \subsection cub 4.6.2. CUB Reader

CUBIT is a toolkit for generating tetrahedral and hexahedral finite element meshes from solid model geometry [16].  This tool saves and restores data in a custom “.cub” file, which stores both mesh and geometry (and data relating the two).  The CUB reader in MOAB can import and interpret much of the meta-data information saved in .cub files.  Ref. [4] describes the conventions used to store this meta-data in the MOAB data model.  The information read from .cub files, and stored in the MOAB data model, includes:

- Geometric model entities and topology
- Model entity names and ids
- Groups, element blocks, nodesets, and sidesets, including model entities stored in them
- Mesh scheme and interval size information assigned to model entities
.

Note that although information about model entities is recovered, MOAB by default does not depend on a solid modeling engine; this information is stored in the form of entity sets and parent/child relations between them.  See Ref. [4] for more information.

 \ref contents

  \subsection cgm 4.6.3. CGM Reader

The Common Geometry Module (CGM) [17] is a library for representing solid model and other types of solid geometry data.  The CUBIT mesh generation toolkit uses CGM for its geometric modeling support, and CGM can restore geometric models in the exact state in which they were represented in CUBIT.  MOAB contains a CGM reader, which can be enabled with a configure option.  Using this reader, MOAB can read geometric models, and represent their model topology using entity sets linked by parent/child relations.  The mesh in these models comes directly from the modeling engine faceting routines; these are the same facets used to visualize solid models in other graphics engines.  When used in conjunction with the VisIt visualization tool (see Section 4.1), this provides a solution for visualizing geometric models.  The figure below  shows a model imported using MOAB’s CGM reader and visualized with VisIt.

\image html simple.png

\ref contents

 \subsection fourseven 4.7. AHF Representation

Currently, the upward (vertex to entities) adjacencies are created and stored the first time a query requiring the adjacency is performed. Any non-vertex entity to entity adjacencies are performed using boolean operations on vertex-entity adjacencies. Because of this approach, such adjacency queries might become expensive with increasing dimension. We have added an alternative approach for obtaining adjacencies using the Array-based Half-Facet (AHF) representation[23]. The AHF uses sibling half-facets as a core abstraction which are generalizations of the opposite half-edge and half-face data structure for 2D and 3D manifold meshes. The AHF data structure consists of two essential maps: 1) the mapping between all sibling half-facets (sibhfs) and, 2) the mapping from each vertex to some incident half-facet (v2hf). The entire range of adjacencies (higher-, same- and lower-dimension) are computed using these two maps.

The easiest way to avail this feature is to configure MOAB with " --enable-ahf " option. The adjacency queries can then be performed through calls to the preserved interface function "get_adjacencies" returning the values in a standard vector. Currently, returning adjacent entityhandles in MOAB::Range is not supported for AHF-based queries. There is one key difference between the native MOAB (adjacency-list based) and AHF based adjacency calls using the "get_adjacencies" interface. In native MOAB adjacency calls, the same-dimensional queries return the query entities whereas for AHF it would return the same-dimensional entities connected to the query entities via a lower dimensional facet. Thus the entire range ( higher-dimensional, same-dimensional, lower-dimensional) of adjacencies can be obtained using the same interface. Similar to MOAB's native adjacency lists, the AHF maps are created during the first adjacency call which will make the first adjacency call expensive.

In the current release, AHF based adjacencies calls do not support the following cases:
  - polygon/polyhedral meshes,
  - mixed entity type meshes,
  - meshsets,
  - create_if_missing option set to true, and
  - modified meshes.

The support for these would be gradually added in the next releases. In these cases, any adjacency call would revert back to MOAB's native adjacency list based queries.

If for some reason, the user does not want to configure MOAB with AHF but would still like to use the AHF-based adjacencies for certain queries, they could use the following three interface functions provided in the HalfFacetRep class which implements the AHF maps and adjacency queries:
  - initialize : This function creates all the necessary AHF maps for the input mesh and hence should be called before any adjacency calls are made.
  - get_adjacencies: Function for adjacency calls.
  - deinitialize: This function deletes all the AHF maps and should be called after all AHF-based adjacency calls have been performed.

TODO:: Other features to be added
  - obtain ring neighborhoods with support for half-rings
  - efficient extraction of boundaries

  \ref contents

  \subsection foureight 4.8. Uniform Mesh Refinement
  Many applications require a hierarchy of successively refined meshes for a number of purposes such as convergence studies, to use multilevel methods like multigrid, generate large meshes in parallel computing to support increasing mesh sizes with increase in number of processors, etc. Uniform mesh refinement provides a simple and efficient way to generate such hierarchies via successive refinement of the mesh at a previous level. It also provides a natural hierarchy via parent and child type of relationship between entities of meshes at different levels. Generally, the standard nested refinement patterns used are the subdivision schemes from 1 to 4 for 2D (triangles, quads) and 1 to 8 for 3D (tets, hexes) entity types. However, many applications might require degree 3 or more for p-refinements.

 MOAB supports generation of a mesh hierarchy i.e., a sequence of meshes with user specified degrees for each level of refinement, from an initial unstructured mesh with support for higher degrees of refinement (supported degrees are listed later).  Thus MOAB supports multi-degree and multi-level mesh generation via uniform refinement. The following figure shows the initial and most refined mesh for four simple meshes to illustrate the multi-degree capability.

  \image html uref_allEtype.png "Uniform Refinement of 2D and 3D meshes"

  Applications using mesh hierarchies require two types of mesh access: intralevel and interlevel. The intralevel access involves working with the mesh at a particular level whereas interlevel access involves querying across different levels. In order to achieve data locality with reduced cache misses for efficient intralevel mesh access, old vertices in the previous i.e. immediate parent mesh are duplicated in the current level. All the entities thus created for the current level use the new entityhandles of old vertices along with the handles of the new vertices. This design makes mesh at each level of the hierarchy independent of those at previous levels. For each mesh in the hierarchy, a MESHSET is created and all entities of the mesh are added to this MESHSET. Thus the meshes of the hierarchy are accessible via these level-wise MESHSET handles.

 For interlevel queries, separate interface functions are defined to allow queries across different levels. These queries mainly allow obtaining the parent at some specified parent level for a child from later levels and vice versa. The child-parent queries are not restricted to a level difference of one as the internal array-based layout of the memory allows traversing between different levels via index relations.

 The hierarchy generation capability is implemented in NestedRefine class. In Table 4, the user interface functions are briefly described. The main hierarchy generating function takes as input a sequence of refinement degrees to be used for generating mesh at each level from the previous level, the total number of levels in the hierarchy. It returns EntityHandles for the meshsets created for the mesh at each level. The number of levels in the hierarchy is prefixed (by the user) and cannot change during hierarchy generation. An example of how to generate a hierarchy can be found  under examples/UniformRefinement.cpp.

 The next three functions for getting the coordinates, connectivity and adjacencies are standard operations to access and query a mesh. The coordinates and connectivity functions, which are similar to their counterparts in MOAB Interface class, allows one to query the mesh at a specific level via its level index. The reason to provide such similar functions is to increase computational efficiency by utilizing the direct access to memory pointers to EntitySequences created during hierarchy generation under the NestedRefine class instead of calling the standard interfaces where a series of calls have to made before the EntitySequence of the requested entity is located. It is important to note that any calls to the standard interfaces available under Interface class should work.

 The underlying mesh data structure used for uniform refinement is the AHF datastructure implemented in HalfFacetRep class. This direct dependence on HalfFacetRep class removes the necessity to configure MOAB with --enable-ahf flag in order to use the uniform refinement capability. During the creation of an object of the NestedRefine class, it will internally initialize all the relevant AHF maps for the mesh in memory. During mesh hierarchy generation, after the creation of mesh at each level all the relevant AHF maps are updated to allow query over it.

 For interlevel queries, currently three kinds of queries are supported. The child to parent function allows querying for a parent in any of the previous levels (including the initial mesh). The parent to child, on the other hand, returns all its children from a requested child level. These two types of queries are only supported for entities, not vertices. A separate vertex to entities function is provided which returns entities from the previous level that are either incident or contain this vertex.

  \subsection tablefour Table 4: User interface functions NestedRefine class.

<table border="1">
<tr>
<th>Function group</th>
<th>Function</th>
<th>Description</th>
</tr>
<tr>
<td>Hierarchy generation</td>
<td>generate_mesh_hierarchy</td>
<td>Generate a mesh hierarchy with a given sequence of refinement degree for each level. </td>
</tr>
<tr>
<td>Vertex coordinates</td>
<td>get_coordinates</td>
<td>Get vertex coordinates</td>
</tr>
<tr>
<td>Connectivity</td>
<td>get_connectivity</td>
<td>Get connectivity of non-vertex entities</td>
</tr>
<tr>
<td>Adjacencies</td>
<td>get_adjacencies</td>
<td>Get topologically adjacent entities</td>
</tr>
<tr>
<td>Interlevel Queries</td>
<td>child_to_parent</td>
<td>Get the parent entity of a child from a specified parent level</td>
</tr>
<tr>
<td>Interlevel Queries</td>
<td>parent_to_child</td>
<td>Get all children of the parent from a specified child level</td>
</tr>
<tr>
<td>Interlevel Queries</td>
<td>vertex_to_entities</td>
<td>Get all the entities in the previous level incident on or containing the vertex</td>
</tr>
</table>

  \subsection tablefive Table 5: The refinement degrees currently supported.

<table border="1">
<tr>
<th>Mesh Dimension</th>
<th>EntityTypes</th>
<th>Degree</th>
<th>Number of children</th>
</tr>
<tr>
<td>1</td>
<td>MBEDGE</td>
<td>2, 3, 5</td>
<td>2, 3, 5 </td>
</tr>
<tr>
<td>2</td>
<td>MBTRI, MBQUAD</td>
<td>2, 3, 5</td>
<td>4, 9, 25 </td>
</tr>
<tr>
<td>3</td>
<td>MBTET, MBHEX</td>
<td>2, 3</td>
<td>8, 27 </td>
</tr>
</table>


In Table 5, the currently supported degrees of refinement for each dimension is listed along with the number of children created for each such degree for a single entity. The following figure shows the cpu times(serial run) for generating hierarchies with various degrees of refinement for each dimension and can be used by the user to guide in choosing the degrees of refinement for the hierarchy. For example, if a multilevel hierarchy is required, a degree 2 refinement per level would give a gradually increasing mesh with more number of levels. If a very refined mesh is desired quickly, then a small hierarchy with high-order refinement should be generated.

\image html uref_timeEtype.png "Mesh sizes Vs. Time"


  Current support:
  - Single dimension(i.e, no curves in surface meshes or curves/surfaces in volume meshes)
  - Linear point projection
  - Serial

  TODO:
   - Mixed-dimensional meshes
   - Mixed-entity types
   - High-Order point projection
   - Parallel


 \ref contents

  \section parallel 5.Parallel Mesh Representation and Query

A parallel mesh representation must strike a careful balance between providing an interface to mesh similar to that of a serial mesh, while also allowing the discovery of parallel aspects of the mesh and performance of parallel mesh-based operations efficiently.  MOAB supports a spatial domain-decomposed view of a parallel mesh, where each subdomain is assigned to a processor, lower-dimensional entities on interfaces between subdomains are shared between processors, and ghost entities can be exchanged with neighboring processors.  Locally, each subdomain, along with any locally-represented ghost entities, are accessed through a local MOAB instance.  Parallel aspects of the mesh, e.g. whether entities are shared, on an interface, or ghost entities, are embedded in the same data model (entities, sets, tags, interface) used in the rest of MOAB.  MOAB provides a suite of parallel functions for initializing and communicating with a parallel mesh, along with functions to query the parallel aspects of the mesh.

  \ref contents

  \subsection fiveone 5.1. Nomenclature & Representation

Before discussing how to access parallel aspects of a mesh, several terms need to be defined:  

<B>Shared entity:</B> An entity shared by one or several other processors.

<B>Multi-shared entity:</B> An entity shared by more than two processors.

<B>Owning Processor:</B> Each shared entity is owned by exactly one processor.  This processor has the right to set tag values on the entity and have those values propagated to any sharing processors.  

<B>Part:</B> The collection of entities assigned to a given processor.  When reading mesh in parallel, the entities in a Part, along with any lower-dimensional entities adjacent to those, will be read onto the assigned processor.

<B>Partition:</B> A collection of Parts which take part in parallel collective communication, usually associated with an MPI communicator.

<B>Interface:</B> A collection of mesh entities bounding entities in multiple parts.  Interface entities are owned by a single processor, but are represented on all parts/processors they bound.

<B>Ghost entity:</B> A shared, non-interface, non-owned entity.

<B>Parallel status:</B> A characteristic of entities and sets represented in parallel. The parallel status, or “pstatus”, is represented by a bit field stored in an unsigned character, with bit values as described in Table 6.

  \subsection tablesix Table 6: Bits representing various parallel characteristics of a mesh.  Also listed are enumerated values that can be used in bitmask expressions; these enumerated variables are declared in MBParallelConventions.h.

<table border="1">
<tr>
<th>Bit</th>
<th>Name</th>
<th>Represents</th>
</tr>
<tr>
<td>0x1</td>
<td>PSTATUS_NOT_OWNED</td>
<td>Not owned by the local processor</td>
</tr>
<tr>
<td>0x2</td>
<td>PSTATUS_SHARED</td>
<td>Shared by exactly two processorstd>
</tr>
<tr>
<td>0x4</td>
<td>PSTATUS_MULTISHARED</td>
<td>Shared by three or more processors</td>
</tr>
<tr>
<td>0x8</td>
<td>PSTATUS_INTERFACE</td>
<td>Part of lower-dimensional interface shared by multiple processors</td>
</tr>
<tr>
<td>0x10</td>
<td>PSTATUS_GHOST</td>
<td>Non-owned, non-interface entities represented locally</td>
</tr>
</table>

Parallel functionality is described in the following sections.  First, methods to load a mesh into a parallel representation are described; next, functions for accessing parallel aspects of a mesh are described; functions for communicating mesh and tag data are described.

  \ref contents

  \subsection fivetwo 5.2. Parallel Mesh Initialization

Parallel mesh is initialized in MOAB in several steps:

-#  Establish a local mesh on each processor, either by reading the mesh into that representation from disk, or by creating mesh locally through the normal MOAB interface.  
-#  Find vertices, then other entities, shared between processors, based on a globally-consistent vertex numbering stored on the GLOBAL_ID tag.  
-#  Exchange ghost or halo elements within a certain depth of processor interfaces with neighboring processors.  
.

These steps can be executed by a single call to MOAB’s load_file function, using the procedure described in Section 5.2.1.  Or, they can be executed in smaller increments calling specific functions in MOAB’s ParallelComm class, as described in Section 5.2.2.  Closely related to the latter method is the handling of communicators, described in more detail in Section.

  \ref contents

  \subsection initialization 5.2.1. Parallel Mesh Initialization by Loading a File

In the file reading approach, a mesh must contain some definition of the partition (the assignment of mesh, usually regions, to processors).  Partitions can be derived from other set structures already on the mesh, or can be computed explicitly for that purpose by tools like mbzoltan (see Section 4.2).  For example, geometric volumes used to generate the mesh, and region-based material type assignments, are both acceptable partitions (see Ref. [4] for information about this and other meta-data often accompanying mesh).  In addition to defining the groupings of regions into parts, the assignment of specific parts to processors can be done implicitly or using additional data stored with the partition.

MOAB implements several specific methods for loading mesh into a parallel representation:

- READ_PART: each processor reads only the mesh used by its part(s).

- READ_DELETE: each processor reads the entire mesh, then deletes the mesh not used by its part(s).

- BCAST_DELETE: the root processor reads and broadcasts the mesh; each processor then deletes the mesh not used by its part(s).

The READ_DELETE and BCAST_DELETE methods are supported for all file types MOAB is able to read, while READ_PART is only implemented for MOAB’s native HDF5-based file format.

Various other options control the selection of part sets or other details of the parallel read process.  For example, the application can specify the tags, and optionally tag values, which identify parts, and whether those parts should be distributed according to tag value or using round-robin assignment.

The options used to specify loading method, the data used to identify parts, and other parameters controlling the parallel read process, are shown in Table 7.
  \subsection tableseven Table 7: Options passed to MOAB’s load_file function identifying the partition and other parameters controlling the parallel read of mesh data.  Options and values should appear in option string as “option=val”, with a delimiter (usually “;”) between options.

<table border="1">
<tr>
<th>Option</th>
<th>Value</th>
<th>Description</th>
</tr>
<tr>
<td>PARTITION</td>
<td><tag_name></td>
<td>Sets with the specified tag name should be used as part sets</td>
</tr>
<tr>
<td>PARTITION_VAL</td>
<td><val1, val2-val3, ...></td>
<td>Integer values to be combined with tag name, with ranges input using val1, val2-val3.  Not meaningful unless PARTITION option is also given.</td>
</tr>
<tr>
<td>PARTITION_DISTRIBUTE</td>
<td>(none)</td>
<td>If present, or values are not input using PARTITION_VAL, sets with tag indicated in PARTITION option are partitioned across processors in round-robin fashion.</td>
</tr>
<tr>
<td>PARALLEL_RESOLVE_SHARED_ENTS</td>
<td><pd.sd></td>
<td>Resolve entities shared between processors, where partition is made up of pd- dimensional entities, and entities of dimension sd and lower should be resolved.</td>
</tr>
<tr>
<td>PARALLEL_GHOSTS</td>
<td><gd.bd.nl[.ad]></td>
<td>Exchange ghost elements at shared inter-processor interfaces.  Ghost elements of dimension gd will be exchanged.  Ghost elements are chosen going through bd-dimensional interface entities.  Number of layers of ghost elements is specified in nl.  If ad is present, lower-dimensional entities bounding exchanged ghost entities will also be exchanged; allowed values for ad are 1 (exchange bounding edges), 2 (faces), or 3 (edges and faces).</td>
</tr>
<td>PARALLEL_COMM</td>
<td><id></td>
<td>Use the ParallelComm with index <id>.  Index for a ParallelComm object can be checked with ParallelComm::get_id(), and a ParallelComm with a given index can be retrieved using ParallelComm::get_pcomm(id).</td>
</tr>
<tr>
<tr>
<td>CPUTIME</td>
<td>(none)</td>
<td>Print cpu time required for each step of parallel read & initialization.</td>
</tr>
<tr>
<td>MPI_IO_RANK</td>
<td><r></td>
<td>If read method requires reading mesh onto a single processor, processor with rank r is used to do that read.</td>
</tr>
</table>

Several example option strings controlling parallel reading and initialization are:

<B>“PARALLEL=READ_DELETE; PARTITION=MATERIAL_SET; PARTITION_VAL=100, 200, 600-700”:</B> The whole mesh is read by every processor; this processor keeps mesh in sets assigned the tag whose name is “MATERIAL_SET” and whose value is any one of 100, 200, and 600-700 inclusive.

<B>“PARALLEL=READ_PART; PARTITION=PARALLEL_PARTITION, PARTITION_VAL=2”:</B> Each processor reads only its mesh; this processor, whose rank is 2, is responsible for elements in a set with the PARALLEL_PARTITION tag whose value is 2.  This would by typical input for a mesh which had already been partitioned with e.g. Zoltan or Parmetis.

<B>“PARALLEL=BCAST_DELETE; PARTITION=GEOM_DIMENSION, PARTITION_VAL=3, PARTITION_DISTRIBUTE”:</B> The root processor reads the file and broadcasts the whole mesh to all processors.  If a list is constructed with entity sets whose GEOM_DIMENSION tag is 3, i.e. sets corresponding to geometric volumes in the original geometric model, this processor is responsible for all elements with index R+iP, i >= 0 (i.e. a round-robin distribution).

 \ref contents

  \subsection functions 5.2.2. Parallel Mesh Initialization Using Functions

After creating the local mesh on each processor, an application can call the following functions in ParallelComm to establish information on shared mesh entities.  See the [http://ftp.mcs.anl.gov/pub/fathom/moab-docs/HelloParMOAB_8cpp-example.html example] in the MOAB source tree for a complete example of how this is done from an application.

- ParallelComm::resolve_shared_entities (collective): Resolves shared entities between processors, based on GLOBAL_ID tag values of vertices.  Various forms are available, based on entities to be evaluated and maximum dimension for which entity sharing should be found.

- ParallelComm::exchange_ghost_cells (collective): Exchange ghost entities with processors sharing an interface with this processor, based on specified ghost dimension (dimension of ghost entities exchanged), bridge dimension, number of layers, and type of adjacencies to ghost entities.  An entity is sent as a ghost if it is within that number of layers, across entities of the bridge dimension, with entities owned by the receiving processor, or if it is a lower-dimensional entity adjacent to a ghost entity and that option is requested.
.

 \ref contents

  \subsection communicator 5.2.3. Communicator Handling

The ParallelComm constructor takes arguments for an MPI communicator and a MOAB instance.  The ParallelComm instance stores the MPI communicator, and registers itself with the MOAB instance.  Applications can specify the ParallelComm index to be used for a given file operation, thereby specifying the MPI communicator for that parallel operation.  For example:

\code
using namespace moab;
// pass a communicator to the constructor, getting back the index
MPI_Comm my_mpicomm;
int pcomm_index;
ParallelComm my_pcomm(moab, my_mpicomm, &pcomm_index);

// write the pcomm index into a string option
char load_opt[32];
sprintf(load_opt, "PARALLEL=BCAST_DELETE;PARALLEL_COMM=%d", 
   pcomm_index);

// specify that option in a parallel read operation
ErrorCode rval = moab->load_file(load_opt, fname, ...)
\endcode

In the above code fragment, the ParallelComm instance with index pcomm_index will be used in the parallel file read, so that the operation executes over the specified MPI communicator.  If no ParallelComm instance is specified for a parallel file operation, a default instance will be defined, using MPI_COMM_WORLD.

Applications needing to retrieve a ParallelComm instance created previously and stored with the MOAB instance, e.g. by a different code component, can do so using a static function on ParallelComm:

\code
ParallelComm *my_pcomm = ParallelComm::get_pcomm(moab, pcomm_index);
\endcode

ParallelComm also provides the ParallelComm::get_all_pcomm function, for retrieving all ParallelComm instances stored with a MOAB instance.  For syntax and usage of this function, see the MOAB online documentation for ParallelComm.hpp [8].

 \ref contents

  \subsection fivethree 5.3. Parallel Mesh Query Functions

Various functions are commonly used in parallel mesh-based applications.  Functions marked as being collective must be called collectively for all processors that are members of the communicator associated with the ParallelComm instance used for the call.

<B>ParallelComm::get_pstatus:</B>  Get the parallel status for the entity.

<B>ParallelComm::get_pstatus_entities:</B> Get all entities whose pstatus satisfies (pstatus & val).

<B>ParallelComm::get_owner:</B> Get the rank of the owning processor for the specified entity.

<B>ParallelComm::get_owner_handle:</B> Get the rank of the owning processor for the specified entity, and the entity's handle on the owning processor.

<B>ParallelComm::get_sharing_data:</B> Get the sharing processor(s) and handle(s) for an entity or entities.  Various overloaded versions are available, some with an optional “operation” argument, where Interface::INTERSECT or Interface::UNION can be specified.  This is similar to the operation arguments to Interface::get_adjacencies.

<B>ParallelComm::get_shared_entities:</B> Get entities shared with the given processor, or with all processors.  This function has optional arguments for specifying dimension, whether interface entities are requested, and whether to return just owned entities.

<B>ParallelComm::get_interface_procs:</B> Return all processors with whom this processor shares an interface.

<B>ParallelComm::get_comm_procs:</B> Return all processors with whom this processor communicates.

  \ref contents

  \subsection fivefour 5.4. Parallel Mesh Communication

Once a parallel mesh has been initialized, applications can call the ParallelComm::exchange_tags function for exchanging tag values between processors.  This function causes the owning processor to send the specified tag values for all shared, owned entities to other processors sharing those entities.  Asynchronous communication is used to hide latency, and only point-to-point communication is used in these calls.

  \ref contents

  \section applications 6.Building MOAB-Based Applications

There are two primary mechanisms supported by MOAB for building applications, one based on MOAB-defined make variables, and the other based on the use of libtool and autoconf.  Both assume the use of a “make”-based build system.  

The easiest way to incorporate MOAB into an application’s build process is to include the “moab.make” file into the application’s Makefile, adding the make variables MOAB_INCLUDES and MOAB_LIBS_LINK to application’s compile and link commands, respectively.  MOAB_INCLUDES contains compiler options specifying the location of MOAB include files, and any preprocessor definitions required by MOAB.  MOAB_LIBS_LINK contains both the options telling where libraries can be found, and the link options which incorporate those libraries into the application.  Any libraries depended on by the particular configuration of MOAB are included in that definition, e.g. the HDF5 library.  Using this method to incorporate MOAB is the most straightforward; for example, the following Makefile is used to build one of the example problems packaged with the MOAB source:
\code
include ${MOAB_LIB_DIR}/moab.make

GetEntities : GetEntities.o
	${CXX} $< ${MOAB_LIBS_LINK} -o $@

.cpp.o : 
	${CXX} ${MOAB_INCLUDES} -c $<
\endcode

Here, the MOAB_LIB_DIR environment variable or make argument definition specifies where the MOAB library is installed; this is also the location of the moab.make file.  Once that file has been included, MOAB_INCLUDES and MOAB_LIBS_LINK can be used, as shown.

Other make variables are defined in the moab.make file which simplify building applications:

- MOAB_LIBDIR, MOAB_INCLUDEDIR: the directories into which MOAB libraries and include files will be installed, respectively.  Note that some include files are put in a subdirectory named “moab” below that directory, to reflect namespace naming conventions used in MOAB.

- MOAB_CXXFLAGS, MOAB_CFLAGS, MOAB_LDFLAGS: Options passed to the C++ and C compilers and the linker, respectively.

- MOAB_CXX, MOAB_CC, MOAB_FC: C++, C, and Fortran compilers specified to MOAB at configure time, respectively.
.

The second method for incorporating MOAB into an application’s build system is to use autoconf and libtool.  MOAB is configured using these tools, and generates the “.la” files that hold information on library dependencies that can be used in application build systems also based on autoconf and libtool.  Further information on this subject is beyond the scope of this User’s Guide; see the “.la” files as installed by MOAB, and contact the MOAB developer’s mailing list [6] for more details.

  \ref contents

  \section implementation  7.iMesh (ITAPS Mesh Interface) Implementation in MOAB

iMesh is a common API to mesh data developed as part of the Interoperable Tools for Advanced Petascale Simulations (ITAPS) project [18].  Applications using the iMesh interface can operate on any implementation of that interface, including MOAB.  MOAB-based applications can take advantage of other services implemented on top of iMesh, including the MESQUITE mesh improvement toolkit [19] and the GRUMMP mesh improvement library [20].

MOAB’s native interface is accessed through the Interface abstract C++ base class.  Wrappers are not provided in other languages; rather, applications wanting to access MOAB from those languages should do so through iMesh.  In most cases, the data models and functionality available through MOAB and iMesh are identical.  However, there are a few differences, subtle and not-so-subtle, between the two:

<B>SPARSE tags used by default:</B> MOAB’s iMesh implementation creates SPARSE tags by default, because of semantic requirements of other tag-related functions in iMesh.  To create DENSE tags through iMesh, use the iMesh_createTagWithOptions extension function (see below).

<B>Higher-order elements:</B> ITAPS currently handles higher-order elements (e.g. a 10-node tetrahedron) usi[21]<sup>5</sup>.  As described in [sec-entities], MOAB supports higher-order entities by allowing various numbers of vertices to define topological entities like quadrilateral or tetrahedron.  Applications can specify flags to the connectivity and adjacency functions specifying whether corner or all vertices are requested.

<B>Self-adjacencies:</B> In MOAB’s native interface, entities are always self-adjacent<sup>6</sup>; that is, adjacencies of equal dimension requested from an entity will always include that entity, while from iMesh will not include that entity.

<B>Option strings:</B> The iMesh specification requires that options in the options string passed to various functions (e.g. iMesh_load) be prepended with the implementation name required to parse them, and delimited with spaces.  Thus, a MOAB-targeted option would appear as “moab:PARALLEL=READ_PART moab:PARTITION=MATERIAL_SET”.

To provide complete MOAB support from other languages through iMesh, a collection of iMesh extension functions are also available.  A general description of these extensions appears below; for a complete description, see the online documentation for iMesh-extensions.h [8].

- Recursive get_entities functions: There are many cases where sets include other sets (see [4] for more information).  MOAB provides iMesh_getEntitiesRec, and other recursive-supporting functions, to get all non-set entities of a given type or topology accessible from input set(s).  Similar functions are available for number of entities of a given type/topology.

- Get entities by tag, and optionally tag value: It is common to search for entities with a given tag, and possibly tag value(s); functions like iMesh_getEntitiesByTag are provided for this purpose.

- Options to createTag: To provide more control over the tag type, the iMesh_createTagWithOptions is provided.  The storage type is controlled with the “

- MBCNType: Canonical numbering evaluations are commonly needed by applications, e.g. to apply boundary conditions locally.  The MBCN package provides these evaluations in terms of entity types defined in MOAB [9]; the getMBCNType is required to translate between iMesh_Topology and MBCN type.

- Iterator step: Step an iterator a specified number of entities; allows advancement of an iterator without needing to allocate memory to hold the entity handles stepped over.

- Direct access to tag storage: The Interface::tag_iterate function allows an application get a pointer to the memory used to store a given tag.  For dense tags on contiguous ranges of entities, this provides more efficient access to tags.  The iMesh functionn iMesh_tagIterate provides access to this functionlity.  See examples/TagIterateC.c and examples/TagIterateF.F for examples of how to use this from C and Fortran, respectively. 
.

As required by the iMesh specification, MOAB generates the “iMesh-Defs.inc” file and installs it with the iMesh and MOAB libraries.  This file defines make variables which can be used to build iMesh-based applications.  The method used here is quite similar to that used for MOAB itself (see Section 6).  In particular, the IMESH_INCLUDES and IMESH_LIBS make variables can be used with application compile and link commands, respectively, with other make variables similar to those provided in moab.make also available.

Note that using the iMesh interface from Fortran-based applications requires a compiler that supports Cray pointers, along with the pass-by-value (%VAL) extension.  Almost all compilers support those extensions; however, if using the gcc series of compilers, you must use gfortran 4.3 or later.

<sup>5</sup>There are currently no implementations of this interface.

<sup>6</sup>iMesh and MOAB both define adjacencies using the topological concept of closure.  Since the closure of an entity includes the entity itself, the d-dimensional entities on the closure of a given entity should include the entity itself.

 \ref contents

  \section representation 8.Structured Mesh Representation

A structured mesh is defined as a D-dimensional mesh whose interior vertices have 2D connected edges.   Structured mesh can be stored without connectivity, if certain information is kept about the parametric space of each structured block of mesh.  MOAB can represent structured mesh with implicit connectivity, saving approximately 57% of the storage cost compared to an unstructured representation<sup>7</sup>.  Since connectivity must be computed on the fly, these queries execute a bit slower than those for unstructured mesh.  More information on the theory behind MOAB's structured mesh representation can be found in “MOAB-SD: Integrated structured and unstructured mesh representation”[17].

Currently, MOAB's structured mesh representation can only be used by creating structured mesh at runtime; that is, structured mesh is saved/restored in an unstructured format in MOAB's HDF5-based native save format.  For more details on how to use MOAB's structured mesh representation, see the scdseq_test.cpp source file in the test/ directory.

<sup>7</sup> This assumes vertex coordinates are represented explicitly, and that there are approximately the same number of vertices and hexahedra in a structured hex mesh.

 \ref contents

  \section element 9.Spectral Element Meshes

The Spectral Element Method (SEM) is a high-order method, using a polynomial Legendre interpolation basis with Gauss-Lobatto quadrature points, in contrast to the Lagrange basis used in (linear) finite elements [20].  SEM obtains exponential convergence with decreasing mesh characteristic sizes, and codes implementing this method typically have high floating-point intensity, making the method highly efficient on modern CPUs.  Most Nth-order SEM codes require tensor product cuboid (quad/hex) meshes, with each d-dimensional element containing (N+1)d degrees of freedom (DOFs).  There are various methods for representing SEM meshes and solution fields on them; this document discusses these methods and the tradeoffs between them.  The mesh parts of this discussion are given in terms of the iMesh mesh interface and its implementation by the MOAB mesh library [21].

The figure above shows a two-dimensional 3rd-order SEM mesh consisting of four quadrilaterals.  For this mesh, each quadrilateral has (N+1)^2=16 DOFs, with corner and edge degrees of freedom shared between neighboring quadrilaterals.

  \ref contents

  \subsection nineone 9.1. Representations

There are various representations of this mesh in a mesh database like MOAB, depending on how DOFs are related to mesh entities and tags on those entities.  We mention several possible representations:

-# Corner vertices, element-based DOFs: Each quadrilateral is defined by four vertices, ordered in CCW order typical of FE meshes.  DOFs are stored as tags on quadrilaterals, with size (N+1)^2 values, ordered lexicographically (i.e. as a 2D array tag(i,j) with i varying faster than j.)  In the figure above, the connectivity for face 1 would be (1, 4, 16, 13), and DOFs would be ordered (1..16).  Note that in this representation, tag values for DOFs shared by neighboring elements must be set multiple times, since there are as many copies of these DOFs as elements sharing them.
-# High-order FE-like elements: Each DOF is represented by a mesh vertex. Quadrilaterals each have (N+1)^2 vertices, ordered as they would be for high-order finite elements (corner vertices first, then mid-edge and mid-face elements; see [22]).  Mid -face, -edge, and -region vertices for a given edge/face/region would be ordered lexicographically, according to positive direction in a corresponding reference element.  In the figure above, the connectivity array for face 1 would be (1, 4, 16, 13, 2, 3, 8, 12, 14, 15, 5, 9, 6, 7, 10, 11).  DOF values are stored as tags on vertices.  Since DOFs are uniquely associated with vertices and vertices are shared by neighboring elements, tag values only need to be set once.  Full vertex-quadrilateral adjacencies are available, for all vertices.
-# Linear FE-like elements, one vertex per DOF, array with DOF vertices: Each quadrilateral is defined by four (corner) vertices, with additional vertices representing mid-edge and mid-face DOFs.  An additional “DOF array” tag is assigned to each quadrilateral, storing the array of vertices representing the (N+1)^2 DOFs for the quadrilateral, ordered lexicographically.  For the figure above, the connectivity array for face 1 would be (1, 4, 16, 13), and the DOF array would be (1..16), assuming that vertex handles are integers as shown in the figure.  DOF values are stored as tags on vertices, and lexicographically-ordered arrays of DOFs can be retrieved using the DOF array tag as input to the tag_get_data function in MOAB.  Adjacency functions would only be meaningful for corner vertices, but tag values would only need to be set once per DOF.
-# High-order FE-like elements, array with DOF vertices: This is a combination of options 2 and 3.  The advantage would be full vertex-quad adjacency support and direct availability of lexicographically-ordered vertex arrays, at the expense of more memory.
-# Convert to linear mesh: Since a spectral element is a cuboid with higher-order vertices, it can always be converted to N^2 linear cuboids using the high-order vertices as corners of the finer quads/hexes.  This is how readers in ParaView and VisIt typically import spectral meshes (CAM-SE also exports connectivity in this form).

As a convenience for applications, functions could also be provided for important tasks, like assembling the vertex handles for an entity in lexographic order (useful for option 2 above), and getting an array of tag values in lexicographic order (for option 3 above).

  \ref contents

  \subsection ninetwo 9.2. Tradeoffs

There are various competing tradeoffs in the various representation types.  These include:

- Adjacencies: being able to retrieve the element(s) using a given (corner or higher-order) vertex.
- Connectivity list: being able to retrieve the connectivity of a given element, consisting of all (corner + higher-order) vertices in the element, usually in lexicographical order.  This is closely linked with being able to access the connectivity list as a const*, i.e. using the list straight from memory without needing to copy it.
- Memory vs. time: There is a memory vs. execution time tradeoff between duplicating interface vertex solution/tag variables in neighboring elements (more memory but more time-efficient and allows direct access to tag storage by applications) versus using vertex-based tags (less memory but requires assembly of variables into lexicographically-ordered arrays, and prevents direct access from applications).
.

The lower-memory option (storing variables on vertices and assembling into lexicographically-ordered arrays for application use) usually ends up costing more in memory anyway, since applications must allocate their own storage for these arrays.  On the other hand, certain applications will always choose to do that, instead of sharing storage with MOAB for these variables.  In the case where applications do share memory with MOAB, other tools would need to interpret the lexicographically-ordered field arrays specially, instead of simply treating the vertex tags as a point-based field.

  \ref contents

  \subsection ninethree 9.3. MOAB Representation
In choosing the right MOAB representation for spectral meshes, we are trying to balance a) minimal memory usage, b) access to properly-ordered and -aligned tag storage, and c) maximal compatibility with tools likely to use MOAB.  The solution we propose is to use a representation most like option 2) above, with a few optional behaviors based on application requirements.  

In brief, we propose to represent elements using the linear, FE-ordered connectivity list (containing only corner vertices from the spectral element), with field variables written to either vertices, lexicographically-ordered arrays on elements, or both, and with a lexicographically-ordered array (stored on tag SPECTRAL_VERTICES) of all (corner+higher-order) vertices stored on elements.  In the either/or case, the choice will be evident from the tag size and the entities on which the tag is set.  In the both case, the tag name will have a “-LEX” suffix for the element tags, and the size of the element tag will be (N+1)^2 times that of the vertex-based tag.  Finally, the file set containing the spectral elements (or the root set, if no file set was input to the read) will contain a “SPECTRAL_ORDER” tag whose value is N.  These conventions are described in the “Metadata Information” document distributed with the MOAB source code.

  \ref contents

  \section performance 10.Performance and Using MOAB Efficiently from Applications

MOAB is designed to operate efficiently on groups of entities and for large meshes.  Applications will be most efficient when they operate on entities in groups, especially groups which are close in their order of creation.  The MOAB API is structured to encourage operations on groups of entities.  Conversely, MOAB will not perform as well as other libraries if there are frequent deletion and creation of entities.  For those types of applications, a mesh library using a C++ object-based representation is more appropriate.  In this section, performance of MOAB when executing a variety of tasks is described, and compared to that of other representations.  Of course, these metrics are based on the particular models and environments where they are run, and may or may not be representative of other application types.

One useful measure of MOAB performance is in the representation and query of a large mesh.  MOAB includes a performance test, located in the test/perf directory, in which a single rectangular region of hexahedral elements is created then queried; the following steps are performed:

- Create the vertices and hexes in the mesh
- For each vertex, get the set of connected hexahedra
- For each hex, get the connected vertices, their coordinates, average them, and assign them as a tag on the hexes
.

This test can be run on your system to determine the runtime and memory performance for these queries in MOAB.

  \ref contents

  \section error-handling 11.Error Handling

Errors are handled through the routine MBError(). This routine calls MBTraceBackErrorHandler(), the default error handler which tries to print a traceback.

The arguments to MBTraceBackErrorHandler() are the line number where the error occurred, the function where error was detected, the file in which
the error was detected, the corresponding directory, the error message, and the error type.

A small set of macros is used to make the error handling lightweight. These macros are used throughout
the MOAB libraries and can be employed by the application programmer as well. When an error is first
detected, one should set it by calling\n
\code
MB_SET_ERR(err_code, err_msg);
\endcode
Note, err_msg can be a string literal, or a C++ style output stream with << operators, such that the error message string is formated, like\n
\code
MB_SET_ERR(MB_FAILURE, "Failed " << n << " times");
\endcode

The user should check the return codes for all MOAB routines (and possibly user-defined routines as well) with\n
\code
ErrorCode rval = MOABRoutine(...);MB_CHK_ERR(rval);
\endcode
To pass back a new error message (if rval is not MB_SUCCESS), use\n
\code
ErrorCode rval = MOABRoutine(...);MB_CHK_SET_ERR(rval, "User specified error message string (or stream)");
\endcode
If this procedure is followed throughout all of the user’s libraries and codes, any error will by default generate
a clean traceback of the location of the error.

In addition to the basic macros mentioned above, there are some variations, such as (for more information, refer to src/moab/ErrorHandler.hpp):
- MB_SET_GLB_ERR() to set a globally fatal error (for all processors)
- MB_SET_ERR_RET() for functions that return void type
- MB_SET_ERR_RET_VAL() for functions that return any data type
- MB_SET_ERR_CONT() to continue execution instead of returning from current function

The error control mechanism are enabled by default if a MOAB Core instance is created. Otherwise, the user needs to call
MBErrorHandler_Init() and MBErrorHandler_Finalize() at the application level in the main function.
For example code on error handling, please refer to examples/TestErrorHandling.cpp, examples/TestErrorHandlingPar.cpp and examples/ErrorHandlingSimulation.cpp.

  \ref contents

  \section conclusions 12.Conclusions and Future Plans

MOAB, a Mesh-Oriented datABase, provides a simple but powerful data abstraction to structured and unstructured mesh, and makes that abstraction available through a function API.  MOAB provides the mesh representation for the VERDE mesh verification tool, which demonstrates some of the powerful mesh metadata representation capabilities in MOAB.  MOAB includes modules that import mesh in the ExodusII, CUBIT .cub and Vtk file formats, as well as the capability to write mesh to ExodusII, all without licensing restrictions normally found in ExodusII-based applications.  MOAB also has the capability to represent and query structured mesh in a way that optimizes storage space using the parametric space of a structured mesh; see Ref. [17] for details.

Initial results have demonstrated that the data abstraction provided by MOAB is powerful enough to represent many different kinds of mesh data found in real applications, including geometric topology groupings and relations, boundary condition groupings, and inter-processor interface representation.  Our future plans are to further explore how these abstractions can be used in the design through analysis process.

  \ref contents

  \section references 13.References

[1]	M. Fatenejad and G.A. Moses, “Cooper radiation hydrodynamics code..”

[2]	T.J. Tautges and A. Caceres, “Scalable parallel solution coupling for multiphysics reactor simulation,” Journal of Physics: Conference Series,  vol. 180, 2009.

[3]	T.J. Tautges, MOAB Meta-Data Information, 2010.

[4]	T.J. Tautges, “MOAB - ITAPS – Trac.”, http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB

[5]	“MOAB Developers Email List.”, moab-dev@mcs.anl.gov.

[6]	“MOAB Users Email List.”, moab@mcs.anl.gov.

[7]	“MOAB online documentation.”, http://ftp.mcs.anl.gov/pub/fathom/moab-docs/index.html

[8]	T.J. Tautges, “Canonical numbering systems for finite-element codes,” Communications in Numerical Methods in Engineering,  vol. Online, Mar. 2009.

[9]	L.A. Schoof and V.R. Yarberry, EXODUS II: A Finite Element Data Model,  Albuquerque, NM: Sandia National Laboratories, 1994.

[10]	M. PATRAN, “PATRAN User’s Manual,” 2005.

[11]	VisIt User's Guide.

[12]	K. Devine, E. Boman, R. Heaphy, B. Hendrickson, and C. Vaughan, “Zoltan Data Management Services for Parallel Dynamic Applications,” Computing in Science and Engineering,  vol. 4, 2002, pp. 90–97.

[13]	T.J. Tautges, P.P.H. Wilson, J. Kraftcheck, B.F. Smith, and D.L. Henderson, “Acceleration Techniques for Direct Use of  CAD-Based Geometries in Monte Carlo Radiation Transport,” International Conference on Mathematics, Computational Methods & Reactor Physics (M&C 2009),  Saratoga Springs, NY: American Nuclear Society, 2009.

[14]	H. Kim and T. Tautges, “EBMesh: An Embedded Boundary Meshing Tool,” in preparation.

[15]	G.D. Sjaardema, T.J. Tautges, T.J. Wilson, S.J. Owen, T.D. Blacker, W.J. Bohnhoff, T.L. Edwards, J.R. Hipp, R.R. Lober, and S.A. Mitchell, CUBIT mesh generation environment Volume 1: Users manual, Sandia National Laboratories, May 1994, 1994.

[16]	T.J. Tautges, “CGM: A geometry interface for mesh generation, analysis and other applications,” Engineering with Computers,  vol. 17, 2001, pp. 299-314.

[17]	T. J. Tautges, MOAB-SD: Integrated structured and unstructured mesh representation, Engineering with Computers, vol. 20, no. 3, pp. 286-293, 2004.

[18]	“Interoperable Technologies for Advanced Petascale Simulations (ITAPS),” Interoperable Technologies for Advanced Petascale Simulations (ITAPS).

[19]	P. Knupp, “Mesh quality improvement for SciDAC applications,” Journal of Physics: Conference Series,  vol. 46, 2006, pp. 458-462.

[20]	M. O. Deville, P. F. Fischer, and E. H. Mund, High-order methods for incompressible fluid flow. Cambridge, UK; New York: Cambridge University Press, 2002.

[21]	T. J. Tautges, “MOAB Wiki.” [Online]. Available: http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB. [Accessed: 30-Oct-2012].

[22]	T. J. Tautges, “Canonical numbering systems for finite-element codes,” International Journal for Numerical Methods in Biomedical Engineering, vol. 26, no. 12, pp. 1559–1572, 2010.

[23]    V. Dyedov, N. Ray, D.Einstein, X. Jiao, T. Tautges, “AHF: Array-based half-facet data structure for mixed-dimensional and non-manifold meshes”, In Proceedings of 22nd International Meshing Roundtable, Orlando, Florida, October 2013.


  \ref contents

  \page differences Differences Between iMesh and MOAB

  The data models used in MOAB and iMesh are quite similar, but not identical.The most significant differences are the following:

- Tags: MOAB differentiates between DENSE, SPARSE, and BIT tags, using different storage models for each, while iMesh uses a single tag concept.  iMesh allows application to query whether an entity has been given a tag of a specified type; this query is incompatible with the concept of a DENSE tag with a default value.  Thus, MOAB’s iMesh implementation creates SPARSE tags by default, and tags created and accessed through this interface will use more memory than DENSE tags created through MOAB’s native interface.  To mitigate this problem, MOAB implements an extension of the iMesh_createTag function which allows specification of the tag type (DENSE, SPARSE, etc.) to be created.  See later in this section for more information.

- Higher-order nodes: ITAPS currently handles higher-order elements (e.g. a 10-node tetrahedron) using a special “Shape” interface.  In this interface, higher-order nodes are only accessible through the AEntities which they resolve.  MOAB’s iMesh implementation provides access to higher-order nodes in the same manner described in Section  , by varying the number of vertices defining each entity.  As a result, if higher-order entities are used in a model, the functions returning connectivity and vertex adjacencies always return all vertices, rather than providing an option to return just corner vertices.

- Self-adjacencies: iMesh specifies that entities are not self-adjacent; that is, requesting adjacencies of the same dimension/type results in an error.  MOAB does not consider this an error, returning the entity itself.

- Adjacency table and AEntities: iMesh uses the concept of an “adjacency table” to determine which AEntities are available and created by default.  MOAB uses input arguments to the get_adjacencies functions to control whether AEntities are created.  These flags provide finer-grained control over AEntities, but make it slightly less convenient to ensure that AEntities of a given dimension are always created.
.

 \page figures List of Figures
 
  This page is intended to be empty.
 
  \page tables List of Tables
 
     \ref tableone

     \ref tabletwo

     \ref tablethree

     \ref tablefour

     \ref tablefive

     \ref tablesix

     \ref tableseven
 
  \page building Building & Installing
 
  MOAB uses an autoconf and libtool-based build process by default.  The procedure used to build MOAB from scratch depends on whether the source code was obtained from a “tarball” or directly from the Subversion repository.  Assuming the latter, the following steps should be executed for building and installing MOAB:
  - Locate and build any required dependencies.  MOAB can be built with no dependencies on other libraries; this may be useful for applications only needing basic mesh representation and not needing to export mesh to formats implemented in other libraries.  MOAB’s native save/restore capability is built on HDF5-based files; applications needing to save and restore files from MOAB reliably should use this library.  MOAB also uses ExodusII, a netCDF-based file format developed at Sandia National Laboratories [10].  Applications needing to execute these tests should also build netCDF.  Note that MOAB uses netCDF’s C++ interface, which is not enabled by default in netCDF but can be enabled using the “–enable-cxx” option to netCDF’s configure script.
  - Unpack source code into <moab>, and change current working directory to that location.
  - Execute “autoreconf –fi”.
  - Run configure script, by executing “./configure <options>”.  Recommended options:
       -# –prefix=<install_dir>: directory below which MOAB library and include files will be installed; can either be the directory used for MOAB source (<moab> from step 1), or a different directory.
       -# –hdf5-dir=<hdf5_dir>: directory whose “include” and “lib” subdirectories hold HDF5 include and library, respectively.  MOAB uses HDF5 for its native save/restore format (see Section 4.6.1).
       -# –netcdf-dir=: directory whose “include” and “lib” subdirectories hold netCDF include and library, respectively.  MOAB uses netCDF-based files for many of its build tests.  If the location of netCDF cannot be found, MOAB’s build tests will not function properly, but MOAB will still be usable.
       .
  - Run “make check”; this runs a series of build tests, to verify that the MOAB build was successful.  Note this check will fail if netCDF is not used, but MOAB itself will still be usable from applications.
  - Run “make install”; this copies include files and libraries to subdirectories of the directory specified in the “prefix” option.
  .

These steps are sufficient for building MOAB against HDF5 and netCDF.  By default, a small number of standard MOAB-based applications are also built, including mbconvert (a utility for reading and writing files), mbsize (for querying basic information about a mesh), and the iMesh interface (see Section 7).  Other utilities can be enabled using various other options to the configure script; for a complete list of build options, execute “./configure –help”.
 
 */
