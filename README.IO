This document describes file formats supported by MOAB and options
controlling file input and output.


Supported File Formats
----------------------

Some of the file formats listed below may not be supported by a particular
build of MOAB depending on the availability of external libraries.  An 
up-to-date list of file formats supported in a particular build of MOAB
can be obtained programatically using the MBReaderWriterSet API or as 
a simple list using the '-l' option of the mbconvert utility.

Format              Name    Read         Write   File name suffixes
------------------  ------  ------------ ------  --------------------
MOAB native         MOAB    yes          yes     h5m mhdf
Exodus II           EXODUS  yes          yes     exo exoII exo2 g gen
Kitware VTK         VTK     up to v3.0   v3.0    vtk
Cubit               CUBIT   yes          no      cub
SLAC                SLAC    no           yes     slac
GMV                 GMV     no           yes     gmv
Ansys               ANSYS   no           yes     ans
Gmsh                GMSH    v1.0, v2.0   v2.0    msh gmsh
Stereo Lithography  STL     yes          yes     stl

Any of the values from the 'Name' column may be passed as an option
to the file IO methods to request a particular file.  If no file 
format is specified, the default is to choose the write format using
the file extension and to try all file readers until one succeeds.


---------------
File IO Options
---------------

An options list as passed to MOAB file IO routines is a single C-style
string containing the concatenation of a list of string options, were
individual options are separated by a designated separator character.
The default separator character is a semicolon (;).  To specify an alternate
separator character, begin the options string with a semicolon followed
by the desired separator.  Options are not case sensitive.


---------------
Common Options
---------------

  PRECISION=<N>

Specify the precision to use when writing float and double values
(such as node coordinates) to text-based file formats.


  CREATE
  
Do not overwrite existing file.

-------------------
Parallel IO Options
-------------------

MOAB must be built with parallel support enabled before these options
can be used.  

Parallel Read Options:
----------------------

  PARALLEL={NONE|BCAST|BCAST_DELETE|READ_DELETE|READ_PARALLEL|FORMAT}

Set parallel read mode.  Options are:
  - NONE  - force serial read/write on each processor (default)
  - BCAST - read on one processor and broadcast a copy to all others
  - BCAST_DELETE - ?
  - READ_DELETE - read file on all processors and partition by deleting
                  mesh from non-local partitions.
  - READ_PARALLEL - ?
  - FORMAT- use format-specific parallel IO capabilities, if available

  
  PARTITION
  PARTITION=<tag_name>

Specify that mesh should be partitioned using partition IDs stored
in a tag.  If the tag name is not specified, the default ("PARTITION")
is used.

  PARTITION_VALS=<int_list>

Specify, for each processor, the list of partition IDs.

  PARTITION_DISTRIBUTE

?

  MPI_IO_RANK=<RANK>

For IO modes in which a single processor handles all disk access, the
MPI rank of the processor to use.  Default is 0.


-----------------------
Format-specific options
-----------------------

Stereo Lithography (STL) files
------------------------------
  
   BINARY|ASCII
   
Write binary or text STL file.  Default is text.

   BIG_ENDIAN|LITTLE_ENDIAN
   
Force byte ordering of binary data.  Default is BIG_ENDIAN for writing
and autodetection for reading (BIG_ENDIAN if autodetect fails).
