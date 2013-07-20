/*!
  \page styleguide Coding Style Guide
Code developed in %MOAB should follow the coding styles described here.  Any deviations from this style
guide will result in severe berating and other verbal abuse.

\section dirstructure MOAB Directory Structure
%MOAB source code is organized in the following directory structure: \n
 - doc: Documentation is put here, along with the input file for Doxygen.  Most %MOAB documentation is doxygen-processed.
 - examples: Examples of %MOAB usage, both small and large.  These programs are not meant to be used as unit tests, but
     rather as further documentation on %MOAB usage.
 - src : Mesh related source codes. It includes:
   - io: %MOAB Input/Output classes.
   - moab: %MOAB core classes.
   - lotte: Computational Meshing basics.
   - parallel: Parallel mesh computation, i/o data processing methods.
 - test: All unit test programs should go below this directory.
   Please put the unit tests into their related subdirectories based on the test's
   purpose if possible.
If you're designing a new class or other code for %MOAB and are not sure where to put it, try to find something similar
and put it there.  Otherwise, email the %MOAB email list for pointers.  <em> In general, you should not need to create new
subdirectories in the %MOAB source code, except when implementing a new algorithm with more than about 2 files.</em>

\section sourcestyle Source Code Style and Best Practices
%MOAB code should abide by the following general rules:
 - Names:
   - Class names should be in the CamelBack style, e.g. EdgeMesh or VertexMesher.
   - Class member variables should be camelBack, e.g. EdgeMesh::schemeType; each member variable, e.g. int memberVariable, 
   should have set/get functions void member_variable(int newval) and int member_variable(), respectively.
   - Enumeration values should be all captitalized, with underscores avoided if possible (the enumeration name indicates
     the general purpose of the enumeration, so e.g. we use EQUAL, not EQUAL_MESH)
 - Source code should not contain tabs or MS-DOS newlines; tabs and other indentations should be set to a width of 2 spaces.
   For general tips on how to set your editor for this, see the %MOAB-dev discussion starting with <a href="https://lists.mcs.anl.gov/mailman/private/moab-dev/2011/000519.html">this message</a>.
 - Each class header should be fully commented; that includes:
   - A \\file comment block at the top of the file; DO NOT include things like Author and Date blocks; this stuff is available
     from subversion if we really need to know.
   - A \\class comment block, formatted like those in the %MOAB core classes.  THE FIELDS AFTER THE CLASS NAME ARE VERY IMPORTANT,
     as they tell developers how to include the class declaration in their code.  This information goes into the "Detailed
     Description" portion of the class documentation.  This block should include any features or limitations of the class.
     Eventually, we'll impose some standard boilerplate that each meshing class should use.
     Until then, please keep this block to around a paragraph.
   - Each function in both the public and private interfaces should be commented, INCLUDING ANY ARGUMENTS AND RETURN VALUES.
     See the %MOAB classes for examples of how to format these comments.  As a rule of thumb, your code should run through
     Doxygen without generating any warnings; in fact, Doxygen is sometimes helpful at pointing out inconsistencies in your
     class declaration.
 - Developers should avoid using \#include in header files, as they propagate dependencies more widely than necessary.  The only
   cases where other includes are needed are to import the declaration for a parent class, and to declare types used as
   non-pointer and non-reference function arguments.  In most cases, a forward-declaration statement (e.g. 'class MKCore') 
   will suffice.
 - Naming classes and other top-level constructs:
   - No names should be added to the global namespace.  Everything should be
     in the MOAB namespace.  An exception can be made for names with a static
     scope declared in a .cpp file, but class member functions never have a
     static scope.
   - Names should be kept as private as possible.  If declaring a struct or 
     utility class that is used internally by some other class, consider 
     defining it in the .cpp file of the main class or a separate header 
     only included in that .cpp file and using (if necessary) only forward
     delcarations (e.g. \c struct \c Point3D;) in the header file used
     by other code.  If that is not possible, then consider nesting the
     definitions such that the scope of the name is limited to that of the
     class using it. 
   - Any names introduced into the top-level MOAB namespace should be
     sufficiently unique to avoid conflicts with other code.  If you must 
     introduce a class to the top-level Meshkit namespace, don't choose
     an overly genereric name like \c Point3D .
 - Constants and Macros
   - Don't use a pre-processor macro where a const variable or an inline or
     template function will suffice.
     There is absolutely benefit to the former over the later with modern 
     compilers.  Further, using  macros bypasses typechecking that the compiler
     would otherwise do for you and if used in headers, introduce names into
     the global rather than Meshkit namespace.
   - Don't define constants that are already provided by standard libraries.
     For example, use \c M_PI as defined in \c math.h rather than defining
     your own constant.
\section commits Making Repository Commits
As a general rule, developers should update frequently, and commit changes often.  However, the repository should always remain
in a state where the code can be compiled.  Most of the time, the code should also successfully execute "make check" run from the
top-level directory.  If you commit code that violates this principal, it should be your first priority to return the repository
code to a compilable state, and your second priority to make sure "make check" runs without errors.

Commits to the repository should also come with a non-trivial, useful, non-verbose log message.  Oftentimes the best way to generate
this message is to run 'svn diff > diffs', and edit the diffs file to remove specific line changes but include a comment on 
each file that changed.  Many times it is helpful to state that 'make check runs successfully' at the end of the log message.
Although it would be possible and many software projects do it, we prefer not to force successful execution of the test suite 
before every commit.  Developers should make every effort to avoid having to impose this constraint, by running a make check
before every commit.

Top: \ref index 
  
 */
