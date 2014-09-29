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
     introduce a class to the top-level MOAB namespace, don't choose
     an overly genereric name like \c Point3D .
 - Constants and Macros
   - Don't use a pre-processor macro where a const variable or an inline or
     template function will suffice.
     There is absolutely benefit to the former over the later with modern 
     compilers.  Further, using  macros bypasses typechecking that the compiler
     would otherwise do for you and if used in headers, introduce names into
     the global rather than MOAB namespace.
   - Don't define constants that are already provided by standard libraries.
     For example, use \c M_PI as defined in \c math.h rather than defining
     your own constant.
\section commits Making Repository Commits
As a general rule, developers should update frequently, and commit changes often.  However, the repository should always remain
in a state where the code can be compiled.  Most of the time, the code should also successfully execute "make check" run from the
top-level directory.  If you commit code that violates this principal, it should be your first priority to return the repository
code to a compilable state, and your second priority to make sure "make check" runs without errors.

Commits to the repository should also come with a non-trivial, useful, non-verbose log message.  Oftentimes the best way to generate
this message is to run 'commit -a', and include a comment on 
each file that changed, then Ctrl+O to write out, followed by 'Enter' and Ctrl+X.  Many times it is helpful to state that 'make check runs successfully' at the end of the log message.
Although it would be possible and many software projects do it, we prefer not to force successful execution of the test suite 
before every commit.  Developers should make every effort to avoid having to impose this constraint, by running a make check
before every commit.

\section git Git Repository Practices
As most of our code repositories uses git as the revision control system, it is important to decide on a workflow that can be followed by the individual developer. The way that any individual developer interact with the upstream git repository can have an important impact on other developers and the ability to identify and manage individual changes.  This set of guidelines and practices attempts to establish some standards for how developers will interact with the upstream git repository.
The Atlassian website <a href="https://www.atlassian.com/git/workflows">describes a number of git workflows</a> , and provides good background reading for a few standard models of such workflows.  More important than choosing any one of these workflows precisely are the adoption of some of the concepts embodied within them and understanding the implications of those concepts.

\subsection outside-master Working Outside the Master Branch
A critical concept is that all changes shall be developed outside of the master<sup>1</sup> branch.  Whether they are in a different branch of the upstream<sup>2</sup> repository (gitflow) or a branch of an entirely different fork (forking workflow) is secondary.  This is a well-established concept regardless of the workflow being adopted, and allows a number of other benefits as described below.

\subsubsection fork Working on a Different Fork
There are a number of benefits of working on a different fork rather than a branch of the upstream repo, although not strictly technical:
- Developers, particularly new developers, are liberated from the constant oversight of others as they explore new code options.  The impact of this may depend on an individual developer’s personality, but for some it creates a refuge where they can be more free and creative.
- Similarly, assuming that all changesets in the upstream repo are communicated to the entire development team, the team is spared a noisy stream of notifications and can focus their attention on the rarer occurrence of a pull request notification.

\subsubsection pr All Changes are Committed by Pull Request
Although this can be imposed technically by limiting the authority to alter the upstream repo (as in the forking workflow), a healthy developer community can also simply rely on convention.  The advantage of doing it by convention rather than by restriction is that it is easier to distribute the load of reviewing and accepting changes.
A critical consequence of this decision is that all code is reviewed before it is committed to the upstream master branch.  This has benefits to overall quality in two related ways:
- the code under review will improve due to the review itself, and
- those involved in the review will maintain a broad awareness of the code base resulting in better contributions from them.

This practice does, however, place a substantial burden on the developers to perform timely reviews of the pull requested (PR’ed) code.  PR’s that languish without sufficient review have a number of negative consequences:
- they need to be refreshed simply to keep them up-to-date with the possibly advancing upstream/master
- they may delay further development on similar or related features
- they breed frustration in the original developer, undermining the community as a whole.
Bitbucket provides powerful collaboration tools that greatly facilitate this process.

<sup>1</sup> Although a repository may choose a different name for its main development branch, this document will refer to that as the “master” branch.

<sup>2</sup> For this discussion, the “upstream” repo will refer to the centralized authoritative repository used to synchronize changes.

\subsection Some Git Mechanics to Keep it Clean
Given the above practices, there are some mechanical details that can help ensure that the upstream/master repository is always in a state that facilitates all repository actions and interactions.

-# Feature branches being used for development should be kept up-to-date with the upstream/master by rebase only.  When a feature branch is rebased against the upstream/master, all changes in the upstream/master are inserted into the feature branch at a point in its history that is prior to any of the changes of the feature branch.  This can require conflict resultion as the feature branch changes are “replayed” on top of the new upstream/master in its current state.  The primary advantage of this policy is that it keeps all of the feature branch changes contiguous.  If, by contrast, the upstream/master is merged into the feature branch, the recent changes in the upstream/master become woven into the set of changes in the feature branch.  This can make it more difficult to isolate changes later on.

    Strict adoption of this practice is important since a single merge into a feature branch that is then merged back into the upstream/master can make it nearly impossible for others to rebase.

    A typical workflow with pull-request might look like this, all using the command-line, except for submitting the final pull request.  Note that there is never a merge operation.
    -# synchronize your local `master` branch before anything else
    \code
     %> git checkout master
     %> git fetch upstream
     %> git rebase upstream/master
     \endcode

    -# now create a new feature branch from master
    \code
     %> git checkout -b my_feature_branch master
    \endcode

    -# now make changes, editing A.cpp, B.hpp, C.cpp

    -# now add/commit your changes to your local feature branch
    \code
     %> git add A.cpp B.hpp C.cpp
     %> git commit -m “Make sure you have a good commit message”
    \endcode
    -# push your changes to your feature branch on your fork (often called `origin`)
    \code
     %> git push origin my_feature_branch
    \endcode
    -# make more changes, editing B.hpp, D.hpp, E.cpp

    -# add/commit your changes to your local feature branch
    \code
    %> git add B.hpp D.hpp E.cpp
    %> git commit -m “Be sure you have another good commit message”
    \endcode
    -# push your changes to your freature branch on your fork (often called `origin`)
    \code
    %> git push origin my_feature_ranch
    \endcode
    -# When you are ready to submit a pull request, be sure that your feature branch is up-to-date. This first step may seem redundant but is here to be clear which branch we are acting on
    \code
    %> git checkout my_feature_branch
    %> git fetch upstream
    %> git rebase upstream/master
    \endcode
      This may generate conflicts that can be addressed at this point.

      NOTE: This step can be performed at any time and should be performed as often as practical to reduce the scope of potential conflicts.

    -# push your updated feature branch on your fork (often called `origin`)
     \code
     %> git push origin my_feature_branch
     \endcode
      This may require the ‘-f’ option to force the push.  (It is frequently necessary to force this push because the act of rebasing will “replay” the commits from the feature branch on top of the master, leading to different commit hashes.  Each of the commits will contain the same actual information, but because it has a different set of hashes, git will think there is an inconsistency and ask you to force the change.)

    -# Submit a pull request on Bitbucket, from your fork to the fathomteam fork.


-# When ready to be adopted into the upstream/master, feature branches should be combined by merge only.  This adds the changeset to the end of the upstream/master as a set of individual commits but in a contiguous block.

   A typical workflow to merge a pull-request might look like this, all using the command-line.
   -# synchronize your local `master` branch before anything else (just because it’s never a bad idea!)
   \code
   %> git checkout master
   %> git fetch upstream
   %> git rebase upstream/master
   \endcode
   -# add a remote for the user with the pull-request, perhaps the user is ‘other_user’
   \code
   %> git remote add other_user \
         git@bitbucket.org:other_user/moab.git
   \endcode
   -# fetch the other users repo
   \code
   %> git fetch other_user
   \endcode
   -# check out their feature branch
   \code
   %> git checkout -b pr_feature_branch \
       other_user/feature_branch
   \endcode
   -# confirm that it is up-to-date with the master. This first step may seem redundant but is here to be clear which branch we are acting on
   \code
   %> git checkout pr_feature_branch
   %> git fetch upstream
   %> git rebase upstream/master
   \endcode
   This may generate conflicts that can be addressed at this point.  You may want to request the original author (other_user) take care of these.
   -# once confirmed that it’s up-to-date with master, review this branch including:
      -reading the code
      -building the code
      -running tests
   -# once satisfied that the code meets the necessary standards and that all required/requested changes are fully incorporated into other_users’s feature branch, merge it into master
    \code
    %> git checkout master
    \endcode
    The next two steps may seem redundant but provide some QA
    \code
    %> git fetch upstream
    %> git rebase upstream/master
    %> git merge other_user/feature_branch
    \endcode
   -# push those changes into the master branch on bitbucket
    \code
    %> git push upstream/master
    \endcode

-# When a pull request is open for review, any changes to the feature branch will automatically update the pull request.  This is the appropriate way for a developer to respond to requests for changes that occur through the PR process.


-# If a developer has ongoing work that is based on a feature branch that is under consideration in an open PR, a new feature branch (B) should be created that is based on the previous feature branch (A).  Moreover, as changes are made to the original feature branch (A) due to the review process, the new feature branch (B) should be kept up-to-date by rebase against feature branch (A).  This keeps all subsequent changes of (B) downstream from the changes in (A).  Once feature branch (A) has been adopted into the upstream/master, the new feature branch (B) can start being rebased against the upstream/master instead.


-# When a repo is forked, its branches are not automatically synchronized with the corresponding branches on the upstream repo.  This requires a manual process of synchronization via a local clone.  Assuming that the local repo’s branch has the same name as the upstream branch (<branch>), and that the fork is known as “origin”:
    \code
     %> git fetch upstream
     %> git checkout <branch>
     %> git rebase upstream/<branch>
     %> git push origin <branch>
    \endcode
The decision of which branches to keep up-to-date is up to the developers.  Developers may choose to delete some branches from their own fork to avoid (a) the need to update it and (b) accidentally assuming that it is up-to-date.


-# When rebasing, it is not uncommon to encounter conflicts.  This will interrupt the rebasing process, and each conflicted file will contain conflict blocks.  You will be offered three choices:
 - manually resolve the conflicts, add the resolved files with git add,  and then git rebase --continue (do not commit the resolved files!)
 - abort the rebasing process entirely with git rebase --abort
 - skip the commit that causes the conflict, assuming that you are sure that this is the right thing to do, with git rebase --skip


-# Bitbucket offers a number of buttons/tools to manage changes between branches and forks.  Some of these operate in ways that are contrary to the practices recommended here, and others are consistent with these practices.  In general, it is best to know how to do all of those operations with the command-line instead of relying on Bitbucket, as it gives you full control over important details.


-# During code development, it might be necessary to work on the same branch on different machines. The workflow to update the local branch is to first fetch the remote changes and then perform a hard reset.
     \code
     %> git fetch origin
     %> git reset --hard origin/branch_name
     \endcode
One should be careful with the branch name as a hard reset would overwrite all changes in the working directory.


Top: \ref index 
  
 */
