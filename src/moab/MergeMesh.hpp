#ifndef MERGEMESH_HPP
#define MERGEMESH_HPP

#include "moab/Interface.hpp"
#include "moab/Range.hpp"

namespace moab {

class MergeMesh 
{
public:
    /* \brief Constructor
     */
  MergeMesh(moab::Interface *mbImpl);
  
    /* \brief Destructor
     */
  virtual ~MergeMesh();

    /* \brief Merge vertices in elements passed in
     */
  moab::ErrorCode merge_entities(moab::EntityHandle *elems,
                      int elems_size,
                      const double merge_tol,
                      const int do_merge = true,
                      const int update_sets = false,
                      moab::Tag merge_tag = 0);

      //- perform the actual merge
  moab::ErrorCode perform_merge(moab::Tag merged_to);
private:
  //iMesh_Instance imeshImpl;

  double mergeTol, mergeTolSq;

  moab::Tag mergeTag;

    //- given a kdtree, set tag on vertices in leaf nodes with vertices
    //- to which they should be merged
  moab::ErrorCode find_merged_to(moab::EntityHandle &tree_root,
				 moab::Tag merged_to);
  
  moab::ErrorCode merge_entities(moab::Range &elems,
                             const int do_merge,
                             const int update_sets,
                             moab::Tag merge_tag);

  moab::Interface *mbImpl;

    //- the tag pointing to the entity to which an entity will be merged
  moab::Tag mbMergeTag;

    //- entities which will go away after the merge
  moab::Range deadEnts;
  
};

inline MergeMesh::MergeMesh(Interface *impl) 
  : mbImpl(impl)
{
}

inline MergeMesh::~MergeMesh() 
{
  if (mbMergeTag) mbImpl->tag_delete(mbMergeTag);
}

}

#endif

