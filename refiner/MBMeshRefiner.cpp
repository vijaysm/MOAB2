#include "MBMeshRefiner.hpp"

#include "MBEdgeSizeEvaluator.hpp"
#include "MBEntityRefiner.hpp"
#include "MBInterface.hpp"
#include "MBRefinerTagManager.hpp"
#include "MBMeshOutputFunctor.hpp"

#ifdef USE_MPI
#include "MBParallelComm.hpp"
#include "MBmpi.h"
#else // USE_MPI
typedef int MPI_Comm;
#endif // USE_MPI

/**\brief Construct a mesh refiner.
  * The input and output mesh pointers may be identical.
  * Existing elements will <b>not</b> be removed from the input mesh
  * as they are refined, so the adjacencies for entitities may appear
  * strange after refinement.
  */
MBMeshRefiner::MBMeshRefiner( MBInterface* imesh, MBInterface* omesh )
{  
  this->mesh_in = imesh;
  this->mesh_out = omesh;
  this->tag_manager = new MBRefinerTagManager( this->mesh_in, this->mesh_out );
  this->output_functor = new MBMeshOutputFunctor( this->tag_manager );
  this->entity_refiner = 0;
  this->comm = MBParallelComm::get_pcomm( this->mesh_out, 0 );
}

/**\brief Destroy a mesh refiner.
  *
  * Note that any MBEntityRefiner object associated with the mesh refiner will be deleted inside this destructor.
  * Destruction is virtual so subclasses may clean up after refinement.
  */
MBMeshRefiner::~MBMeshRefiner()
{
  delete this->tag_manager;
  delete this->output_functor;
  if ( this->entity_refiner )
    delete this->entity_refiner;
}

/**\brief Specify which techniqe will be used to refine entities in the mesh.
  * The entity refiner object you specify is ``owned'' by the mesh refiner after this call;
  * the entity refiner will be deleted when this mesh refiner is destroyed.
  * You should not delete the entity refiner yourself.
  */
bool MBMeshRefiner::set_entity_refiner( MBEntityRefiner* er )
{
  if ( ! er || er == this->entity_refiner )
    return false;

  this->entity_refiner = er;
  return true;
}

/**\brief A convenience method to reset the list of tags to be copied to output vertices.
  * This simply calls the method of the same name on the tag manager.
  */
void MBMeshRefiner::reset_vertex_tags()
{
  this->tag_manager->reset_vertex_tags();
}

/**\brief A convenience method to add a tag to be copied/interpolated from input vertices to output vertices.
  * This simply calls the method of the same name on the tag manager.
  */
int MBMeshRefiner::add_vertex_tag( MBTag tag_handle )
{
  return this->tag_manager->add_vertex_tag( tag_handle );
}

struct MBMeshRefinerIterator {
  MBRange subset;
  MBEntityHandle destination_set;
};

/**\brief Refine entities in a mesh set.
  * This will recursively descend any mesh sets contained in the \a range.
  * It returns false when not able to refine (because no entity refiner is
  * set or no edge size evaluator has been set on the entity refiner) and
  * true otherwise.
  */
bool MBMeshRefiner::refine( MBRange& range )
{
  this->tag_manager->create_output_tags();
  if ( ! this->entity_refiner->prepare( this->tag_manager, this->output_functor ) )
    { // Oops, perhaps the edge_size_evaluator was not set?
    return false;
    }

  MBMeshRefinerIterator entry;
  std::vector<MBMeshRefinerIterator> work;

  entry.subset = range;
  entry.destination_set = 0;
  work.push_back( entry );

  while ( ! work.empty() )
    {
    entry = work.back();
    work.pop_back();
    this->output_functor->destination_set = entry.destination_set;
    for ( MBRange::const_iterator it = entry.subset.begin(); it != entry.subset.end(); ++ it )
      {
      MBEntityType etyp = this->mesh_in->type_from_handle( *it );
      if ( etyp == MBENTITYSET )
        {
        MBRange set_ents;
        if ( this->mesh_in->get_entities_by_handle( *it, set_ents, false ) == MB_SUCCESS )
          {
          // Create a matching set on the output mesh.
          MBMeshRefinerIterator set_work;
          unsigned int set_work_opts;
          this->mesh_in->get_meshset_options( *it, set_work_opts );
          this->mesh_out->create_meshset( set_work_opts, set_work.destination_set );
          set_work.subset = set_ents;
          work.push_back( set_work );
          // Copy any per-element tag values the user has requested to the output set.
          this->tag_manager->set_element_tags_from_ent( *it );
          this->tag_manager->assign_element_tags( set_work.destination_set );
          // Copy the global ID to the new set (assuming it exists).
          this->tag_manager->copy_gid( *it, set_work.destination_set );
          }
        }
      else
        {
        this->tag_manager->set_element_tags_from_ent( *it );
        this->tag_manager->set_element_procs_from_ent( *it );
        this->entity_refiner->refine_entity( etyp, *it );
        }
      }
    }
  this->output_functor->assign_global_ids( this->comm );

  return true;
}

