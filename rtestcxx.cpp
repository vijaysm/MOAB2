/**
 * \file rtestcxx.cpp
 *
 * \a unit test for the TSTT association interface
 *
 */
#include "iRel.hh"
#include "iBase.hh"
#include "iMesh.hh"
#include "iGeom.hh"
#include "iBase_SNL_SIDL_defs.h"
#include "sidl_BaseInterface_IOR.h"
#include "MBTagConventions.hpp"
#include <iostream>
#include <assert.h>

const bool debug = true;

class RefEntity;
using namespace std;

typedef void* iBase_TagHandle;
typedef void* iRel_EntityHandle;

#define ARRAY_PTR(array, type) reinterpret_cast<type*>(array._get_ior()->d_firstElement)

#define CLASSIFICATION_TAG_NAME "MESH_GEOM_CLASSIFICATION"

/*!
  prints out a result string based on the value of error_code
*/
void handle_error_code(const bool result,
                       int &number_failed,
                       int &/*number_not_implemented*/,
                       int &number_successful)
{
  if (result) {
    cout << "Success";
    number_successful++;
  }
  else {
    cout << "Failure";
    number_failed++;
  }
}

#define CHECK_SIZE(array, size)  \
          if (array._is_nil() || ARRAY_SIZE(array) == 0) array = array.create1d(size); \
          else if (ARRAY_SIZE(array) < size) { \
             std::cerr << "   Array passed in is non-zero but too short." << std::cerr; \
           return false; }


bool print_geom_info(iGeom::Geometry &geom,
                     iRel_EntityHandle *gentities,
                     int gentities_size)
{
    // print information about this entity
  iGeom::EntityType ent_type;
  CAST_iGeom_INTERFACE(geom, geom_topo, Topology, false);
  CAST_iBase_INTERFACE(geom, geom_etag, EntTag, false);
  iBase_TagHandle tag_handle = 0;
  int gid = -1;
  try {
    tag_handle = geom_etag.getTagHandle("GLOBAL_ID");
  }
  catch (iBase::Error err) {
    return 0;
  }

  for (int g = 0; g < gentities_size; g++) {
    try {
      ent_type = geom_topo.getEntType(gentities[g]);
      gid = geom_etag.getIntData(gentities[g], tag_handle);
    }
    catch (iBase::Error err) {
      std::cerr << "Trouble getting entity adjacencies or types." << std::endl;
      std::cerr << "Description: " << err.getDescription() << std::endl;
      return 0;
    }

    const char *type_names[] = {"Vertex", "Edge", "Face", "Region"};
    printf("%s %d\n", type_names[ent_type], gid);
  }

  return 1;
}

bool print_mesh_info(iMesh::Mesh &mesh, iRel_EntityHandle *mentities,
                     int mentities_size, bool is_set)
{
    // print information about these entities
  CAST_iMesh_INTERFACE(mesh, mesh_ent, Entity, false);
  CAST_iMesh_INTERFACE(mesh, mesh_arr, Arr, false);
  CAST_iBase_INTERFACE(mesh, mesh_etag, EntTag, false);
  CAST_iBase_INTERFACE(mesh, mesh_stag, SetTag, false);
  CAST_iBase_INTERFACE(mesh, mesh_eset, EntSet, false);

  for (int m = 0; m < mentities_size; m++) {

      // get adjacencies first; assume not more than 50
    sidl::array<iRel_EntityHandle> adj_ents = adj_ents.create1d(50);
    sidl::array<iMesh::EntityType> ent_types;
    int adj_ents_size, ent_types_size = 0;
    if (!is_set) {
      try {
        mesh_ent.getEntAdj(mentities[m], iMesh::EntityType_ALL_TYPES,
                           adj_ents, adj_ents_size);

        assert(adj_ents_size < 50);

          // put this ent on the end, then get types
        adj_ents.set(adj_ents_size, mentities[m]);
        mesh_arr.getEntArrType(adj_ents, adj_ents_size+1,
                               ent_types, ent_types_size);
      }
      catch (iBase::Error err) {
        cerr << "Trouble getting entity adjacencies or types." << endl;
        cerr << "Description: " << err.getDescription() << endl;
        return false;
      }
    }

      // get tags on mentities[m]
    sidl::array<iBase_TagHandle> ment_tags;
    int ment_tags_size;
    std::vector<string> tag_names;
    if (!is_set) {
      try {
        mesh_etag.getAllTags(mentities[m], ment_tags, ment_tags_size);

          // while we're at it, get all the tag names
        for (int i = 0; i < ment_tags_size; i++)
          tag_names.push_back(mesh_etag.getTagName(ment_tags.get(i)));
      }
      catch (iBase::Error err) {
        cerr << "Trouble getting tags on an entity or their names." << endl;
        cerr << "Description: " << err.getDescription() << endl;
        return false;
      }
    }
    else {
      try {
        mesh_stag.getAllEntSetTags(mentities[m], ment_tags, ment_tags_size);

          // while we're at it, get all the tag names
        for (int i = 0; i < ment_tags_size; i++)
          tag_names.push_back(mesh_etag.getTagName(ment_tags.get(i)));
      }
      catch (iBase::Error err) {
        cerr << "Trouble getting tags on an entity or their names." << endl;
        cerr << "Description: " << err.getDescription() << endl;
        return false;
      }
    }

      // now print the information
    const char *type_names[] = {"Vertex", "Edge", "Face", "Region"};
    if (!is_set)
      cout << type_names[ent_types.get(ent_types_size)];
    else
      cout << "Set";

    cout << " " << (unsigned long)mentities[m] << ": " << endl;
    if (!is_set) {
      cout << "Adjacencies: ";
      for (int i = 0; i < adj_ents_size; i++) {
        if (i > 0) cout << ", ";
        cout << type_names[ent_types.get(ent_types_size)] << " " << (long)adj_ents.get(i);
      }
    }

    std::vector<std::string> opaque_tags;

    for (int i = 0; i < ment_tags_size; i++) {
      iBase::TagValueType tag_type;
      try {
        tag_type = mesh_etag.getTagType(ment_tags.get(i));
      }
      catch (iBase::Error) {
        cout << "(trouble getting type...)" << endl;
      }
      if (tag_type != iBase::TagValueType_BYTES)
        cout << tag_names[i] << " ";
      else opaque_tags.push_back(tag_names[i]);

      try {
        switch (tag_type) {
          case iBase::TagValueType_INTEGER:
            cout << "(Int value=";
            if (!is_set) cout << mesh_etag.getIntData(mentities[m], ment_tags.get(i)) << ")";
            else cout << mesh_stag.getEntSetIntData(mentities[m], ment_tags.get(i)) << ")";
            cout << endl;
            break;
          case iBase::TagValueType_DOUBLE:
            cout << "(Dbl value=";
            if (!is_set) cout << mesh_etag.getDblData(mentities[m], ment_tags.get(i)) << ")";
            else cout << mesh_stag.getEntSetDblData(mentities[m], ment_tags.get(i)) << ")";
            cout << endl;
            break;
          case iBase::TagValueType_ENTITY_HANDLE:
            cout << "(EH value=";
            if (!is_set) cout << mesh_etag.getEHData(mentities[m], ment_tags.get(i)) << ")";
            else cout << mesh_stag.getEntSetEHData(mentities[m], ment_tags.get(i)) << ")";
            cout << endl;
            break;
          default:
            break;
        }
      }
      catch (iBase::Error) {
        cout << "(trouble getting value...)";
      }
    }

    if (!opaque_tags.empty()) {
      cout << "Opaque tags: ";
      for (std::vector<std::string>::iterator vit = opaque_tags.begin();
           vit != opaque_tags.end(); vit++) {
        if (vit != opaque_tags.begin()) cout << ", ";
        cout << *vit;
      }
      cout << endl;
    }

    cout << endl;
  }

  return true;
}

bool print_entity_info(sidl::BaseInterface &iface,
                       iRel_EntityHandle *entities,
                       int entities_size, bool is_set)
{
  try {
    iMesh::Mesh mesh = iface;
    if (!mesh._is_nil())
      // if successful, it's mesh
      return print_mesh_info(mesh, entities, entities_size, is_set);
  }
  catch (iBase::Error err) {
  }

  try {
    iGeom::Geometry geom = iface;
    if (!geom._is_nil())
      // if successful, it's geom
    return print_geom_info(geom, entities, entities_size);
  }
  catch (iBase::Error err) {
  }

  return false;
}

bool get_gentities(iGeom::Geometry &geom,
                   const int low_dim, const int high_dim,
                   sidl::array<iRel_EntityHandle> &gents,
                   int &gents_size)
{
  assert (ARRAY_SIZE(gents) == 0);
  CAST_iGeom_INTERFACE(geom, geom_topo, Topology, false);

    // get all the gentities in the model

  std::vector<iRel_EntityHandle> gents_vec;
  for (int dim = high_dim; dim >= low_dim; dim--) {
    sidl::array<iRel_EntityHandle> these_gents;
    int these_gents_size = 0;
    try {
        // get all entities of this dimension
      geom_topo.getEntities(0, (iGeom::EntityType)dim,
                            these_gents, these_gents_size);
    }
    catch (iBase::Error err) {
      std::cerr << "Trouble getting gentities of dimension" << dim << std::endl;
      cerr << "Description: " << err.getDescription() << endl;
      return false;
    }
    std::copy(ARRAY_PTR(these_gents, iRel_EntityHandle),
              ARRAY_PTR(these_gents, iRel_EntityHandle)+these_gents_size,
              std::back_inserter(gents_vec));
  }

  CHECK_SIZE(gents, (int) gents_vec.size());
  std::copy(gents_vec.begin(), gents_vec.end(), ARRAY_PTR(gents, iRel_EntityHandle));
  gents_size = (int) gents_vec.size();

  return true;
}

bool project_to_geom(iGeom::Geometry &geom,
                     iMesh::Mesh &mesh,
                     iRel::Associate &lasso,
                     sidl::array<iRel_EntityHandle> &gents,
                     int &gents_size)
{
  CAST_iRel_INTERFACE(lasso, lasso_create, createEnt, false);
  for (int i = 0; i < gents_size; i++) {
    try {
      lasso_create.moveTo(geom, mesh, gents.get(i));
    }
    catch (iBase::Error err) {
      cerr << "ERROR : can not move meshes." << endl;
      cerr << "Description: " << err.getDescription() << endl;
      return false;
    }
  }

  return true;
}

/*!
  @test
  Load Mesh
  @li Load a geom and a mesh file
*/
bool load_geom_mesh_test(const std::string geom_filename,
                         const std::string mesh_filename,
                         iGeom::Geometry &geom,
                         iMesh::Mesh &mesh)
{
  CAST_iGeom_INTERFACE(geom, geom_cq, CoreQuery, false);

    // load a geom
    // load a mesh
  sidl::array<std::string> options;
  int options_size = 0;

  try {
    geom_cq.load(geom_filename, options, options_size);
  }
  catch (iBase::Error err) {
    cerr << "ERROR : can not load a geometry" << endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

    // load a mesh
  try {
    mesh.load(0, mesh_filename);
  }
  catch (iBase::Error err) {
    cerr << "ERROR : can not load a mesh" << endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

  if (debug) {
    cout << "Geometry: " << endl;
    sidl::array<iRel_EntityHandle> gentities;
    int gentities_size = 0;
    CAST_iGeom_INTERFACE(geom, geom_topo, Topology, false);
    try {
      geom_topo.getEntities(NULL, iGeom::EntityType_ALL_TYPES,
                            gentities, gentities_size);
    } catch (iBase::Error err) {
      return false;
    }

    print_geom_info(geom, ARRAY_PTR(gentities, iRel_EntityHandle),
                    gentities_size);

    cout << endl << endl << "Mesh: " << endl;
    sidl::array<iRel_EntityHandle> mentities;
    int mentities_size;
    CAST_iBase_INTERFACE(mesh, mesh_eset, EntSet, false);
    try {
      mesh_eset.getEntSets(NULL, 1, mentities, mentities_size);
    } catch (iBase::Error err) {
      return false;
    }
    print_mesh_info(mesh, ARRAY_PTR(mentities, iRel_EntityHandle),
                    mentities_size, true);
  }

  return true;
}

bool check_inferred_associations(iRel::Associate &lasso,
                                 sidl::BaseInterface iface1,
                                 sidl::BaseInterface iface2,
                                 sidl::array<iRel_EntityHandle> &ents1,
                                 int ents1_size,
                                 bool is_set1, bool is_set2,
                                 const char *assoc_string)
{
  sidl::array<iRel_EntityHandle> other_entities;
  sidl::array<int> offset;
  int offset_size, other_entities_size;
  try {
    lasso.getArrAssociation(iface1, iface2, ents1, ents1_size,
                            is_set1, is_set2,
                            other_entities, other_entities_size,
                            offset, offset_size);
  } catch (iBase::Error err) {
    cerr << "Failed to check " << assoc_string
         << " vertex associations." << std::endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }
  if (other_entities_size != ents1_size) {
    cerr << assoc_string << " association check returned wrong # entities."
         << endl;
    return false;
  }

  if (debug) {
    cerr << assoc_string << " associations:" << endl;
    print_entity_info(iface2, ARRAY_PTR(other_entities, iRel_EntityHandle),
                      other_entities_size, is_set2);
  }

  return true;
}

/*!
  @test
  TSTTLasso relate geom and mesh Test
  @li Check relation between geom and mesh
*/
bool relate_geom_mesh_test(iRel::Associate &lasso,
                           iGeom::Geometry &geom, iMesh::Mesh &mesh)
{
  CAST_iGeom_INTERFACE(geom, geom_topo, Topology, false);
  CAST_iBase_INTERFACE(mesh, mesh_eset, EntSet, false);
  CAST_iBase_INTERFACE(mesh, mesh_tag, Tag, false);
  CAST_iBase_INTERFACE(mesh, mesh_stag, SetTag, false);
  CAST_iBase_INTERFACE(mesh, mesh_atag, ArrTag, false);
  CAST_iRel_INTERFACE(lasso, assoc_infer, Infer, false);


    // create an association pair between geometry and mesh interfaces
  try {
    lasso.createAssociation(geom, 0, mesh, 2);
  } catch (iBase::Error err) {
    cerr << "Failed to relate create association in relate_geom_mesh_test." << std::endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

    // TEST 1: infer associations for geom vertices
    // get geom vertices
  sidl::array<iRel_EntityHandle> gentities;
  int gentities_size = 0;
  try {
    geom_topo.getEntities
      (NULL, iGeom::EntityType_VERTEX, gentities, gentities_size);
  } catch (iBase::Error err) {
    std::cerr << "Failed to get gentities by type in relate_geom_mesh_test."
              << std::endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

    // infer associations for geom vertices
  try {
    assoc_infer.inferArrAssociations(geom, mesh, gentities, gentities_size, 0);
  } catch (iBase::Error err) {
    cerr << "Failed to relate geom entities in relate_geom_mesh_test." << std::endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

    // check inferred associations
  bool success = check_inferred_associations(lasso, geom, mesh,
                                             gentities, gentities_size,
                                             false, true, "vertex");
  if (!success) return success;

    // TEST 2: infer associations for 1-d mesh sets
    // get 1-dimensional mesh entitysets
  sidl::array<iRel_EntityHandle> mentity_handles;
  int mentity_handles_size;
  try {
    mesh_eset.getEntSets(NULL, 1, mentity_handles, mentity_handles_size);
  } catch (iBase::Error err) {
    cerr << "Problem to get all entity sets." << endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

    // get geom dimension tags for mesh entitysets
  std::string dim_tag_name(GEOM_DIMENSION_TAG_NAME);
  iBase_TagHandle dim_tag_mesh;
  try {
    mesh_tag.createTag(dim_tag_name, sizeof(int), iBase::TagValueType_INTEGER, dim_tag_mesh);
  }
  catch (iBase::Error err) {
    if (err.getErrorType() != iBase::ErrorType_TAG_ALREADY_EXISTS) {
      std::cerr << "Couldn't create geom dim tag for mesh entities." << std::endl;
      cerr << "Description: " << err.getDescription() << endl;
      return false;
    }
  }

    // get 1-dimensional mesh entitysets
  std::vector<iRel_EntityHandle> mentities_vector;
  int n_mentities = 0;
  for (int i = 0; i < mentity_handles_size; i++) { // test
    sidl::array<int> dim_val_mesh = dim_val_mesh.create1d(1);
    void *mentity = mentity_handles.get(i);
    int dim = -1;
    try {
      dim = mesh_stag.getEntSetIntData(mentity, dim_tag_mesh);
    }
    catch (iBase::Error err) {
      continue;
    }

    if (dim == 1) {
      mentities_vector.push_back(mentity);
      n_mentities++;
    }
  }
  if (debug) {
    cerr << "Inferring associations for 1d entity sets:" << endl;
    print_mesh_info(mesh, &mentities_vector[0], mentities_vector.size(),
                    true);
  }

    // copy to sidl array, then infer associations for them
  sidl::array<iRel_EntityHandle> mentities = mentities.create1d(n_mentities);
  int mentities_size = n_mentities;
  for (int i = 0; i < n_mentities; i++)
    mentities.set(i, mentities_vector[i]);

  CAST_iRel_INTERFACE(lasso, lasso_infer, Infer, false);
  try {
    lasso_infer.inferArrAssociations(mesh, geom, mentities, mentities_size, 1);
  } catch (iBase::Error err) {
    cerr << "Failed to relate mesh entities in relate_geom_mesh_test."
         << endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

    // check inferred associations
  success = check_inferred_associations(lasso, mesh, geom,
                                        mentities, mentities_size,
                                        true, false,
                                        "1d mentity set");
  if (!success) return success;

    // relate all geometry and mesh entities
  try {
    lasso_infer.inferAllAssociations(geom, mesh);
  } catch (iBase::Error err) {
    cerr << "Failed to relate all geom and mesh entities in relate_geom_mesh_test."
         << endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

    // check by getting all gentities, getting associated mentities,
    // then getting gentities related to those; two gentity lists
    // should be the same
  sidl::array<iRel_EntityHandle> gentities2, mentities2, gentities3;
  int gentities2_size, mentities2_size, gentities3_size;
  sidl::array<int> offsets, offsets2;
  int offsets_size, offsets2_size;
  try {
    geom_topo.getEntities
      (NULL, iGeom::EntityType_ALL_TYPES, gentities2, gentities2_size);
    lasso.getArrAssociation(geom, mesh, gentities2, gentities2_size, 0, 1,
                            mentities2, mentities2_size,
                            offsets, offsets_size);
    lasso.getArrAssociation(mesh, geom, mentities2, mentities2_size, 1, 0,
                            gentities3, gentities3_size,
                            offsets2, offsets2_size);
  } catch (iBase::Error err) {
    std::cerr << "Failed to reverse-check all associations in relate_geom_mesh_test."
              << std::endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

  if (gentities3_size != gentities2_size) {
    cerr << "Number of reverse-checked geom entities doesn't agree with orig number" << endl;
    return false;
  }

    // check the lists
  int num_disagree = 0;
  for (int i = 0; i < gentities3_size; i++)
    if (gentities2[i] != gentities3[i]) num_disagree++;

  if (num_disagree) {
    cerr << num_disagree << " entities don't agree between 2 geom entity lists." << endl;
    return false;
  }

  return true;
}

/*!
  @test
  TSTTLasso move to test
  @li Move meshes onto the given geometry
*/
bool query_relations_test(iRel::Associate &lasso,
                          iGeom::Geometry &geom, iMesh::Mesh &mesh)
{
  CAST_iBase_INTERFACE(geom, geom_tag, Tag, false);
  CAST_iGeom_INTERFACE(geom, geom_topo, Topology, false);
  CAST_iBase_INTERFACE(mesh, mesh_eset, EntSet, false);
  CAST_iBase_INTERFACE(mesh, mesh_tag, Tag, false);

    // get all the geom entities, and find relation to some mesh entity
  sidl::array<iRel_EntityHandle> gentities;
  int gentities_size;
  try {
    geom_topo.getEntities(NULL, iGeom::EntityType_ALL_TYPES,
                          gentities, gentities_size);
  } catch (iBase::Error err) {
    cerr << "Problem to get all geom entities." << endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

    // get corresponding mesh entities
  sidl::array<iRel_EntityHandle> out_mentities =
    out_mentities.create1d(gentities_size);
  int out_mentities_size;
  sidl::array<int> offsets;
  int offsets_size = 0;

  try {
    lasso.getArrAssociation(geom, mesh,
                            gentities, gentities_size, 0, 1,
                            out_mentities, out_mentities_size,
                            offsets, offsets_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get mesh entities related to geom entities in query_relations_test."
         << endl;
    cerr << "Description: " << err.getDescription();

    for (int i = 0; i < gentities_size; i++) {
      if (out_mentities.get(i) == 0) {
        print_geom_info(geom, ARRAY_PTR(gentities, iRel_EntityHandle)+i, 1);
      }
    }

    return false;
  }

    // check the numbers
  if (out_mentities_size != gentities_size) {
    cerr << "Number of geom & related mesh entities don't match." << endl;
    return false;
  }

    // check to make sure they're mesh sets; how to do that?
  bool is_list;
  for (int i = 0; i < out_mentities_size; i++) {
    if (0 == out_mentities[i]) {
      std::cerr << "NULL entity set returned for a geom entity."
                << std::endl;
      return false;
    }

    try {
      is_list = mesh_eset.isList(out_mentities[i]);
    }
    catch (iBase::Error err) {
      if (iBase::ErrorType_SUCCESS != err.getErrorType()) {
        std::cerr << "Entity set returned from classification wasn't valid."
                  << std::endl;
        std::cerr << "Description: " << err.getDescription() << std::endl;
        return false;
      }
    }
  }

    // now turn around and check classification of those mesh entities
  sidl::array<iRel_EntityHandle> out_gentities;
  int out_gentities_size;

  try {
    lasso.getArrAssociation(mesh, geom,
                            out_mentities, out_mentities_size, 1, 0,
                            out_gentities, out_gentities_size,
                            offsets, offsets_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get geom entities related to mesh entities in query_relations_test."
         << endl;
    cerr << "Description: " << err.getDescription() << endl;
    return false;
  }

    // check that they're all non-null
  if (out_mentities_size != out_gentities_size) {
    cerr << "Number of mesh & related geom entities don't match." << endl;
    return false;
  }

    // ok, we're done
  return true;
}

/*!
  @test
  TSTTLasso move to test
  @li Move meshes onto the given geometry
*/
bool move_to_test(iRel::Associate &lasso,
                  iGeom::Geometry &geom, iMesh::Mesh &mesh)
{
  CAST_iGeom_INTERFACE(geom, geom_topo, Topology, false);
  CAST_iBase_INTERFACE(mesh, mesh_eset, EntSet, false);
  CAST_iBase_INTERFACE(mesh, mesh_tag, Tag, false);
  CAST_iBase_INTERFACE(mesh, mesh_stag, SetTag, false);

    // get all gvertices, gedges, gfaces
  sidl::array<iRel_EntityHandle> geom_ents;
  int geom_ents_size;
  bool success = get_gentities(geom, 0, 2, geom_ents, geom_ents_size);
  if (!success) return false;

    // project mesh down to geometry for those entities
  success = project_to_geom(geom, mesh, lasso, geom_ents, geom_ents_size);

  if (!success) return false;

  return true;
}


int main( int argc, char *argv[] )
{
    // Check command line arg
  std::string geom_filename;
  std::string mesh_filename;

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <geom_filename>" << " <mesh_filename>" << std::endl;
    return 1;
  }

  geom_filename = argv[1];
  mesh_filename = argv[2];
  bool result;
  int number_tests = 0;
  int number_tests_successful = 0;
  int number_tests_not_implemented = 0;
  int number_tests_failed = 0;

    // initialize the Geometry
  iGeom::Geometry geom = iGeom::Factory::newGeom("");

    // initialize the Mesh
  iMesh::Mesh mesh = iMesh::Factory::newMesh("");

    // initialize the Associate
  iRel::Associate lasso = iRel::Factory::newAssoc("");

    // Print out Header information
  cout << "\n\nTSTT TEST PROGRAM:\n\n";

    // load_geom_mesh test
  cout << "   load_geom_mesh: ";
  result = load_geom_mesh_test(geom_filename, mesh_filename, geom, mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

    // relate_geom_mesh test
  cout << "   relate_geom_mesh: ";
  result = relate_geom_mesh_test(lasso, geom, mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

    // query_relations test
  cout << "   query_relations: ";
  result = query_relations_test(lasso, geom, mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

    // move_to test
  cout << "   move_to: ";
  result = move_to_test(lasso, geom, mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

    // summary

  cout << "\nTSTT TEST SUMMARY: \n"
       << "   Number Tests:           " << number_tests << "\n"
       << "   Number Successful:      " << number_tests_successful << "\n"
       << "   Number Not Implemented: " << number_tests_not_implemented << "\n"
       << "   Number Failed:          " << number_tests_failed
       << "\n\n" << endl;

  return 0;
}

