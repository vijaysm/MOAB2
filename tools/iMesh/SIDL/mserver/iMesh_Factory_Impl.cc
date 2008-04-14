// 
// File:          iMesh_Factory_Impl.cc
// Symbol:        iMesh.Factory-v0.7
// Symbol Type:   class
// Babel Version: 0.10.12
// sidl Created:  20080414 17:40:04 GMT
// Generated:     20080414 17:40:09 GMT
// Description:   Server-side implementation for iMesh.Factory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// source-line   = 296
// source-url    = file:/home/tautges/MOAB/tools/iMesh/SIDL/iMesh.sidl
// 
#include "iMesh_Factory_Impl.hh"

// DO-NOT-DELETE splicer.begin(iMesh.Factory._includes)
// Insert-Code-Here {iMesh.Factory._includes} (additional includes or code)
#include "iMesh_SIDL.hh"
// DO-NOT-DELETE splicer.end(iMesh.Factory._includes)

// user-defined constructor.
void iMesh::Factory_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(iMesh.Factory._ctor)
  // Insert-Code-Here {iMesh.Factory._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(iMesh.Factory._ctor)
}

// user-defined destructor.
void iMesh::Factory_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(iMesh.Factory._dtor)
  // Insert-Code-Here {iMesh.Factory._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(iMesh.Factory._dtor)
}

// static class initializer.
void iMesh::Factory_impl::_load() {
  // DO-NOT-DELETE splicer.begin(iMesh.Factory._load)
  // Insert-Code-Here {iMesh.Factory._load} (class initialization)
  // DO-NOT-DELETE splicer.end(iMesh.Factory._load)
}

// user-defined static methods:
/**
 * Method:  newMesh[]
 */
void
iMesh::Factory_impl::newMesh (
  /* in */ const ::std::string& option,
  /* out */ ::iMesh::Mesh& new_mesh ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh.Factory.newMesh)
  // Insert-Code-Here {iMesh.Factory.newMesh} (newMesh method)
  try {
    new_mesh = iMesh_SIDL::MeshSidl::_create();
  } catch (iBase::Error err) {
    throw err;
  }
  // DO-NOT-DELETE splicer.end(iMesh.Factory.newMesh)
}


// user-defined non-static methods: (none)

// DO-NOT-DELETE splicer.begin(iMesh.Factory._misc)
// Insert-Code-Here {iMesh.Factory._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(iMesh.Factory._misc)

