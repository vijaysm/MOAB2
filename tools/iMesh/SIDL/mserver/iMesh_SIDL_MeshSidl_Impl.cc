// 
// File:          iMesh_SIDL_MeshSidl_Impl.cc
// Symbol:        iMesh_SIDL.MeshSidl-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// sidl Created:  20070927 14:57:59 CDT
// Generated:     20070927 14:58:07 CDT
// Description:   Server-side implementation for iMesh_SIDL.MeshSidl
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// source-line   = 5
// source-url    = file:/home/tautges/MOAB/tools/iMesh/SIDL/iMesh_SIDL.sidl
// 
#include "iMesh_SIDL_MeshSidl_Impl.hh"

// DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl._includes)
// Insert-Code-Here {iMesh_SIDL.MeshSidl._includes} (additional includes or code)
#include "iMesh.h"
extern iBase_Error iMesh_LAST_ERROR;

#define PROCESS_ERROR do { \
  if (imeshError != iBase_SUCCESS) { \
    iMesh_LAST_ERROR.error_type = iBase_ErrorType(imeshError); \
    sprintf(iMesh_LAST_ERROR.description, "Undescribed error type"); \
    this->processError(); \
  } \
} while (false)

#define PROCESS_ERROR_MSG(a,b) do { \
   iMesh_LAST_ERROR.error_type = iBase_ErrorType(a); \
   sprintf(iMesh_LAST_ERROR.description, "%s", b);\
   imeshError = a; \
   if (imeshError != iBase_SUCCESS) \
     this->processError(); \
 } while (false)

// need this for definitions in TSTTB_SNL_SIDL_defs.h
#define LOCAL_IBASE_ERROR iMesh_LAST_ERROR
#include "iBase_SIDL_defs.h"
#include <iostream>

// DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl._includes)

// user-defined constructor.
void iMesh_SIDL::MeshSidl_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl._ctor)
  iMesh_newMesh( 0, &imeshInstance, &imeshError, 0 );
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl._ctor)
}

// user-defined destructor.
void iMesh_SIDL::MeshSidl_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl._dtor)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl._dtor} (destructor)
  iMesh_dtor(imeshInstance, &imeshError);
  imeshInstance = 0;
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl._dtor)
}

// static class initializer.
void iMesh_SIDL::MeshSidl_impl::_load() {
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl._load)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl._load} (class initialization)
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  createTag[]
 */
void
iMesh_SIDL::MeshSidl_impl::createTag (
  /* in */ const ::std::string& tag_name,
  /* in */ int32_t number_of_values,
  /* in */ ::iBase::TagValueType tag_type,
  /* out */ void*& tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.createTag)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.createTag} (createTag method)
  iMesh_createTag (imeshInstance,
                   tag_name.c_str(), number_of_values, 
                   tag_type, (iBase_TagHandle*)&tag_handle, &imeshError, tag_name.length());
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.createTag)
}

/**
 * Method:  destroyTag[]
 */
void
iMesh_SIDL::MeshSidl_impl::destroyTag (
  /* in */ void* tag_handle,
  /* in */ bool forced ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.destroyTag)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.destroyTag} (destroyTag method)
  int tmp_forced = forced;
  iMesh_destroyTag (imeshInstance, (iBase_TagHandle)tag_handle, tmp_forced, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.destroyTag)
}

/**
 * Method:  getTagName[]
 */
void
iMesh_SIDL::MeshSidl_impl::getTagName (
  /* in */ void* tag_handle,
  /* out */ ::std::string& tag_name ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getTagName)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getTagName} (getTagName method)
  char tmp_name[120];
  iMesh_getTagName (imeshInstance, (iBase_TagHandle)tag_handle, tmp_name, &imeshError, 120);
  PROCESS_ERROR;
  tag_name = tmp_name;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getTagName)
}

/**
 * Method:  getTagSizeValues[]
 */
void
iMesh_SIDL::MeshSidl_impl::getTagSizeValues (
  /* in */ void* tag_handle,
  /* out */ int32_t& size_values ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getTagSizeValues)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getTagSizeValues} (getTagSizeValues method)
  iMesh_getTagSizeValues (imeshInstance, (iBase_TagHandle)tag_handle, &size_values, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getTagSizeValues)
}

/**
 * Method:  getTagSizeBytes[]
 */
void
iMesh_SIDL::MeshSidl_impl::getTagSizeBytes (
  /* in */ void* tag_handle,
  /* out */ int32_t& size_bytes ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getTagSizeBytes)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getTagSizeBytes} (getTagSizeBytes method)
  iMesh_getTagSizeBytes (imeshInstance, (iBase_TagHandle)tag_handle, &size_bytes, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getTagSizeBytes)
}

/**
 * Method:  getTagHandle[]
 */
void
iMesh_SIDL::MeshSidl_impl::getTagHandle (
  /* in */ const ::std::string& tag_name,
  /* out */ void*& tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getTagHandle)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getTagHandle} (getTagHandle method)
  iMesh_getTagHandle (imeshInstance, tag_name.c_str(),
                      (iBase_TagHandle*)&tag_handle, &imeshError, tag_name.length());
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getTagHandle)
}

/**
 * Method:  getTagType[]
 */
void
iMesh_SIDL::MeshSidl_impl::getTagType (
  /* in */ void* tag_handle,
  /* out */ ::iBase::TagValueType& tag_data_type ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getTagType)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getTagType} (getTagType method)
  iMesh_getTagType (imeshInstance, (iBase_TagHandle) tag_handle, 
                    (int*)&tag_data_type, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getTagType)
}

/**
 * Method:  getData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<char>& tag_value,
  /* out */ int32_t& tag_value_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getData} (getData method)
  CREATE_TEMP_TAG_ARRAY(tag_value);

  iMesh_getData (imeshInstance,
                 (iBase_EntityHandle)entity_handle, 
                 (iBase_TagHandle) tag_handle,
                 TEMP_ARRAY_INOUT(tag_value),
                 &imeshError);

  PROCESS_ERROR;

  ASSIGN_TAG_ARRAY(tag_value);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getData)
}

/**
 * Method:  getIntData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getIntData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* out */ int32_t& int_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getIntData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getIntData} (getIntData method)
  iMesh_getIntData (imeshInstance, (iBase_EntityHandle) entity_handle, (iBase_TagHandle) tag_handle, &int_data, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getIntData)
}

/**
 * Method:  getDblData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getDblData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* out */ double& dbl_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getDblData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getDblData} (getDblData method)
  iMesh_getDblData (imeshInstance, (iBase_EntityHandle) entity_handle, (iBase_TagHandle) tag_handle, &dbl_data, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getDblData)
}

/**
 * Method:  getEHData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEHData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* out */ void*& eh_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEHData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEHData} (getEHData method)
  iMesh_getEHData (imeshInstance, (iBase_EntityHandle) entity_handle, (iBase_TagHandle) tag_handle, 
                   (iBase_TagHandle*)&eh_data, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEHData)
}

/**
 * Method:  setData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<char> tag_value,
  /* in */ int32_t tag_value_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setData} (setData method)
  iMesh_setData (imeshInstance, (iBase_EntityHandle) entity_handle,
                 (iBase_TagHandle) tag_handle, 
                 ARRAY_PTR(tag_value, const char),
                 tag_value_size, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setData)
}

/**
 * Method:  setIntData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setIntData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* in */ int32_t tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setIntData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setIntData} (setIntData method)
  iMesh_setIntData (imeshInstance, (iBase_EntityHandle) entity_handle,
                    (iBase_TagHandle) tag_handle, tag_value, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setIntData)
}

/**
 * Method:  setDblData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setDblData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* in */ double tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setDblData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setDblData} (setDblData method)
  iMesh_setDblData (imeshInstance, (iBase_EntityHandle) entity_handle,
                    (iBase_TagHandle) tag_handle, tag_value, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setDblData)
}

/**
 * Method:  setEHData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setEHData (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle,
  /* in */ void* tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setEHData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setEHData} (setEHData method)
  iMesh_setEHData (imeshInstance, reinterpret_cast<iBase_EntityHandle>(entity_handle),
                   (iBase_TagHandle) tag_handle, 
                   reinterpret_cast<iBase_EntityHandle>(tag_value), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setEHData)
}

/**
 * Method:  getAllTags[]
 */
void
iMesh_SIDL::MeshSidl_impl::getAllTags (
  /* in */ void* entity_handle,
  /* inout */ ::sidl::array<void*>& tag_handles,
  /* out */ int32_t& tag_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getAllTags)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getAllTags} (getAllTags method)
  CREATE_TEMP_EH_ARRAY(tag_handles);

  iMesh_getAllTags (imeshInstance, (iBase_EntityHandle) entity_handle,
                    TEMP_ARRAY_INOUT(tag_handles), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, tag_handles);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getAllTags)
}

/**
 * Method:  rmvTag[]
 */
void
iMesh_SIDL::MeshSidl_impl::rmvTag (
  /* in */ void* entity_handle,
  /* in */ void* tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.rmvTag)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.rmvTag} (rmvTag method)
  iMesh_rmvTag (imeshInstance, (iBase_EntityHandle) entity_handle, (iBase_TagHandle) tag_handle, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.rmvTag)
}

/**
 * Method:  getArrData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<char>& value_array,
  /* out */ int32_t& value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getArrData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getArrData} (getArrData method)
  CREATE_TEMP_TAG_ARRAY(value_array);

  iMesh_getArrData (imeshInstance,
                    TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                    (iBase_TagHandle) tag_handle,
                    TEMP_ARRAY_INOUT(value_array), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TAG_ARRAY(value_array);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getArrData)
}

/**
 * Method:  getIntArrData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getIntArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<int32_t>& value_array,
  /* out */ int32_t& value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getIntArrData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getIntArrData} (getIntArrData method)
  CREATE_TEMP_ARRAY(int32_t, value_array);

  iMesh_getIntArrData (imeshInstance, TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles),
                       (iBase_TagHandle) tag_handle, TEMP_ARRAY_INOUT(value_array), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(int, value_array);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getIntArrData)
}

/**
 * Method:  getDblArrData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getDblArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<double>& value_array,
  /* out */ int32_t& value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getDblArrData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getDblArrData} (getDblArrData method)
  CREATE_TEMP_ARRAY(double, value_array);

  iMesh_getDblArrData (imeshInstance,
                       TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                       (iBase_TagHandle) tag_handle,
                       TEMP_ARRAY_INOUT(value_array), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(double, value_array);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getDblArrData)
}

/**
 * Method:  getEHArrData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEHArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<void*>& value_array,
  /* out */ int32_t& value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEHArrData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEHArrData} (getEHArrData method)
  CREATE_TEMP_EH_ARRAY(value_array);

  iMesh_getEHArrData (imeshInstance,
                      TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                      (iBase_TagHandle) tag_handle,
                      TEMP_ARRAY_INOUT(value_array), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, value_array);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEHArrData)
}

/**
 * Method:  setArrData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<char> value_array,
  /* in */ int32_t value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setArrData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setArrData} (setArrData method)
  iMesh_setArrData (imeshInstance,
                    TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                    (iBase_TagHandle) tag_handle,
                    TEMP_TAG_ARRAY_IN(value_array), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setArrData)
}

/**
 * Method:  setIntArrData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setIntArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<int32_t> value_array,
  /* in */ int32_t value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setIntArrData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setIntArrData} (setIntArrData method)
  iMesh_setIntArrData (imeshInstance,
                       TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                       (iBase_TagHandle) tag_handle,
                       TEMP_TYPED_ARRAY_IN(int32_t, value_array), &imeshError);
  PROCESS_ERROR;

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setIntArrData)
}

/**
 * Method:  setDblArrData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setDblArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<double> value_array,
  /* in */ int32_t value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setDblArrData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setDblArrData} (setDblArrData method)
  iMesh_setDblArrData (imeshInstance,
                       TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                       (iBase_TagHandle) tag_handle,
                       TEMP_TYPED_ARRAY_IN(double, value_array), &imeshError);
  PROCESS_ERROR;

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setDblArrData)
}

/**
 * Method:  setEHArrData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setEHArrData (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle,
  /* in */ ::sidl::array<void*> value_array,
  /* in */ int32_t value_array_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setEHArrData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setEHArrData} (setEHArrData method)
  iMesh_setEHArrData (imeshInstance,
                      TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                      (iBase_TagHandle) tag_handle,
                      TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, value_array), 
                      &imeshError);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setEHArrData)
}

/**
 * Method:  rmvArrTag[]
 */
void
iMesh_SIDL::MeshSidl_impl::rmvArrTag (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ void* tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.rmvArrTag)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.rmvArrTag} (rmvArrTag method)
  iMesh_rmvArrTag (imeshInstance,
                   TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                   (iBase_TagHandle) tag_handle, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.rmvArrTag)
}

/**
 * Method:  setEntSetData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setEntSetData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<char>& tag_value,
  /* in */ int32_t tag_value_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setEntSetData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setEntSetData} (setEntSetData method)
  iMesh_setEntSetData (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), (iBase_TagHandle) tag_handle,
                       ARRAY_PTR(tag_value, const char), tag_value_size, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setEntSetData)
}

/**
 * Method:  setEntSetIntData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setEntSetIntData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* in */ int32_t tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setEntSetIntData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setEntSetIntData} (setEntSetIntData method)
  iMesh_setEntSetIntData (imeshInstance,
                          reinterpret_cast<iBase_EntitySetHandle>(entity_set), (iBase_TagHandle) tag_handle, tag_value, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setEntSetIntData)
}

/**
 * Method:  setEntSetDblData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setEntSetDblData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* in */ double tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setEntSetDblData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setEntSetDblData} (setEntSetDblData method)
  iMesh_setEntSetDblData (imeshInstance,
                          reinterpret_cast<iBase_EntitySetHandle>(entity_set), (iBase_TagHandle) tag_handle, tag_value, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setEntSetDblData)
}

/**
 * Method:  setEntSetEHData[]
 */
void
iMesh_SIDL::MeshSidl_impl::setEntSetEHData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* in */ void* tag_value ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setEntSetEHData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setEntSetEHData} (setEntSetEHData method)
  iMesh_setEntSetEHData (imeshInstance,
                         reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                         (iBase_TagHandle) tag_handle, 
                         reinterpret_cast<iBase_EntityHandle>(tag_value), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setEntSetEHData)
}

/**
 * Method:  getEntSetData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntSetData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* inout */ ::sidl::array<char>& tag_value,
  /* out */ int32_t& tag_value_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntSetData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntSetData} (getEntSetData method)
  CREATE_TEMP_TAG_ARRAY(tag_value);
  
  iMesh_getEntSetData (imeshInstance,
                       reinterpret_cast<iBase_EntitySetHandle>(entity_set), (iBase_TagHandle) tag_handle, 
                       TEMP_ARRAY_INOUT(tag_value), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TAG_ARRAY(tag_value);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntSetData)
}

/**
 * Method:  getEntSetIntData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntSetIntData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* out */ int32_t& int_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntSetIntData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntSetIntData} (getEntSetIntData method)
  iMesh_getEntSetIntData (imeshInstance,
                          reinterpret_cast<iBase_EntitySetHandle>(entity_set), (iBase_TagHandle) tag_handle, 
                          &int_data, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntSetIntData)
}

/**
 * Method:  getEntSetDblData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntSetDblData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* out */ double& dbl_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntSetDblData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntSetDblData} (getEntSetDblData method)
  iMesh_getEntSetDblData (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                          (iBase_TagHandle) tag_handle, &dbl_data, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntSetDblData)
}

/**
 * Method:  getEntSetEHData[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntSetEHData (
  /* in */ void* entity_set,
  /* in */ void* tag_handle,
  /* out */ void*& eh_data ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntSetEHData)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntSetEHData} (getEntSetEHData method)
  iMesh_getEntSetEHData (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                         (iBase_TagHandle) tag_handle, 
                         reinterpret_cast<iBase_EntityHandle*>(&eh_data), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntSetEHData)
}

/**
 * Method:  getAllEntSetTags[]
 */
void
iMesh_SIDL::MeshSidl_impl::getAllEntSetTags (
  /* in */ void* entity_set,
  /* inout */ ::sidl::array<void*>& tag_handles,
  /* out */ int32_t& tag_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getAllEntSetTags)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getAllEntSetTags} (getAllEntSetTags method)

  CREATE_TEMP_EH_ARRAY(tag_handles);

  iMesh_getAllEntSetTags (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set),
                          TEMP_ARRAY_INOUT(tag_handles), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, tag_handles);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getAllEntSetTags)
}

/**
 * Method:  rmvEntSetTag[]
 */
void
iMesh_SIDL::MeshSidl_impl::rmvEntSetTag (
  /* in */ void* entity_set,
  /* in */ void* tag_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.rmvEntSetTag)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.rmvEntSetTag} (rmvEntSetTag method)
  iMesh_rmvEntSetTag (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), (iBase_TagHandle) tag_handle, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.rmvEntSetTag)
}

/**
 * Method:  createEntSet[]
 */
void
iMesh_SIDL::MeshSidl_impl::createEntSet (
  /* in */ bool isList,
  /* out */ void*& entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.createEntSet)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.createEntSet} (createEntSet method)
  int tmp_islist = isList;
  iMesh_createEntSet (imeshInstance, tmp_islist, reinterpret_cast<iBase_EntitySetHandle*>(&entity_set), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.createEntSet)
}

/**
 * Method:  destroyEntSet[]
 */
void
iMesh_SIDL::MeshSidl_impl::destroyEntSet (
  /* in */ void* entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.destroyEntSet)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.destroyEntSet} (destroyEntSet method)
  iMesh_destroyEntSet (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.destroyEntSet)
}

/**
 * Method:  isList[]
 */
void
iMesh_SIDL::MeshSidl_impl::isList (
  /* in */ void* entity_set,
  /* out */ int32_t& is_list ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.isList)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.isList} (isList method)
  iMesh_isList (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                &is_list, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.isList)
}

/**
 * Method:  getNumEntSets[]
 */
void
iMesh_SIDL::MeshSidl_impl::getNumEntSets (
  /* in */ void* entity_set,
  /* in */ int32_t num_hops,
  /* out */ int32_t& num_sets ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getNumEntSets)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getNumEntSets} (getNumEntSets method)
  iMesh_getNumEntSets (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                       num_hops, &num_sets, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getNumEntSets)
}

/**
 * Method:  getEntSets[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntSets (
  /* in */ void* entity_set,
  /* in */ int32_t num_hops,
  /* inout */ ::sidl::array<void*>& contained_entset_handles,
  /* out */ int32_t& contained_entset_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntSets)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntSets} (getEntSets method)
  CREATE_TEMP_EH_ARRAY(contained_entset_handles);

  iMesh_getEntSets (imeshInstance,
                    reinterpret_cast<iBase_EntitySetHandle>(entity_set), num_hops,
                    TEMP_ARRAY_INOUT(contained_entset_handles), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, contained_entset_handles);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntSets)
}

/**
 * Method:  addEntToSet[]
 */
void
iMesh_SIDL::MeshSidl_impl::addEntToSet (
  /* in */ void* entity_handle,
  /* inout */ void*& entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.addEntToSet)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.addEntToSet} (addEntToSet method)
  iMesh_addEntToSet (imeshInstance, (iBase_EntityHandle) entity_handle, reinterpret_cast<iBase_EntitySetHandle*>(&entity_set), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.addEntToSet)
}

/**
 * Method:  rmvEntFromSet[]
 */
void
iMesh_SIDL::MeshSidl_impl::rmvEntFromSet (
  /* in */ void* entity_handle,
  /* inout */ void*& entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.rmvEntFromSet)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.rmvEntFromSet} (rmvEntFromSet method)
  iMesh_rmvEntFromSet (imeshInstance, (iBase_EntityHandle) entity_handle, reinterpret_cast<iBase_EntitySetHandle*>(&entity_set), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.rmvEntFromSet)
}

/**
 * Method:  addEntArrToSet[]
 */
void
iMesh_SIDL::MeshSidl_impl::addEntArrToSet (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* inout */ void*& entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.addEntArrToSet)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.addEntArrToSet} (addEntArrToSet method)
  iMesh_addEntArrToSet (imeshInstance,
                        TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                        reinterpret_cast<iBase_EntitySetHandle*>(&entity_set), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.addEntArrToSet)
}

/**
 * Method:  rmvEntArrFromSet[]
 */
void
iMesh_SIDL::MeshSidl_impl::rmvEntArrFromSet (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* inout */ void*& entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.rmvEntArrFromSet)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.rmvEntArrFromSet} (rmvEntArrFromSet method)
  iMesh_rmvEntArrFromSet (imeshInstance,
                          TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), 
                          reinterpret_cast<iBase_EntitySetHandle*>(&entity_set), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.rmvEntArrFromSet)
}

/**
 * Method:  addEntSet[]
 */
void
iMesh_SIDL::MeshSidl_impl::addEntSet (
  /* in */ void* entity_set_to_add,
  /* inout */ void*& entity_set_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.addEntSet)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.addEntSet} (addEntSet method)
  iMesh_addEntSet (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set_to_add), 
                   reinterpret_cast<iBase_EntitySetHandle*>(&entity_set_handle), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.addEntSet)
}

/**
 * Method:  rmvEntSet[]
 */
void
iMesh_SIDL::MeshSidl_impl::rmvEntSet (
  /* in */ void* entity_set_to_remove,
  /* inout */ void*& entity_set_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.rmvEntSet)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.rmvEntSet} (rmvEntSet method)
  iMesh_rmvEntSet(imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set_to_remove),
                  reinterpret_cast<iBase_EntitySetHandle*>(&entity_set_handle), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.rmvEntSet)
}

/**
 * Method:  isEntContained[]
 */
void
iMesh_SIDL::MeshSidl_impl::isEntContained (
  /* in */ void* containing_entity_set,
  /* in */ void* entity_handle,
  /* out */ int32_t& is_contained ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.isEntContained)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.isEntContained} (isEntContained method)
  iMesh_isEntContained(imeshInstance,
                       reinterpret_cast<iBase_EntitySetHandle>(containing_entity_set), (iBase_EntityHandle) entity_handle, &is_contained, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.isEntContained)
}

/**
 * Method:  isEntSetContained[]
 */
void
iMesh_SIDL::MeshSidl_impl::isEntSetContained (
  /* in */ void* containing_entity_set,
  /* in */ void* contained_entity_set,
  /* out */ int32_t& is_contained ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.isEntSetContained)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.isEntSetContained} (isEntSetContained method)
  iMesh_isEntSetContained(imeshInstance,
                          reinterpret_cast<iBase_EntitySetHandle>(containing_entity_set), 
                          reinterpret_cast<iBase_EntitySetHandle>(contained_entity_set),
                          &is_contained, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.isEntSetContained)
}

/**
 * Method:  addPrntChld[]
 */
void
iMesh_SIDL::MeshSidl_impl::addPrntChld (
  /* inout */ void*& parent_entity_set,
  /* inout */ void*& child_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.addPrntChld)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.addPrntChld} (addPrntChld method)
  iMesh_addPrntChld (imeshInstance,
                     reinterpret_cast<iBase_EntitySetHandle*>(&parent_entity_set), 
                     reinterpret_cast<iBase_EntitySetHandle*>(&child_entity_set), 
                     &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.addPrntChld)
}

/**
 * Method:  rmvPrntChld[]
 */
void
iMesh_SIDL::MeshSidl_impl::rmvPrntChld (
  /* inout */ void*& parent_entity_set,
  /* inout */ void*& child_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.rmvPrntChld)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.rmvPrntChld} (rmvPrntChld method)
  iMesh_rmvPrntChld (imeshInstance,
                       /*inout*/ reinterpret_cast<iBase_EntitySetHandle*>(&parent_entity_set),
                     /*inout*/ reinterpret_cast<iBase_EntitySetHandle*>(&child_entity_set), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.rmvPrntChld)
}

/**
 * Method:  isChildOf[]
 */
void
iMesh_SIDL::MeshSidl_impl::isChildOf (
  /* in */ void* parent_entity_set,
  /* in */ void* child_entity_set,
  /* out */ int32_t& is_child ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.isChildOf)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.isChildOf} (isChildOf method)
  iMesh_isChildOf (imeshInstance,
                   reinterpret_cast<iBase_EntitySetHandle>(parent_entity_set), 
                   reinterpret_cast<iBase_EntitySetHandle>(child_entity_set),
                   &is_child, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.isChildOf)
}

/**
 * Method:  getNumChld[]
 */
void
iMesh_SIDL::MeshSidl_impl::getNumChld (
  /* in */ void* entity_set,
  /* in */ int32_t num_hops,
  /* out */ int32_t& num_child ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getNumChld)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getNumChld} (getNumChld method)
  iMesh_getNumChld (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                    num_hops, &num_child, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getNumChld)
}

/**
 * Method:  getNumPrnt[]
 */
void
iMesh_SIDL::MeshSidl_impl::getNumPrnt (
  /* in */ void* entity_set,
  /* in */ int32_t num_hops,
  /* out */ int32_t& num_parent ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getNumPrnt)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getNumPrnt} (getNumPrnt method)
  iMesh_getNumPrnt (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                    num_hops, &num_parent, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getNumPrnt)
}

/**
 * Method:  getChldn[]
 */
void
iMesh_SIDL::MeshSidl_impl::getChldn (
  /* in */ void* from_entity_set,
  /* in */ int32_t num_hops,
  /* inout */ ::sidl::array<void*>& child_handles,
  /* out */ int32_t& child_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getChldn)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getChldn} (getChldn method)

  CREATE_TEMP_EH_ARRAY(child_handles);

  iMesh_getChldn (imeshInstance, 
                  reinterpret_cast<iBase_EntitySetHandle>(from_entity_set), 
                  num_hops,
                  TEMP_ARRAY_INOUT(child_handles), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, child_handles);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getChldn)
}

/**
 * Method:  getPrnts[]
 */
void
iMesh_SIDL::MeshSidl_impl::getPrnts (
  /* in */ void* from_entity_set,
  /* in */ int32_t num_hops,
  /* inout */ ::sidl::array<void*>& parent_handles,
  /* out */ int32_t& parent_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getPrnts)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getPrnts} (getPrnts method)

  CREATE_TEMP_EH_ARRAY(parent_handles);

  iMesh_getPrnts (imeshInstance, 
                  reinterpret_cast<iBase_EntitySetHandle>(from_entity_set), 
                  num_hops,
                  TEMP_ARRAY_INOUT(parent_handles), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, parent_handles);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getPrnts)
}

/**
 * Method:  subtract[]
 */
void
iMesh_SIDL::MeshSidl_impl::subtract (
  /* in */ void* entity_set_1,
  /* in */ void* entity_set_2,
  /* out */ void*& result_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.subtract)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.subtract} (subtract method)
  iMesh_subtract (imeshInstance, 
                  reinterpret_cast<iBase_EntitySetHandle>(entity_set_1),
                  reinterpret_cast<iBase_EntitySetHandle>(entity_set_2), 
                  reinterpret_cast<iBase_EntitySetHandle*>(&result_entity_set), 
                  &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.subtract)
}

/**
 * Method:  intersect[]
 */
void
iMesh_SIDL::MeshSidl_impl::intersect (
  /* in */ void* entity_set_1,
  /* in */ void* entity_set_2,
  /* out */ void*& result_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.intersect)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.intersect} (intersect method)
  iMesh_intersect (imeshInstance, 
                   reinterpret_cast<iBase_EntitySetHandle>(entity_set_1),
                   reinterpret_cast<iBase_EntitySetHandle>(entity_set_2), 
                   reinterpret_cast<iBase_EntitySetHandle*>(&result_entity_set), 
                   &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.intersect)
}

/**
 * Method:  unite[]
 */
void
iMesh_SIDL::MeshSidl_impl::unite (
  /* in */ void* entity_set_1,
  /* in */ void* entity_set_2,
  /* out */ void*& result_entity_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.unite)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.unite} (unite method)
  iMesh_unite (imeshInstance, 
               reinterpret_cast<iBase_EntitySetHandle>(entity_set_1),
               reinterpret_cast<iBase_EntitySetHandle>(entity_set_2), 
               reinterpret_cast<iBase_EntitySetHandle*>(&result_entity_set), 
               &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.unite)
}

/**
 * Method:  load[]
 */
void
iMesh_SIDL::MeshSidl_impl::load (
  /* in */ void* entity_set_handle,
  /* in */ const ::std::string& name,
  /* in */ const ::std::string& options ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.load)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.load} (load method)

  iMesh_load (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set_handle), 
              name.c_str(), options.c_str(),
              &imeshError, name.length(), options.length());
  PROCESS_ERROR;
 
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.load)
}

/**
 * Method:  save[]
 */
void
iMesh_SIDL::MeshSidl_impl::save (
  /* in */ void* entity_set_handle,
  /* in */ const ::std::string& name,
  /* in */ const ::std::string& options ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.save)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.save} (save method)

  iMesh_save (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set_handle), 
              name.c_str(), options.c_str(),
              &imeshError, name.length(), options.length());
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.save)
}

/**
 * Method:  getRootSet[]
 */
void
iMesh_SIDL::MeshSidl_impl::getRootSet (
  /* out */ void*& root_set ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getRootSet)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getRootSet} (getRootSet method)

  iMesh_getRootSet (imeshInstance, 
                    reinterpret_cast<iBase_EntityHandle*>(&root_set), 
                    &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getRootSet)
}

/**
 * Method:  getGeometricDim[]
 */
void
iMesh_SIDL::MeshSidl_impl::getGeometricDim (
  /* out */ int32_t& dim ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getGeometricDim)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getGeometricDim} (getGeometricDim method)

  iMesh_getGeometricDimension(imeshInstance, &dim, &imeshError);
  PROCESS_ERROR;

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getGeometricDim)
}

/**
 * Method:  getDfltStorage[]
 */
void
iMesh_SIDL::MeshSidl_impl::getDfltStorage (
  /* out */ ::iBase::StorageOrder& dflt_storage ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getDfltStorage)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getDfltStorage} (getDfltStorage method)
  iMesh_getDfltStorage (imeshInstance, (int*)&dflt_storage, &imeshError);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getDfltStorage)
}

/**
 * Method:  getAdjTable[]
 */
void
iMesh_SIDL::MeshSidl_impl::getAdjTable (
  /* inout */ ::sidl::array< ::iMesh::AdjacencyInfo>& adjacency_table,
  /* out */ int32_t& adjacency_table_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getAdjTable)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getAdjTable} (getAdjTable method)
  CREATE_TEMP_ENUM_ARRAY(int, adjacency_table);

  iMesh_getAdjTable (imeshInstance,
                     TEMP_ARRAY_INOUT(adjacency_table), &imeshError);
  PROCESS_ERROR;
  ASSIGN_ENUM_ARRAY(adjacency_table);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getAdjTable)
}

/**
 * Method:  areEHValid[]
 */
void
iMesh_SIDL::MeshSidl_impl::areEHValid (
  /* in */ int32_t reset,
  /* out */ int32_t& are_valid ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.areEHValid)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.areEHValid} (areEHValid method)
  iMesh_areEHValid(imeshInstance, reset, &are_valid, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.areEHValid)
}

/**
 * Method:  getNumOfType[]
 */
void
iMesh_SIDL::MeshSidl_impl::getNumOfType (
  /* in */ void* entity_set_handle,
  /* in */ ::iBase::EntityType entity_type,
  /* out */ int32_t& num_type ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getNumOfType)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getNumOfType} (getNumOfType method)
  iMesh_getNumOfType (imeshInstance,
                      reinterpret_cast<iBase_EntitySetHandle>(entity_set_handle), (iBase_EntityType)entity_type, 
                      &num_type, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getNumOfType)
}

/**
 * Method:  getNumOfTopo[]
 */
void
iMesh_SIDL::MeshSidl_impl::getNumOfTopo (
  /* in */ void* entity_set_handle,
  /* in */ ::iMesh::EntityTopology entity_topology,
  /* out */ int32_t& num_topo ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getNumOfTopo)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getNumOfTopo} (getNumOfTopo method)
  iMesh_getNumOfTopo (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set_handle),
                      (iMesh_EntityTopology) entity_topology, 
                      &num_topo, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getNumOfTopo)
}

/**
 * Method:  getAllVtxCoords[]
 */
void
iMesh_SIDL::MeshSidl_impl::getAllVtxCoords (
  /* in */ void* entity_set,
  /* inout */ ::sidl::array<double>& coords,
  /* out */ int32_t& coords_size,
  /* inout */ ::sidl::array<int32_t>& in_entity_set,
  /* out */ int32_t& in_entity_set_size,
  /* inout */ ::iBase::StorageOrder& storage_order ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getAllVtxCoords)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getAllVtxCoords} (getAllVtxCoords method)
  CREATE_TEMP_ARRAY(double, coords);
  CREATE_TEMP_ARRAY(int32_t, in_entity_set);
  
  iMesh_getAllVtxCoords (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set),
                         TEMP_ARRAY_INOUT(coords),
                         TEMP_ARRAY_INOUT(in_entity_set),
                         (int*)&storage_order, &imeshError);

  PROCESS_ERROR;
  ASSIGN_TYPED_ARRAY(double, coords);
  ASSIGN_TYPED_ARRAY(int32_t, in_entity_set);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getAllVtxCoords)
}

/**
 * Method:  getVtxCoordIndex[]
 */
void
iMesh_SIDL::MeshSidl_impl::getVtxCoordIndex (
  /* in */ void* entity_set,
  /* in */ ::iBase::EntityType requested_entity_type,
  /* in */ ::iMesh::EntityTopology requested_entity_topology,
  /* in */ ::iBase::EntityType entity_adjacency_type,
  /* inout */ ::sidl::array<int32_t>& offset,
  /* out */ int32_t& offset_size,
  /* inout */ ::sidl::array<int32_t>& index,
  /* out */ int32_t& index_size,
  /* inout */ ::sidl::array< ::iMesh::EntityTopology>& entity_topologies,
  /* out */ int32_t& entity_topologies_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getVtxCoordIndex)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getVtxCoordIndex} (getVtxCoordIndex method)

  CREATE_TEMP_ARRAY(int32_t, offset);
  CREATE_TEMP_ARRAY(int32_t, index);
  CREATE_TEMP_ENUM_ARRAY(int, entity_topologies);
  
  iMesh_getVtxCoordIndex (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                          (iBase_EntityType)requested_entity_type, 
                          (iMesh_EntityTopology)requested_entity_topology,
                          (iBase_EntityType)entity_adjacency_type,
                          TEMP_ARRAY_INOUT(offset), 
                          TEMP_ARRAY_INOUT(index), 
                          TEMP_ARRAY_INOUT(entity_topologies), &imeshError);

  PROCESS_ERROR;
  ASSIGN_TYPED_ARRAY(int, offset);
  ASSIGN_TYPED_ARRAY(int, index);
  ASSIGN_ENUM_ARRAY(entity_topologies);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getVtxCoordIndex)
}

/**
 * Method:  getEntities[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntities (
  /* in */ void* entity_set,
  /* in */ ::iBase::EntityType entity_type,
  /* in */ ::iMesh::EntityTopology entity_topology,
  /* inout */ ::sidl::array<void*>& entity_handles,
  /* out */ int32_t& entity_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntities)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntities} (getEntities method)
  
  CREATE_TEMP_EH_ARRAY(entity_handles);
  
  iMesh_getEntities (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                     (iBase_EntityType)entity_type, (iMesh_EntityTopology)entity_topology,
                     TEMP_ARRAY_INOUT(entity_handles), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, entity_handles);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntities)
}

/**
 * Method:  getVtxArrCoords[]
 */
void
iMesh_SIDL::MeshSidl_impl::getVtxArrCoords (
  /* in */ ::sidl::array<void*> vertex_handles,
  /* in */ int32_t vertex_handles_size,
  /* inout */ ::iBase::StorageOrder& storage_order,
  /* inout */ ::sidl::array<double>& coords,
  /* out */ int32_t& coords_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getVtxArrCoords)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getVtxArrCoords} (getVtxArrCoords method)
  CREATE_TEMP_ARRAY(double, coords);
  
  iMesh_getVtxArrCoords (imeshInstance,
                         TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, vertex_handles),
                         (int*)&storage_order,
                         TEMP_ARRAY_INOUT(coords), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(double, coords);
  
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getVtxArrCoords)
}

/**
 * Method:  getAdjEntities[]
 */
void
iMesh_SIDL::MeshSidl_impl::getAdjEntities (
  /* in */ void* entity_set,
  /* in */ ::iBase::EntityType entity_type_requestor,
  /* in */ ::iMesh::EntityTopology entity_topology_requestor,
  /* in */ ::iBase::EntityType entity_type_requested,
  /* inout */ ::sidl::array<void*>& adj_entity_handles,
  /* out */ int32_t& adj_entity_handles_size,
  /* inout */ ::sidl::array<int32_t>& offset,
  /* out */ int32_t& offset_size,
  /* inout */ ::sidl::array<int32_t>& in_entity_set,
  /* out */ int32_t& in_entity_set_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getAdjEntities)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getAdjEntities} (getAdjEntities method)
  
  CREATE_TEMP_EH_ARRAY(adj_entity_handles);
  CREATE_TEMP_ARRAY(int32_t, offset);
  CREATE_TEMP_ARRAY(int32_t, in_entity_set);
  
  iMesh_getAdjEntities (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set), 
                        (iBase_EntityType)entity_type_requestor,
                        (iMesh_EntityTopology)entity_topology_requestor, 
                        (iBase_EntityType)entity_type_requested,
                        TEMP_ARRAY_INOUT(adj_entity_handles),
                        TEMP_ARRAY_INOUT(offset),
                        TEMP_ARRAY_INOUT(in_entity_set), &imeshError);
  PROCESS_ERROR;


  ASSIGN_TYPED_ARRAY(void*, adj_entity_handles);
  ASSIGN_TYPED_ARRAY(int32_t, offset);
  ASSIGN_TYPED_ARRAY(int32_t, in_entity_set);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getAdjEntities)
}

/**
 * Method:  initEntIter[]
 */
void
iMesh_SIDL::MeshSidl_impl::initEntIter (
  /* in */ void* entity_set_handle,
  /* in */ ::iBase::EntityType requested_entity_type,
  /* in */ ::iMesh::EntityTopology requested_entity_topology,
  /* out */ void*& entity_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.initEntIter)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.initEntIter} (initEntIter method)
  iMesh_initEntIter (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set_handle),
                     (iBase_EntityType)requested_entity_type, 
                     (iMesh_EntityTopology)requested_entity_topology,
                     reinterpret_cast<iMesh_EntityIterator*>(&entity_iterator), 
                     &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.initEntIter)
}

/**
 * Method:  getNextEntIter[]
 */
void
iMesh_SIDL::MeshSidl_impl::getNextEntIter (
  /* in */ void* entity_iterator,
  /* out */ void*& entity_handle,
  /* out */ int32_t& at_end ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getNextEntIter)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getNextEntIter} (getNextEntIter method)
  iMesh_getNextEntIter (imeshInstance, reinterpret_cast<iMesh_EntityIterator>(entity_iterator), 
                        reinterpret_cast<iBase_EntityHandle*>(&entity_handle),
                        &at_end, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getNextEntIter)
}

/**
 * Method:  resetEntIter[]
 */
void
iMesh_SIDL::MeshSidl_impl::resetEntIter (
  /* in */ void* entity_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.resetEntIter)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.resetEntIter} (resetEntIter method)
  iMesh_resetEntIter (imeshInstance, 
                      reinterpret_cast<iMesh_EntityIterator>(entity_iterator), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.resetEntIter)
}

/**
 * Method:  endEntIter[]
 */
void
iMesh_SIDL::MeshSidl_impl::endEntIter (
  /* in */ void* entity_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.endEntIter)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.endEntIter} (endEntIter method)
  iMesh_endEntIter (imeshInstance, 
                    reinterpret_cast<iMesh_EntityIterator>(entity_iterator), 
                    &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.endEntIter)
}

/**
 * Method:  getEntTopo[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntTopo (
  /* in */ void* entity_handle,
  /* out */ ::iMesh::EntityTopology& ent_topo ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntTopo)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntTopo} (getEntTopo method)
  iMesh_getEntTopo (imeshInstance, (iBase_EntityHandle) entity_handle, 
                    (int*)&ent_topo, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntTopo)
}

/**
 * Method:  getEntType[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntType (
  /* in */ void* entity_handle,
  /* out */ ::iBase::EntityType& ent_type ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntType)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntType} (getEntType method)
  iMesh_getEntType (imeshInstance, (iBase_EntityHandle) entity_handle, 
                    (int*)&ent_type, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntType)
}

/**
 * Method:  getVtxCoord[]
 */
void
iMesh_SIDL::MeshSidl_impl::getVtxCoord (
  /* in */ void* vertex_handle,
  /* out */ double& x,
  /* out */ double& y,
  /* out */ double& z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getVtxCoord)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getVtxCoord} (getVtxCoord method)

  iMesh_getVtxCoord (imeshInstance, 
                     reinterpret_cast<iBase_EntityHandle>(vertex_handle),
                     &x, &y, &z, &imeshError);
  PROCESS_ERROR;

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getVtxCoord)
}

/**
 * Method:  getEntAdj[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntAdj (
  /* in */ void* entity_handle,
  /* in */ ::iBase::EntityType entity_type_requested,
  /* inout */ ::sidl::array<void*>& adj_entity_handles,
  /* out */ int32_t& adj_entity_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntAdj)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntAdj} (getEntAdj method)

  CREATE_TEMP_EH_ARRAY(adj_entity_handles);

  iMesh_getEntAdj (imeshInstance, (iBase_EntityHandle) entity_handle,  
                   (iBase_EntityType)entity_type_requested,
                   TEMP_ARRAY_INOUT(adj_entity_handles), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, adj_entity_handles);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntAdj)
}

/**
 * Method:  getEnt2ndAdj[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEnt2ndAdj (
  /* in */ void* entity_handle,
  /* in */ ::iBase::EntityType order_adjacent_key,
  /* in */ ::iBase::EntityType entity_type_requested,
  /* inout */ ::sidl::array<void*>& adj_entity_handles,
  /* out */ int32_t& adj_entity_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEnt2ndAdj)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEnt2ndAdj} (getEnt2ndAdj method)
  CREATE_TEMP_EH_ARRAY(adj_entity_handles);

  iMesh_getEnt2ndAdj (imeshInstance, (iBase_EntityHandle) entity_handle,  
                      (iBase_EntityType)order_adjacent_key,
                      (iBase_EntityType)entity_type_requested,
                      TEMP_ARRAY_INOUT(adj_entity_handles), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, adj_entity_handles);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEnt2ndAdj)
}

/**
 * Method:  initEntArrIter[]
 */
void
iMesh_SIDL::MeshSidl_impl::initEntArrIter (
  /* in */ void* entity_set_handle,
  /* in */ ::iBase::EntityType requested_entity_type,
  /* in */ ::iMesh::EntityTopology requested_entity_topology,
  /* in */ int32_t requested_array_size,
  /* out */ void*& entArr_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.initEntArrIter)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.initEntArrIter} (initEntArrIter method)
  iMesh_initEntArrIter (imeshInstance, reinterpret_cast<iBase_EntitySetHandle>(entity_set_handle),
                        (iBase_EntityType)requested_entity_type, 
                        (iMesh_EntityTopology)requested_entity_topology,
                        requested_array_size, 
                        reinterpret_cast<iMesh_EntityArrIterator*>(&entArr_iterator), 
                        &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.initEntArrIter)
}

/**
 * Method:  getNextEntArrIter[]
 */
void
iMesh_SIDL::MeshSidl_impl::getNextEntArrIter (
  /* in */ void* entArr_iterator,
  /* inout */ ::sidl::array<void*>& entity_handles,
  /* out */ int32_t& entity_handles_size,
  /* out */ int32_t& at_end ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getNextEntArrIter)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getNextEntArrIter} (getNextEntArrIter method)
  CREATE_TEMP_EH_ARRAY(entity_handles);
  
  iMesh_getNextEntArrIter (imeshInstance, 
                           reinterpret_cast<iMesh_EntityArrIterator>(entArr_iterator),
                           TEMP_ARRAY_INOUT(entity_handles), 
                           &at_end, &imeshError);


  ASSIGN_TYPED_ARRAY(void*, entity_handles);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getNextEntArrIter)
}

/**
 * Method:  resetEntArrIter[]
 */
void
iMesh_SIDL::MeshSidl_impl::resetEntArrIter (
  /* in */ void* entArr_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.resetEntArrIter)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.resetEntArrIter} (resetEntArrIter method)
  iMesh_resetEntArrIter (imeshInstance, 
                         reinterpret_cast<iMesh_EntityArrIterator>(entArr_iterator), 
                         &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.resetEntArrIter)
}

/**
 * Method:  endEntArrIter[]
 */
void
iMesh_SIDL::MeshSidl_impl::endEntArrIter (
  /* in */ void* entArr_iterator ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.endEntArrIter)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.endEntArrIter} (endEntArrIter method)
  iMesh_endEntArrIter (imeshInstance, 
                       reinterpret_cast<iMesh_EntityArrIterator>(entArr_iterator), 
                       &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.endEntArrIter)
}

/**
 * Method:  getEntArrTopo[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntArrTopo (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* inout */ ::sidl::array< ::iMesh::EntityTopology>& topology,
  /* out */ int32_t& topology_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntArrTopo)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntArrTopo} (getEntArrTopo method)

  CREATE_TEMP_ENUM_ARRAY(int, topology);

  iMesh_getEntArrTopo (imeshInstance,
                       TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles),
                       TEMP_ARRAY_INOUT(topology), &imeshError);
  PROCESS_ERROR;


  ASSIGN_ENUM_ARRAY(topology);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntArrTopo)
}

/**
 * Method:  getEntArrType[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntArrType (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* inout */ ::sidl::array< ::iBase::EntityType>& type,
  /* out */ int32_t& type_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntArrType)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntArrType} (getEntArrType method)


  CREATE_TEMP_ENUM_ARRAY(int, type);

  iMesh_getEntArrType (imeshInstance,
                       TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles),
                       TEMP_ARRAY_INOUT(type), &imeshError);
  PROCESS_ERROR;

  ASSIGN_ENUM_ARRAY(type);
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntArrType)
}

/**
 * Method:  getEntArrAdj[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntArrAdj (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ ::iBase::EntityType entity_type_requested,
  /* inout */ ::sidl::array<void*>& adj_entity_handles,
  /* out */ int32_t& adj_entity_handles_size,
  /* inout */ ::sidl::array<int32_t>& offset,
  /* out */ int32_t& offset_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntArrAdj)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntArrAdj} (getEntArrAdj method)

  CREATE_TEMP_EH_ARRAY(adj_entity_handles);
  CREATE_TEMP_ARRAY(int32_t, offset);

  iMesh_getEntArrAdj (imeshInstance,
                      TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles),
                      (iBase_EntityType)entity_type_requested,
                      TEMP_ARRAY_INOUT(adj_entity_handles),
                      TEMP_ARRAY_INOUT(offset), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, adj_entity_handles);
  ASSIGN_TYPED_ARRAY(int32_t, offset);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntArrAdj)
}

/**
 * Method:  getEntArr2ndAdj[]
 */
void
iMesh_SIDL::MeshSidl_impl::getEntArr2ndAdj (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size,
  /* in */ ::iBase::EntityType order_adjacent_key,
  /* in */ ::iBase::EntityType entity_type_requested,
  /* inout */ ::sidl::array<void*>& adj_entity_handles,
  /* out */ int32_t& adj_entity_handles_size,
  /* inout */ ::sidl::array<int32_t>& offset,
  /* out */ int32_t& offset_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.getEntArr2ndAdj)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.getEntArr2ndAdj} (getEntArr2ndAdj method)
  CREATE_TEMP_EH_ARRAY(adj_entity_handles);
  CREATE_TEMP_ARRAY(int32_t, offset);

  iMesh_getEntArr2ndAdj (imeshInstance,
                         TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles),
                         (iBase_EntityType)order_adjacent_key,
                         (iBase_EntityType)entity_type_requested,
                         TEMP_ARRAY_INOUT(adj_entity_handles),
                         TEMP_ARRAY_INOUT(offset), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, adj_entity_handles);
  ASSIGN_TYPED_ARRAY(int32_t, offset);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.getEntArr2ndAdj)
}

/**
 * Method:  setVtxCoords[]
 */
void
iMesh_SIDL::MeshSidl_impl::setVtxCoords (
  /* in */ void* vertex_handle,
  /* in */ double x,
  /* in */ double y,
  /* in */ double z ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setVtxCoords)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setVtxCoords} (setVtxCoords method)
  iMesh_setVtxCoords (imeshInstance, 
                      reinterpret_cast<iBase_EntityHandle>(vertex_handle),
                      x, y, z, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setVtxCoords)
}

/**
 * Method:  createVtx[]
 */
void
iMesh_SIDL::MeshSidl_impl::createVtx (
  /* in */ double x,
  /* in */ double y,
  /* in */ double z,
  /* out */ void*& new_vertex_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.createVtx)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.createVtx} (createVtx method)
  iMesh_createVtx (imeshInstance,
                   x, y, z, 
                   reinterpret_cast<iBase_EntityHandle*>(&new_vertex_handle), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.createVtx)
}

/**
 * Method:  createEnt[]
 */
void
iMesh_SIDL::MeshSidl_impl::createEnt (
  /* in */ ::iMesh::EntityTopology new_entity_topology,
  /* in */ ::sidl::array<void*> lower_order_entity_handles,
  /* in */ int32_t lower_order_entity_handles_size,
  /* out */ void*& new_entity_handle,
  /* out */ ::iBase::CreationStatus& status ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.createEnt)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.createEnt} (createEnt method)
  iMesh_createEnt (imeshInstance, 
                   (iMesh_EntityTopology)new_entity_topology,
                   TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, lower_order_entity_handles),
                   reinterpret_cast<iBase_EntityHandle*>(&new_entity_handle), 
                   (int*)&status, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.createEnt)
}

/**
 * Method:  deleteEnt[]
 */
void
iMesh_SIDL::MeshSidl_impl::deleteEnt (
  /* in */ void* entity_handle ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.deleteEnt)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.deleteEnt} (deleteEnt method)
  iMesh_deleteEnt (imeshInstance, (iBase_EntityHandle) entity_handle, &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.deleteEnt)
}

/**
 * Method:  setVtxArrCoords[]
 */
void
iMesh_SIDL::MeshSidl_impl::setVtxArrCoords (
  /* in */ ::sidl::array<void*> vertex_handles,
  /* in */ int32_t vertex_handles_size,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> new_coords,
  /* in */ int32_t new_coords_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.setVtxArrCoords)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.setVtxArrCoords} (setVtxArrCoords method)
  iMesh_setVtxArrCoords (imeshInstance,
                         TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, vertex_handles), 
                         (iBase_StorageOrder)storage_order,
                         TEMP_TYPED_ARRAY_IN(double, new_coords), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.setVtxArrCoords)
}

/**
 * Method:  createVtxArr[]
 */
void
iMesh_SIDL::MeshSidl_impl::createVtxArr (
  /* in */ int32_t num_verts,
  /* in */ ::iBase::StorageOrder storage_order,
  /* in */ ::sidl::array<double> new_coords,
  /* in */ int32_t new_coords_size,
  /* inout */ ::sidl::array<void*>& new_vertex_handles,
  /* out */ int32_t& new_vertex_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.createVtxArr)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.createVtxArr} (createVtxArr method)
  CREATE_TEMP_EH_ARRAY(new_vertex_handles);

  iMesh_createVtxArr (imeshInstance,
                      num_verts, 
                      (iBase_StorageOrder)storage_order,
                      TEMP_TYPED_ARRAY_IN(double, new_coords),
                      TEMP_ARRAY_INOUT(new_vertex_handles), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, new_vertex_handles);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.createVtxArr)
}

/**
 * Method:  createEntArr[]
 */
void
iMesh_SIDL::MeshSidl_impl::createEntArr (
  /* in */ ::iMesh::EntityTopology new_entity_topology,
  /* in */ ::sidl::array<void*> lower_order_entity_handles,
  /* in */ int32_t lower_order_entity_handles_size,
  /* inout */ ::sidl::array<void*>& new_entity_handles,
  /* out */ int32_t& new_entity_handles_size,
  /* inout */ ::sidl::array< ::iBase::CreationStatus>& status,
  /* out */ int32_t& status_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.createEntArr)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.createEntArr} (createEntArr method)
  CREATE_TEMP_EH_ARRAY(new_entity_handles);
  CREATE_TEMP_ENUM_ARRAY(int, status);

  iMesh_createEntArr (imeshInstance,
                      (iMesh_EntityTopology)new_entity_topology,
                      TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, lower_order_entity_handles),
                      TEMP_ARRAY_INOUT(new_entity_handles),
                      TEMP_ARRAY_INOUT(status), &imeshError);
  PROCESS_ERROR;

  ASSIGN_TYPED_ARRAY(void*, new_entity_handles);
  ASSIGN_ENUM_ARRAY(status);

  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.createEntArr)
}

/**
 * Method:  deleteEntArr[]
 */
void
iMesh_SIDL::MeshSidl_impl::deleteEntArr (
  /* in */ ::sidl::array<void*> entity_handles,
  /* in */ int32_t entity_handles_size ) 
throw ( 
  ::iBase::Error
){
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl.deleteEntArr)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl.deleteEntArr} (deleteEntArr method)
  iMesh_deleteEntArr (imeshInstance,
                      TEMP_TYPED_ARRAY_IN(iBase_EntityHandle, entity_handles), &imeshError);
  PROCESS_ERROR;
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl.deleteEntArr)
}


// DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl._misc)
// Insert-Code-Here {iMesh_SIDL.MeshSidl._misc} (miscellaneous code)
// call this function when you want to throw an error
void
iMesh_SIDL::MeshSidl_impl::processError() throw(::iBase::Error)
{
  static void *behavior_tag = NULL;
  iBase_ErrorActions action = iBase_THROW_ERROR;

    // save this info before calling tag get function to get behavior
  std::string this_desc(iMesh_LAST_ERROR.description);
  iBase::ErrorType this_err = (iBase::ErrorType)iMesh_LAST_ERROR.error_type;

  if (0 == behavior_tag) {
    iMesh_getTagHandle(imeshInstance, "Error_Behavior", 
                       reinterpret_cast<iBase_TagHandle*>(&behavior_tag),
                       &imeshError, 16);
      // reset the iMesh error indicator to what was passed in here
    iMesh_LAST_ERROR.error_type = (iBase_ErrorType) this_err;
    iMesh_LAST_ERROR.description[0] = '\0';
  }

  if (0 != behavior_tag) {
    int tmp_val;
    iMesh_getIntData (imeshInstance, 0, 
                      reinterpret_cast<iBase_TagHandle>(behavior_tag), 
                      &tmp_val, &imeshError);
    if (iBase_SUCCESS == iMesh_LAST_ERROR.error_type) 
      action = (iBase_ErrorActions) tmp_val;
  }

  iBase::Error this_error = iBase::Error::_create();
  switch (action) {
    case iBase_THROW_ERROR:
      this_error.set(this_err, this_desc);
      throw this_error;
      break;
    case iBase_WARN_ONLY:
      std::cerr << iMesh_LAST_ERROR.description << std::endl;
      break;
    case iBase_SILENT:
      break;
  }
}

// DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl._misc)

