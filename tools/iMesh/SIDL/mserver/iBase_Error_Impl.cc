// 
// File:          iBase_Error_Impl.cc
// Symbol:        iBase.Error-v0.7
// Symbol Type:   class
// Babel Version: 0.10.12
// sidl Created:  20070628 20:55:23 CDT
// Generated:     20070628 20:55:31 CDT
// Description:   Server-side implementation for iBase.Error
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// source-line   = 84
// source-url    = file:/home/tautges/MOAB/tools/iMesh/SIDL/mserver/../iBase.sidl
// 
#include "iBase_Error_Impl.hh"

// DO-NOT-DELETE splicer.begin(iBase.Error._includes)
// Insert-Code-Here {iBase.Error._includes} (additional includes or code)
// DO-NOT-DELETE splicer.end(iBase.Error._includes)

// user-defined constructor.
void iBase::Error_impl::_ctor() {
  // DO-NOT-DELETE splicer.begin(iBase.Error._ctor)
  // Insert-Code-Here {iBase.Error._ctor} (constructor)
  // DO-NOT-DELETE splicer.end(iBase.Error._ctor)
}

// user-defined destructor.
void iBase::Error_impl::_dtor() {
  // DO-NOT-DELETE splicer.begin(iBase.Error._dtor)
  // Insert-Code-Here {iBase.Error._dtor} (destructor)
  // DO-NOT-DELETE splicer.end(iBase.Error._dtor)
}

// static class initializer.
void iBase::Error_impl::_load() {
  // DO-NOT-DELETE splicer.begin(iBase.Error._load)
  // Insert-Code-Here {iBase.Error._load} (class initialization)
  // DO-NOT-DELETE splicer.end(iBase.Error._load)
}

// user-defined static methods: (none)

// user-defined non-static methods:
/**
 * Method:  set[]
 */
void
iBase::Error_impl::set (
  /* in */ ::iBase::ErrorType error,
  /* in */ const ::std::string& description ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(iBase.Error.set)
  // Insert-Code-Here {iBase.Error.set} (set method)
  errorType = error;
  errorDescription = description;
  // DO-NOT-DELETE splicer.end(iBase.Error.set)
}

/**
 * Method:  getErrorType[]
 */
void
iBase::Error_impl::getErrorType (
  /* out */ ::iBase::ErrorType& err_type ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(iBase.Error.getErrorType)
  // Insert-Code-Here {iBase.Error.getErrorType} (getErrorType method)
  err_type = errorType;
  // DO-NOT-DELETE splicer.end(iBase.Error.getErrorType)
}

/**
 * Method:  get[]
 */
void
iBase::Error_impl::get (
  /* out */ ::iBase::ErrorType& err,
  /* out */ ::std::string& description ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(iBase.Error.get)
  // Insert-Code-Here {iBase.Error.get} (get method)
  err = errorType;
  description = errorDescription;
  // DO-NOT-DELETE splicer.end(iBase.Error.get)
}

/**
 * Method:  getDescription[]
 */
void
iBase::Error_impl::getDescription (
  /* out */ ::std::string& description ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(iBase.Error.getDescription)
  // Insert-Code-Here {iBase.Error.getDescription} (getDescription method)
  description = errorDescription;
  // DO-NOT-DELETE splicer.end(iBase.Error.getDescription)
}

/**
 * Method:  echo[]
 */
void
iBase::Error_impl::echo (
  /* in */ const ::std::string& label ) 
throw () 
{
  // DO-NOT-DELETE splicer.begin(iBase.Error.echo)
  // Insert-Code-Here {iBase.Error.echo} (echo method)
  errorDescription = label + errorDescription;
  // DO-NOT-DELETE splicer.end(iBase.Error.echo)
}


// DO-NOT-DELETE splicer.begin(iBase.Error._misc)
// Insert-Code-Here {iBase.Error._misc} (miscellaneous code)
// DO-NOT-DELETE splicer.end(iBase.Error._misc)

