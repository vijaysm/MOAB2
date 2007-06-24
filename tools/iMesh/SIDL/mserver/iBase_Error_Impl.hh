// 
// File:          iBase_Error_Impl.hh
// Symbol:        iBase.Error-v0.7
// Symbol Type:   class
// Babel Version: 0.10.12
// sidl Created:  20070624 16:06:51 CDT
// Generated:     20070624 16:07:38 CDT
// Description:   Server-side implementation for iBase.Error
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// source-line   = 84
// source-url    = file:/home/tautges/MOAB/tools/iMesh/SIDL/iBase.sidl
// xml-url       = /home/tautges/MOAB/tools/iMesh/SIDL/repo/iBase.Error-v0.7.xml
// 

#ifndef included_iBase_Error_Impl_hh
#define included_iBase_Error_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_iBase_Error_IOR_h
#include "iBase_Error_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_iBase_Error_hh
#include "iBase_Error.hh"
#endif
#ifndef included_iBase_ErrorType_hh
#include "iBase_ErrorType.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(iBase.Error._includes)
// Insert-Code-Here {iBase.Error._includes} (includes or arbitrary code)
// DO-NOT-DELETE splicer.end(iBase.Error._includes)

namespace iBase { 

  /**
   * Symbol "iBase.Error" (version 0.7)
   */
  class Error_impl
  // DO-NOT-DELETE splicer.begin(iBase.Error._inherits)
  // Insert-Code-Here {iBase.Error._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(iBase.Error._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Error self;

    // DO-NOT-DELETE splicer.begin(iBase.Error._implementation)
    // Insert-Code-Here {iBase.Error._implementation} (additional details)
    ::iBase::ErrorType errorType;
    std::string errorDescription;
    // DO-NOT-DELETE splicer.end(iBase.Error._implementation)

  private:
    // private default constructor (required)
    Error_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Error_impl( struct iBase_Error__object * s ) : self(s,true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Error_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:

    /**
     * user defined non-static method.
     */
    void
    set (
      /* in */ ::iBase::ErrorType error,
      /* in */ const ::std::string& description
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    getErrorType (
      /* out */ ::iBase::ErrorType& err_type
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    get (
      /* out */ ::iBase::ErrorType& err,
      /* out */ ::std::string& description
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    getDescription (
      /* out */ ::std::string& description
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    echo (
      /* in */ const ::std::string& label
    )
    throw () 
    ;

  };  // end class Error_impl

} // end namespace iBase

// DO-NOT-DELETE splicer.begin(iBase.Error._misc)
// Insert-Code-Here {iBase.Error._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(iBase.Error._misc)

#endif
