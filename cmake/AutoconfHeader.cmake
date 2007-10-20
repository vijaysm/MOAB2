macro ( autoconf_header INFILE OUTFILE )

file( READ ${INFILE} autoconf_HEADER_IN )

string( REGEX MATCHALL "#undef ([^\n]*)" autoconf_VARS "${autoconf_HEADER_IN}" )
set( autoconf_HEADER "${autoconf_HEADER_IN}" )
foreach ( VAR ${autoconf_VARS} )
  string ( REGEX REPLACE "#undef (.*)" "\\1" VAR "${VAR}" )
  # CMake variables should usually be left undefined if they are blank or 0...
  if ( ${VAR} )
    string( REGEX REPLACE "#undef ${VAR}\n" "#define ${VAR} ${${VAR}}\n" autoconf_HEADER "${autoconf_HEADER}" )
  endif ( ${VAR} )
  # ... but always define version numbers, even if they have "0" as a value
  if ( \"${VAR}\" MATCHES ".*VERSION.*" )
    string( REGEX REPLACE "#undef ${VAR}\n" "#define ${VAR} ${${VAR}}\n" autoconf_HEADER "${autoconf_HEADER}" )
  endif ( \"${VAR}\" MATCHES ".*VERSION.*" )
endforeach ( VAR )
string( CONFIGURE "${autoconf_HEADER}"  autoconf_HEADER_OUT )

file( WRITE ${OUTFILE} "${autoconf_HEADER_OUT}" )

endmacro( autoconf_header )
