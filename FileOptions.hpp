/*
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/**\file FileOptions.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2007-08-21
 */

#ifndef FILE_OPTIONS_HPP
#define FILE_OPTIONS_HPP

#include <string>
#include <vector>
#include "MBTypes.h"

/**\brief Parse options string passed to file IO routines
 *
 * This is a utility class used by file-IO-related code to 
 * parse the options string passed to MBCore::load_file and
 * MBCore::write_file
 */
class FileOptions {
public:

  /*\param options_string The concatenation of a list of 
   *          options, separated either by the default separator
   *          (semicolon) with a custom separator specified at
   *          the beginning of the string (semicolon followed by
   *          destired separator character.)
   */
  FileOptions( const char* option_string );
  
  FileOptions( const FileOptions& copy );
  FileOptions& operator=( const FileOptions& copy );
  
  ~FileOptions();
  
  /**\brief Check for option with no value 
   *
   * Check for and remove an option w/out a value.
   *\param name The option name
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but has value
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  MBErrorCode get_null_option( const char* name ) const;
  
  /**\brief Check for option with an integer value.
   *
   * Check for and remove an option with an integer value
   *\param name The option name
   *\param value Output. The value.
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but does not have an integer value
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  MBErrorCode get_int_option( const char* name, int& value ) const;
  
  /**\brief Check for option with a double value.
   *
   * Check for and remove an option with a double value
   *\param name The option name
   *\param value Output. The value.
   *\param remove If true (default), remove the option from the list.
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but does not have a double value
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  MBErrorCode get_real_option( const char* name, double& value ) const;
  
  /**\brief Check for option with any value.
   *
   * Check for and remove an option with any value.
   *\param name The option name
   *\param value Output. The value.
   *\param remove If true (default), remove the option from the list.
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but does not have a value
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  MBErrorCode get_str_option( const char* name, std::string& value ) const;
  
  /**\brief Check for option 
   *
   * Check for and remove and option
   *\param name The option name
   *\param value The option value, or an empty string if no value.
   *\param remove If true (default), remove the option from teh list.
   *\return MB_SUCCESS or MB_ENTITY_NOT_FOUND
   */
  MBErrorCode get_option( const char* name, std::string& value ) const;
   
  /**\brief Check for option for which the value is an ID list
   *
   * Check for and remove an ID list option.  The value is expected to
   * be a comma-separated list of ID ranges, where an ID range can be 
   * either a single integer value or a range of integer values separated
   * by a dash ('-').
   *
   *\param name The option name
   *\param value Output. The value.
   *\param remove If true (default), remove the option from the list.
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but does not contain an ID list
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  //MBErrorCode get_id_list_option( const char* name, std::vector<unsigned>& value, bool remove = true );
  
  /** number of options */
  inline unsigned size() const 
    { return mOptions.size(); }
  
  /** true if no options */
  inline bool empty() const 
    { return mOptions.empty(); }
  
  /** Get list of options */
  void get_options( std::vector<std::string>& list ) const;
  
private:
  
  /**\brief Check for option 
   *
   * Check for and remove and option
   *\param name The option name
   *\param value The option value, or an empty string if no value.
   *\param remove If true (default), remove the option from teh list.
   *\return MB_SUCCESS or MB_ENTITY_NOT_FOUND
   */
  MBErrorCode get_option( const char* name, const char*& value) const;

  char* mData;
  std::vector<const char*> mOptions;

    /** Case-insensitive compare of name with option value. */
  static bool compare( const char* name, const char* option );
};

#endif

