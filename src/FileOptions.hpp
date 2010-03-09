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
   * Check for an option w/out a value.
   *\param name The option name
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but has value
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  MBErrorCode get_null_option( const char* name ) const;
  
  /**\brief Check for option with an integer value.
   *
   * Check for an option with an integer value
   *\param name The option name
   *\param value Output. The value.
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but does not have an integer value
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  MBErrorCode get_int_option( const char* name, int& value ) const;
  
  /**\brief Check for option with a double value.
   *
   * Check for an option with a double value
   *\param name The option name
   *\param value Output. The value.
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but does not have a double value
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  MBErrorCode get_real_option( const char* name, double& value ) const;
  
  /**\brief Check for option with any value.
   *
   * Check for an option with any value.
   *\param name The option name
   *\param value Output. The value.
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but does not have a value
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  MBErrorCode get_str_option( const char* name, std::string& value ) const;
  
  /**\brief Check for option 
   *
   * Check for an option
   *\param name The option name
   *\param value The option value, or an empty string if no value.
   *\return MB_SUCCESS or MB_ENTITY_NOT_FOUND
   */
  MBErrorCode get_option( const char* name, std::string& value ) const;
  
  /**\brief Check the string value of an option
   *
   * Check which of a list of possible values a string option contains.
   *\param name The option name
   *\param values A NULL-terminated array of C-style strings enumerating
   *              the possible option values.
   *\param index  Output: The index into <code>values</code> for the
   *              option value.
   *\return MB_SUCCESS if matched name and value.
   *        MB_ENTITY_NOT_FOUND if the option was not specified
   *        MB_FAILURE if the option value is not in the input <code>values</code> array.
   */
  MBErrorCode match_option( const char* name, const char* const* values, int& index ) const;
  
  /**\brief Check if an option matches a string value
   *
   * Check if the value for an option is the passed string.
   *\param name The option name
   *\param value The expected value.
   *\return MB_SUCCESS if matched name and value.
   *        MB_ENTITY_NOT_FOUND if the option was not specified
   *        MB_FAILURE if the option value doesn't match the passed string/
   */
  MBErrorCode match_option( const char* name, const char* value ) const;
  
  /**\brief Check for option for which the value is a list of ints
   *
   * Check for an option which is an int list.  The value is expected to
   * be a comma-separated list of int ranges, where an int range can be 
   * either a single integer value or a range of integer values separated
   * by a dash ('-').
   *
   *\param name The option name
   *\param values Output. The list of integer values.
   *\return - MB_SUCCESS if option is found
   *        - MB_TYPE_OUT_OF_RANGE if options is found, but does not contain an ID list
   *        - MB_ENTITY_NOT_FOUND if option is not found.
   */
  MBErrorCode get_ints_option( const char* name, std::vector<int>& values) const;
  
  /** number of options */
  inline unsigned size() const 
    { return mOptions.size(); }
  
  /** true if no options */
  inline bool empty() const 
    { return mOptions.empty(); }
  
  /** Get list of options */
  void get_options( std::vector<std::string>& list ) const;
  
  /** Check if all options have been looked at */
  bool all_seen() const;
  
  /** Get first unseen option */
  MBErrorCode get_unseen_option( std::string& value ) const;
  
private:
  
  /**\brief Check for option 
   *
   * Check for an option
   *\param name The option name
   *\param value The option value, or an empty string if no value.
   *\return MB_SUCCESS or MB_ENTITY_NOT_FOUND
   */
  MBErrorCode get_option( const char* name, const char*& value) const;

  char* mData;
  std::vector<const char*> mOptions;
  mutable std::vector<bool> mSeen;

    /** Case-insensitive compare of name with option value. */
  static bool compare( const char* name, const char* option );
};

#endif

