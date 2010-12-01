#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#include "program_opt.hpp"

class ProgOpt{

  enum types{
    FLAG = 0,
    INT, REAL, STRING 
  };

  template <typename T> 
  static types get_type(){ return FLAG; } //specialized for other types at bottom of this file


  std::string shortname, longname;
  std::vector< std::string > args;
  enum types type;
  void* storage;
  int flags;
  ProgOpt* cancel_opt;

  const char* get_argstring() const { 
    switch( type ){
    case INT:
      return "<int>";
    case REAL:
      return "<val>";
    case FLAG:
      return "";
    default:
      return "<arg>";
    }
  }

public:
  ProgOpt( const std::string& longname_p, const std::string& shortname_p, int flags_p, types t = FLAG ):
    shortname( shortname_p ), longname( longname_p ), type(t), 
    storage(NULL), flags(flags_p), cancel_opt(NULL)
  {}
  
  friend class ProgOptions;
};

ProgOptions::ProgOptions( const std::string& helpstring ){
  main_help.push_back( helpstring );
  addOpt<void>( "help,h", "Show full help text", help_flag );
}

ProgOptions::~ProgOptions(){
  for( std::vector<help_line>::iterator i = option_help_strings.begin();
       i != option_help_strings.end(); ++ i )
  {
    if( (*i).first ){ delete (*i).first; }
  }

  for( std::vector<help_line>::iterator i = arg_help_strings.begin();
       i != arg_help_strings.end(); ++ i )
  {
    delete (*i).first;
  }
}


void ProgOptions::get_namestrings( const std::string& namestring, 
                                   std::string* longname, std::string* shortname )
{
  *shortname = "";
  *longname = namestring;

  size_t idx = namestring.find_first_of(',');
  if( idx != namestring.npos ){
    *longname = namestring.substr(0, idx);
    *shortname = namestring.substr( idx+1, namestring.npos );
  }

  
}



template < typename T >
void ProgOptions::addOpt( const std::string& namestring, const std::string& helpstring, 
			  T* value, int flags ){

  std::string shortname, longname;
  get_namestrings( namestring, &longname, &shortname );

  ProgOpt* opt = new ProgOpt( longname, shortname, flags, ProgOpt::get_type<T>() );
  if( value ) opt->storage = value;


  if(longname.length())  long_names[longname] = opt;
  if(shortname.length()) short_names[shortname] = opt;

  help_line help = std::make_pair( opt, helpstring );
  option_help_strings.push_back( help );

  if( flags & add_cancel_opt ){
    std::string flag = "no-" + (longname.length() ? longname : shortname );
    ProgOpt* cancel_opt = new ProgOpt( flag, "", ProgOpt::FLAG );
    cancel_opt->cancel_opt = opt;
    long_names[flag] = cancel_opt;
    std::string clear_helpstring = "Clear previous " + flag.substr(3,flag.npos) + " flag";
    help = std::make_pair( cancel_opt, clear_helpstring );
    option_help_strings.push_back( help );
  }
}


template < typename T >
void ProgOptions::addRequiredArg( const std::string& helpname, const std::string& helpstring, T* value ){
  
  ProgOpt::types type = ProgOpt::get_type<T>();

  ProgOpt* opt = new ProgOpt( helpname, "", 0,  type );
  if( value ) opt->storage = value;
  help_line help = std::make_pair( opt, helpstring );
  arg_help_strings.push_back( help );
  required_args[helpname] = opt;
}


void ProgOptions::addOptionHelpHeading( const std::string& s ){
  option_help_strings.push_back( std::make_pair( (ProgOpt*)NULL, s) );
}

void ProgOptions::printHelp( std::ostream& out ){
  
  /* Print introductory help text */
  for( std::vector<std::string>::iterator i = main_help.begin(); i!= main_help.end(); ++i ){
    if( (*i).length() ){
      out << *i << std::endl;
    }
  }

  printUsage( out );

  // max number of characters to pad argument/option names with
  // options with long names may exceed this, but will appear out of alignment in help text
  const int max_padding = 20;

  /* List required arguments, with help text */
  if( arg_help_strings.size() > 0 ){
    
    int max_arg_namelen = 0;
    
    for( std::vector<help_line>::iterator i = arg_help_strings.begin();
         i != arg_help_strings.end(); ++i )
      {
        max_arg_namelen = std::max( max_arg_namelen, (int)((*i).first->longname.length()) );
      }
    
    max_arg_namelen = std::min( max_arg_namelen+3, max_padding );

    out << "Required Arguments: " << std::endl;
    
    for( std::vector<help_line>::iterator i = arg_help_strings.begin();
         i != arg_help_strings.end(); ++i )
      {
        ProgOpt* option = (*i).first;
        std::string& info = (*i).second;

        std::stringstream s;
        s << "  " << option->longname;
        out << std::setw(max_arg_namelen) << std::left << s.str();
        out << ": " << info << std::endl;
        
      }
  }
    
  /* List options, with help text */
  out << "Options: " << std::endl;
  int max_option_prefix_len = 0;

  for( std::vector<help_line>::iterator i = option_help_strings.begin();
       i != option_help_strings.end(); ++ i )
  {
    ProgOpt* option = (*i).first;
    std::string& info = (*i).second;

    if( option ){

      if( max_option_prefix_len == 0 ){
        // iterate ahead in the option list to determine whitespace padding
        // stop if (*j).first is NULL, which indicates a help header message 
        for( std::vector<help_line>::iterator j = i; j!=option_help_strings.end() && (*j).first; ++j ){
          int len = get_option_usage_prefix( *((*j).first) ).length();
          max_option_prefix_len = std::max (max_option_prefix_len, len);
        }
      }
      max_option_prefix_len = std::min( max_option_prefix_len, max_padding );
      std::string option_prefix = get_option_usage_prefix( *option );

      out << std::setw(max_option_prefix_len) << std::left <<  option_prefix; 
      out << ": ";
    }
    else{ 
      // no option: this is a help header.  Reset max name length.
      max_option_prefix_len = 0;
    }
    out << info << std::endl;
  }
}

std::string ProgOptions::get_option_usage_prefix( const  ProgOpt& option ){
  bool has_shortname = option.shortname.length() > 0;
  bool has_longname  = option.longname.length() > 0;  
  std::string argstr = option.get_argstring();

  std::stringstream s;
  s << "  ";
  if( has_shortname ){
    
    s << "-" << option.shortname;
    if( has_longname ){ s << " "; }
    
  }
  if( has_longname ){
    
    if( has_shortname ) s << "[";
    s << "--" << option.longname; 
    if( has_shortname ) s << "]"; 
    
  }
  
  if( argstr.length() ) s << " " << argstr;
  return s.str();
}

void ProgOptions::printUsage( std::ostream& out ){

  out << "Usage: " << progname << " --help | [options] ";

  if( arg_help_strings.size() > 0 ){
    
    for( std::vector<help_line>::iterator i = arg_help_strings.begin();
         i != arg_help_strings.end(); ++i )
      {
        std::cout << (*i).first->longname << " ";        
      }

  }
    
  out << std::endl;

}


ProgOpt* ProgOptions::lookup( const std::map<std::string, ProgOpt* >& table, const std::string& arg ){
  std::map<std::string, ProgOpt*>::const_iterator it = table.find( arg );
  if ( it == table.end() ) return NULL;
  else return (*it).second;
}

ProgOpt* ProgOptions::lookup_option( const std::string& namestring ){
  std::string longname, shortname;
  get_namestrings( namestring, &longname, &shortname );
  
  ProgOpt* opt = lookup( long_names, longname );
  if( !opt ) opt = lookup( short_names, shortname );
  
  if( !opt ){
    error( "Could not look up option: " + namestring );
  }
  
  return opt;
}

void ProgOptions::error( const std::string& error ){
  std::cerr << "Error: " << error << "\n"<< std::endl;;
  printUsage( std::cerr );
  std::cerr << std::endl;
  std::exit( EXIT_FAILURE );
}

/**
 * Check the input to a given option for correctness, converting it to its expected type (e.g. int)
 * and storing the result to target, if target is non-NULL.
 * @param option Used only in error messages to state which option could not be successfully converted
 * @param arg_idx If non-NULL, evaluate the (*arg_idx)'th item in opt's args list
 */
bool ProgOptions::evaluate( const ProgOpt& opt, void* target, const std::string& option, unsigned* arg_idx ){

  unsigned idx = arg_idx ? *arg_idx : opt.args.size()-1;

  switch( opt.type ){
  case ProgOpt::FLAG:
    error("Cannot evaluate a flag");
    break;
  case ProgOpt::INT:
    {
      int temp;
      int* i = target ? reinterpret_cast<int*>(target) : &temp;
      if( opt.args.size() < 1 ){
	error( "Missing argument to " + option + " option");
      }
      const char* arg = opt.args.at(idx).c_str();
      char* p;
      *i = std::strtol( arg, &p, 10 );
      if( *p != '\0' ){ error("Bad integer argument '" + opt.args.at(idx) + "' to " + option + " option."); }
      return true;
    }
  case ProgOpt::REAL:
    {
      double temp;
      double* i = target ? reinterpret_cast<double*>(target) : &temp;
      if( opt.args.size() < 1 ){
	error( "Missing argument to " + option + " option");
      }
      const char* arg = opt.args.at(idx).c_str();
      char* p;
      *i = std::strtod( arg, &p );
      if( *p != '\0' ){ error("Bad real argument '" + opt.args.at(idx) + "' to " + option + " option."); }
      return true;
    
    }
  
  case ProgOpt::STRING:
    {
      std::string temp;
      std::string* i = target ? reinterpret_cast<std::string*>(target) : &temp;
      if( opt.args.size() < 1 ){
	error( "Missing argument to " + option + " option");
      }
      *i = opt.args.at(idx);
      return true;
    }

  }

  return false;
}

template <typename T>
bool ProgOptions::getOpt( const std::string& namestring, T* t ){
 
  ProgOpt* opt = lookup_option( namestring );

  if( ProgOpt::get_type<T>() != opt->type ){
    error( "Option '" + namestring + "' looked up with incompatible type" );
  }

  // This call to evaluate is inefficient, because opt was already evaluated when it was parsed.
  if( opt->args.size() ){
    evaluate( *opt, t, "" );
    return true;
  }
  else return false;

}

template <typename T>
void ProgOptions::getOptAllArgs( const std::string& namestring, std::vector<T>& values ){
  ProgOpt* opt = lookup_option( namestring );
  
  if( ProgOpt::get_type<T>() != opt->type ){
    error( "Option '" + namestring + "' looked up with incompatible type" );
  }
  
  values.resize( opt->args.size() );

  // These calls to evaluate are inefficient, because the arguments were evaluated when they were parsed
  for( unsigned i = 0; i < opt->args.size(); ++i ){
    evaluate( *opt, &(values[i]), "", &i );
  }

}

int ProgOptions::numOptSet( const std::string& namestring ){
  std::string longname, shortname;
  get_namestrings( namestring, &longname, &shortname );
  
  ProgOpt* opt = lookup( long_names, longname );
  if( !opt ) opt = lookup( short_names, shortname );

  if( !opt ){
    error( "Could not look up option: " + namestring );
  }
  
  return opt->args.size();

}

template <typename T>
T ProgOptions::getReqArg( const std::string& namestring ){
  
  ProgOpt* opt = lookup( required_args, namestring );
  
  if( !opt ){
    error( "Could not look up required arg: " + namestring );
  }
  
  // if parseProgramOptions succeeded, we can assume each required arg has a value,
  // so calling evaluate is valid
  T value;
  evaluate( *opt, &value, "" );
  return value; 

}


void ProgOptions::parseCommandLine( int argc, char* argv[] ){
  
  this->progname = argv[0];


  unsigned int non_opt_args_handled = 0;
  for( int i = 1; i < argc; ++i ){
    std::string arg(argv[i]);

    if( arg.length() > 1 && arg[0] == '-'){
      /* option handling */

      ProgOpt* opt = NULL;
      std::string nextarg;

      if( arg.length() > 2 && arg[1] == '-') { 
	// long-form argument parsing : split at a = symbol if one is found
	size_t idx = arg.find_first_of('=');
	opt = lookup( long_names, arg.substr( 2, idx-2 ) );
	if( idx != arg.npos ){
	  nextarg = arg.substr( idx+1, arg.npos );
	  std::cout << "equals" << std::endl;
	}

      }
      else if( arg.length() == 2 ){ 
	opt = lookup( short_names, arg.substr( 1, arg.npos ) );
      }

      
      if( !opt ){
	error ("Unknown option: " + arg );
      }

      if( opt->flags & help_flag ){
        printHelp( std::cout );
        exit( EXIT_SUCCESS );
      }

      if( opt->type != ProgOpt::FLAG ){

	// need to find a nextarg.
	if( nextarg.length() == 0 ){
	  if( (i+1) < argc ){
	    nextarg = argv[++i]; 
	  } 
	  else{
	    error( arg + " option is missing argument" );
	  }
	}
	
	opt->args.push_back( nextarg );
	evaluate( *opt, opt->storage, arg );

      }
      else{
	if( nextarg.length() != 0 ){
	  // --option=arg was specified, but option is meant to be a flag
	  error( "Unknown flag: " + arg );
	}

	// do flag operations
	if( opt->cancel_opt ){ opt->cancel_opt->args.clear(); }
        if( opt->storage ){
          *static_cast<bool*>(opt->storage) = ( opt->flags & store_false ) ? false : true;            
        }
        opt->args.push_back(""); 
      }
    } /* end option handling */
    else{ 
      /* non-optional arguments */
      
      if( non_opt_args_handled < required_args.size () ){
        ProgOpt* opt = arg_help_strings[non_opt_args_handled++].first;
        opt->args.push_back( arg );
        evaluate( *opt, opt->storage, arg );

      }
      else{ 
        error( "Unexpected argument: " + arg );
      }
    

    } /* End argument handling */
    
  }/* End loop over inputs */

  if( non_opt_args_handled < required_args.size() ){
    const std::string& missed_arg = arg_help_strings[non_opt_args_handled].first->longname; 
    error("Did not find required positional argument: " + missed_arg );
  }

}


/* Ensure g++ instantiates the template types we expect to use, 
   and also specialize the ProgOpt::get_type function for each supported type */

#define DECLARE_OPTION_TYPE(T, PO_TYPE)                                 \
  template void ProgOptions::addOpt<T>( const std::string&, const std::string&, T*, int ); \
  template bool ProgOptions::getOpt<T>( const std::string&, T* );       \
  template<> ProgOpt::types ProgOpt::get_type<T>(){ return PO_TYPE; }

#define DECLARE_VALUED_OPTION_TYPE(T, PO_TYPE)                          \
  DECLARE_OPTION_TYPE(T, PO_TYPE)                                       \
  template void ProgOptions::getOptAllArgs<T> (const std::string&, std::vector<T>& ); \
  template void ProgOptions::addRequiredArg<T>( const std::string&, const std::string&, T* ); \
  template T ProgOptions::getReqArg<T>( const std::string& );     
 
DECLARE_OPTION_TYPE(void, ProgOpt::FLAG)
DECLARE_VALUED_OPTION_TYPE(int,  ProgOpt::INT)
DECLARE_VALUED_OPTION_TYPE(double, ProgOpt::REAL)
DECLARE_VALUED_OPTION_TYPE(std::string, ProgOpt::STRING)
