/*
 * CoreOptions.hpp
 */

#ifndef COREOPTIONS_HPP_
#define COREOPTIONS_HPP_
namespace moab
{
class CoreOptions
{
public:
  CoreOptions(double s) : option_seq(s) {};
  double get_sequence_option() {return option_seq;}

  void set_sequence_option(double factor) {  option_seq = factor;}
private:
  double option_seq;
};

extern CoreOptions coreopts;
}
#endif /* COREOPTIONS_HPP_ */
