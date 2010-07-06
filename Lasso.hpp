#ifndef LASSO_HPP
#define LASSO_HPP

#include <set>
#include <vector>
#include "iRel.h"

class AssocPair;

class Lasso 
{
public:
  Lasso(){}

  virtual ~Lasso();

    //! find a pair equivalent to these ifaces, passed as pointer to
    //! SIDL interface or interface instance
  AssocPair *find_pair(void *iface0, void *iface1,
                       bool *switched = NULL);

  AssocPair *find_pair(IfaceType type1, IfaceType type2,
                       bool *switched = NULL);
  
  void find_pairs(void *iface, std::vector<AssocPair*> &iface_pairs);

  int insert_pair(AssocPair *this_pair);
  int erase_pair(AssocPair *this_pair);
private:
  std::set<AssocPair*> assocPairs;
};

#endif // #ifndef LASSO_HPP

