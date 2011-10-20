#ifndef __CONSTRACC_H__
#define __CONSTRACC_H__
#include <ginac/ginac.h>
using namespace GiNaC;
class constr_acc
{
  lst constraints;
  lst w_lst;
 public:
  constr_acc(lst constr_lst,lst w_lst);
  constr_acc(exset,exset);
  bool add_single(const ex&);
  bool test_single(const ex&);
  size_t add_lst(lst&);
};
#endif // __CONSTRACC_H__
