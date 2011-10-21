#ifndef __CONSTRACC_H__
#define __CONSTRACC_H__
#include <ginac/ginac.h>
using namespace GiNaC;
class constr_acc
{
  lst constraints;
  lst w_lst;
 public:
  constr_acc(lst,lst);
  constr_acc(exset,exset);
  bool add_single(const ex&);
  bool test_single(const ex&);
  size_t add_lst(lst&);
  bool test_lst(lst& cl);
};
#endif // __CONSTRACC_H__
