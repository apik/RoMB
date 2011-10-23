#ifndef __CONSTRACC_H__
#define __CONSTRACC_H__
#include <ginac/ginac.h>
#include "mbintegral.h"
using namespace GiNaC;
class constr_acc
{
  lst constraints;
  lst w_lst;
 public:
  constr_acc(MBintegral::w_lst_type,MBintegral::w_lst_type);
  //  constr_acc(MBintegral::w_lst_type,MBintegral::w_lst_type);
  bool add_single(const ex&);
  bool test_single(const ex&);
  size_t add_lst(lst&);
  bool test_lst(lst& cl);
};
#endif // __CONSTRACC_H__
