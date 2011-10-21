#include "constracc.h"
#include "utils.h"
constr_acc::constr_acc(lst constr_lst,lst w_lst_in) : constraints(constr_lst),w_lst(w_lst_in)
{
}
constr_acc::constr_acc(exset constr_lst,exset w_lst_in)
{
  constraints = set2lst(constr_lst);
  w_lst = set2lst(w_lst_in);
}
bool constr_acc::add_single(const ex& constr)
{
  lst add_lst(constraints);
  add_lst.append(constr);
  if(zero_volume(add_lst,w_lst))return false;
  else
    {
      constraints.append(constr);
      return true;
    }
}
bool constr_acc::test_single(const ex& constr)
{
  lst add_lst(constraints);
  add_lst.append(constr);
  return !zero_volume(add_lst,w_lst);
}
size_t constr_acc::add_lst(lst& constr)
{
}

bool constr_acc::test_lst(lst& cl)
{
  lst add_lst(constraints);
  for(lst::const_iterator lit = cl.begin(); lit != cl.end(); ++lit)
    add_lst.append(*lit);
  return !zero_volume(add_lst,w_lst);
}
