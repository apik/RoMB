#include "constracc.h"
#include "utils.h"
/*constr_acc::constr_acc(lst constr_lst,lst w_lst_in) : constraints(constr_lst),w_lst(w_lst_in)
{
}
*/
constr_acc::constr_acc(MBintegral::w_lst_type constr_lst,MBintegral::w_lst_type w_lst_in)
{
  BOOST_FOREACH(ex ce,constr_lst)
    {
      constraints.append(ce);
    }
  BOOST_FOREACH(ex we,w_lst_in)
    {
      w_lst.append(we);
    }
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
/*size_t constr_acc::add_lst(lst& constr)
{
}
*/
bool constr_acc::test_lst(lst& cl)
{
  lst add_lst(constraints);
  for(lst::const_iterator lit = cl.begin(); lit != cl.end(); ++lit)
    add_lst.append(*lit);
  return !zero_volume(add_lst,w_lst);
}
bool interior_point(MBintegral::p_lst_type ineq_lst, exmap subs_map)
{
  for(lst::const_iterator it = ineq_lst.begin(); it != ineq_lst.end(); ++it)
    {
        cout<<"int point: " << *it << " --> " <<it->subs(subs_map)<<endl;
       if((it->subs(subs_map) <0)) return false;
    }
  return true;
}
