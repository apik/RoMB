#ifndef __ROMB_H__
#define __ROMB_H__
#include "mbintegral.h"
#include "constracc.h"
struct int_map_comparator : public std::binary_function<MBintegral::w_lst_type,MBintegral::w_lst_type, bool>
{
  bool operator()(const MBintegral::w_lst_type& a ,const MBintegral::w_lst_type& b) const
  {
    ex_is_less eilcom;
    if(a.size() > b.size())
      return false;
    else if(a.size() < b.size())
      return true;
    else
      return lexicographical_compare(a.begin(),a.end(),b.begin(),b.end(),eilcom);
  }
};




class RoMB_loop_by_loop
{
  lst nu_sym;
  MBlst int_lst;
  std::map<MBintegral::w_lst_type,ex,int_map_comparator> intmap;
  exmap w_shared;
  ConstrAcc constraints_;
public:
  /**
   *
   *  loop momentums,propagator expressions,
   *  invariants substitutions,propagator powers,number of loops
   *
   \param k_lst loop momentums list
   \param p_lst propagator expressions list
   \param subs_lst invariants substitutions list
   \param nu propagator powers list
   \param l number of loops 
   
   \return 
   *
  */
  RoMB_loop_by_loop(lst,lst,lst,lst,bool subs_U = true);
  void merge();
  void integrate(lst,int exp_order = 1);
  void integrate_map(lst,int exp_order = 1);

  std::pair<ex,ex> expand_and_integrate_map(ex ,MBintegral::w_lst_type,exmap, lst, int ); // up to O(eps^1) 

  void print_mathematica(MBintegral);

  MBtree MBcontinue_tree(MBintegral rootint,ex eps0 = 0);

  MBlst MBcontinue(MBintegral rootint,ex eps0 = 0);

// utils
  NearestPoleParams GetLeadingEps(MBintegral, numeric, numeric eps0 = 0);
};


#endif // __ROMB_H__
