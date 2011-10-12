#ifndef __ROMB_H__
#define __ROMB_H__
#include "mbintegral.h"

class RoMB_loop_by_loop
{
  MBlst int_lst;
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
  RoMB_loop_by_loop(lst,lst,lst,lst);

  void integrate(lst, int exp_order = 1);
};
#endif // __ROMB_H__
