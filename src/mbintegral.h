#ifndef __MBINTEGRAL_H__
#define __MBINTEGRAL_H__
#include <iostream>
#include <iomanip>
#include <string>
#include <numeric>
#include <map>
#include <exception>
#include <ginac/ginac.h>
#include <boost/lexical_cast.hpp>
#include <boost/utility.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/container/map/map_fwd.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <boost/fusion/sequence/intrinsic/at_key.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/container/generation/make_map.hpp>
#include <boost/fusion/include/make_map.hpp>



#include "lemon/lp.h" // linear programming solver
#include <boost/random.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/special_functions/modf.hpp>



//#include <boost/icl/interval.hpp>
//#include <boost/icl/discrete_interval.hpp>
//#include <boost/icl/continuous_interval.hpp>
//#include <boost/icl/right_open_interval.hpp>
//#include <boost/icl/left_open_interval.hpp>
//#include <boost/icl/closed_interval.hpp>
//#include <boost/icl/open_interval.hpp>
//#include <boost/icl/functors.hpp>

#include <iterator>
#include <functional>

#include <gmpxx.h>


#include "romb_excompiler.h"
#include "collect_square.h"
//#include "sim.h"
#include <cuba.h>
// binder
#include <boost/mem_fn.hpp>
#include <boost/bind.hpp>

#include "utils.h"
#include "uf.h"
#include "tree.hh"
#include "constants.h"
#include "compat.h"
using namespace GiNaC;
using namespace boost;
using std::cout;
using std::setw;
using std::left;
using std::endl;
using std::string;
using std::distance;
using std::make_pair;
using std::pair;
using std::accumulate;


class MBintegral
{


  exlist w_lst;
  exlist gamma_poles;
  ex full_int_expr;
  exmap eps_w_current;
  exmap w_current;
  relational eps_current;
  int tree_level;
  ex res_pole;
  bool opt_flag;
public:
  typedef std::list<double>::iterator gamma_iterator;
  //  typedef lst::const_iterator pole_iterator;
  typedef exlist w_lst_type;
  typedef exlist p_lst_type;
  typedef w_lst_type::iterator pole_iterator;
  typedef p_lst_type::iterator w_iterator;

  ex_is_lesseq_degrevlex * comparator;
  void fix_inv();
  std::pair<pole_iterator,pole_iterator> gamma_iters()
  {
    return std::make_pair(gamma_poles.begin(),gamma_poles.end());
  }
  std::pair<w_iterator,w_iterator> w_iters()
  {
    return std::make_pair(w_lst.begin(),w_lst.end());
  }

  exlist::size_type w_size()
    {
      return w_lst.size();
    }
  exlist::size_type p_size()
    {
      return gamma_poles.size();
    }
  void add_pole(const ex& pex_in)
  {
    if(has_w(pex_in).nops() > 0)
      gamma_poles.push_back(pex_in);
    gamma_poles.unique();
  }
  void set_respole(const ex& setr)
  {
    res_pole = setr;
  }
  ex get_respole()
  {
    return res_pole;
  }
  void set_optimizable(const bool& op)
  {
    opt_flag = op;
  }
  bool get_optimizable()
  {
    return opt_flag;
  }
  lst poles_with_eps();
 
  int get_level()
  {
    return tree_level;
  }
  void inc_level()
  {
    tree_level++;
  }
  void set_level(int lev)
  {
    tree_level = lev;
  }
 
  // lst get_pole_lst();

  //void update_poles_from_ex();
  exset poles_from_ex(ex);

  //  lst get_poles();
  
  p_lst_type  get_poles()
  {
    return gamma_poles;
  }
  
  /*
  void set_poles_set(const exset& inset)
  {
    gamma_poles = inset;
  } 
  */
  
  w_lst_type  get_w_lst()
  {
    return w_lst;
  }
  
  /*
  exset get_w_set()
  {
    return w_lst;
  }
  */
  
  w_lst_type get_w_eps()
  {
    w_lst_type new_set(w_lst);
    new_set.push_back(get_symbol("eps"));
    return new_set;
  }
  

  exmap start_point_diff_w(w_lst_type,p_lst_type);
  exmap new_point();

  exmap get_point()
  {
    return eps_w_current;
  }
  exmap get_w()
  {
    return w_current;
  }
  relational get_eps()
  {
    return eps_current;
  }

  MBintegral res(relational,ex,relational);
  static lst has_w(const ex& ,w_lst_type ); 
  lst has_w(const ex&);
  ex get_expr()
  {
    return full_int_expr; 
  }
  
 ex operator*=(ex to_mul)
  {
    full_int_expr*=to_mul;
    return full_int_expr;
  }

 void operator+=(MBintegral& to_add)
  {

    
    w_lst.insert(w_lst.end(),to_add.w_iters().first,to_add.w_iters().second);
    cout<<"in add"<<endl;
    gamma_poles.insert(gamma_poles.end(),to_add.gamma_iters().first,to_add.gamma_iters().second);
  }


  // add terms to w_lst  
  /*void insert_w_lst(lst w_to_add)
  {
    for(lst::const_iterator wtit  = w_to_add.begin();wtit != w_to_add.end(); ++wtit)
      w_lst.insert(*wtit);
      }
  */
  void insert_w(ex w)
  {
    w_lst.push_back(w);
  }
  /*
  // add terms to gamma_poles_lst  
  void insert_pole_lst(lst poles_to_add)
  {
    for(lst::const_iterator wtit  = poles_to_add.begin();wtit != poles_to_add.end(); ++wtit)
      gamma_poles.insert(*wtit);
  }
  */
  // good
  
  bool pole_has(ex ex_has,ex cex)
  {
    return ex_has.has(cex);
  }

    

  void insert_pole(ex pole)
  {
   
    //inserting only poles with W, and exists in expr
    if (0 < std::count_if(w_lst.begin(),w_lst.end(), bind(&MBintegral::pole_has,ref(this),pole,_1)))
      //        && (full_int_expr.has(tgamma(pole)) || full_int_expr.has(psi(pole)) || full_int_expr.has(psi(wild(),pole)) ) )
      gamma_poles.push_back(pole);
  }
  
  void barnes1()
  {
    exmap match_map;
    if( full_int_expr.match(wild(5)*tgamma(wild(1)+wild())*tgamma(wild(2)+wild())*tgamma(wild(3)-wild())*tgamma(wild(4)-wild()),match_map))
      cout<<"BARNES1 MATCHES with coeff:  "<<match_map[wild()]<<endl;
    if( full_int_expr.match(tgamma(wild(1)+wild())*tgamma(wild(2)+wild())*tgamma(wild(3)-wild())*tgamma(wild(4)-wild()),match_map))                    
          cout<<"BARNES1 MATCHES wo coeff:  "<<match_map<<endl;
  }
  void barnes2()
  {
    exmap match_map;
    if( full_int_expr.match(wild(6)*tgamma(wild(1)+wild())*tgamma(wild(2)+wild())*tgamma(wild(3)+wild())*tgamma(wild(4)-wild())*tgamma(wild(5)-wild())/tgamma(wild(1)+wild(2)+wild(3)+wild(4)+wild(5)+wild()),match_map))
      cout<<"BARNES2 MATCHES with coeff:  "<<match_map<<endl;
if( full_int_expr.match(tgamma(wild(1)+wild())*tgamma(wild(2)+wild())*tgamma(wild(3)+wild())*tgamma(wild(4)-wild())*tgamma(wild(5)-wild())/tgamma(wild(1)+wild(2)+wild(3)+wild(4)+wild(5)+wild()),match_map))
      cout<<"BARNES2 MATCHES:wo coeff  "<<match_map<<endl;

  }
  /*
  // Constructor for Residue
  MBintegral(lst w_lst_in,lst pole_lst_in,ex full_int_expr_in,exmap w_current_in,relational eps_in):w_lst(w_lst_in),gamma_poles(pole_lst_in),full_int_expr(full_int_expr_in),w_current(w_current_in),eps_current(eps_in),tree_level(0)
  {
 
  }
*/

  MBintegral()
    {
    }
  
  // Constructor for Residue
 MBintegral(lst w_lst_in,ex full_int_expr_in,exmap w_current_in,relational eps_in,size_t tree_lvl):full_int_expr(full_int_expr_in),w_current(w_current_in),eps_current(eps_in),tree_level(tree_lvl),res_pole(0),opt_flag(false)
  {
    w_lst.assign(w_lst_in.begin(),w_lst_in.end());
    //update_poles_from_ex();
    exset new_poles = poles_from_ex(full_int_expr_in);
    gamma_poles.assign(new_poles.begin(),new_poles.end());
    
  }
  
  // constructor for two MBintegrals concatination
  MBintegral(lst w_lst_in,lst pole_lst_in,ex full_int_expr_in):full_int_expr(full_int_expr_in),tree_level(0),res_pole(0),opt_flag(false)
  {
    w_lst.assign(w_lst_in.begin(),w_lst_in.end());
    gamma_poles.assign(pole_lst_in.begin(),pole_lst_in.end());
  }
 

  MBintegral(UFXmap,lst,numeric,bool, unsigned int displacement = 0); // lst nu is a list of powers of propagators and l is a number of loops
};

typedef std::list<MBintegral> MBlst;
typedef mbtree::tree<MBintegral> MBtree;

MBlst MBcontinue(MBintegral rootint,ex eps0 = 0);
MBtree MBcontinue_tree(MBintegral rootint,ex eps0 = 0);

ex expand_and_integrate(MBintegral& int_in, lst num_subs, int expansion_order = 1); // up to O(eps^1) 

#endif // __MBINTEGRAL_H__
