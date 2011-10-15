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


  exset w_lst;
  exset gamma_poles;
  ex full_int_expr;
  exmap eps_w_current;
  exmap w_current;
  relational eps_current;
  int tree_level;
public:
  typedef std::list<double>::iterator gamma_iterator;
  //  typedef lst::const_iterator pole_iterator;
  typedef exset::iterator pole_iterator;
  typedef exset::iterator w_iterator;

  std::pair<pole_iterator,pole_iterator> gamma_iters()
  {
    return std::make_pair(gamma_poles.begin(),gamma_poles.end());
  }
  std::pair<w_iterator,w_iterator> w_iters()
  {
    return std::make_pair(w_lst.begin(),w_lst.end());
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

  lst get_poles();

  exset  get_poles_set()
  {
    return gamma_poles;
  }
 
  lst get_w_lst()
  {
    return set2lst(w_lst);
  }
  exset get_w_set()
  {
    return w_lst;
  }

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

    
    w_lst.insert(to_add.w_iters().first,to_add.w_iters().second);
    cout<<"in add"<<endl;
    gamma_poles.insert(to_add.gamma_iters().first,to_add.gamma_iters().second);
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
    w_lst.insert(w);
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
    //ex tmp_pole(pole);
    //const ex a;
    //    cout<<bind(&MBintegral::pole_has,ref(this),pole,_1)(a);
    //    std::for_each(w_lst.begin(),w_lst.end(), bind(&ex::has, ref(pole), _1));
    //    std::for_each(w_lst.begin(),w_lst.end(), boost::mem_fn(&ex::has));
    //    bind(&ex::has, ref(ex), _1)(i);
    // w_lst

    //inserting only poles with W, and exists in expr
    if (0 < std::count_if(w_lst.begin(),w_lst.end(), bind(&MBintegral::pole_has,ref(this),pole,_1)))
      //        && (full_int_expr.has(tgamma(pole)) || full_int_expr.has(psi(pole)) || full_int_expr.has(psi(wild(),pole)) ) )
      gamma_poles.insert(pole);
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
 MBintegral(lst w_lst_in,ex full_int_expr_in,exmap w_current_in,relational eps_in,size_t tree_lvl):full_int_expr(full_int_expr_in),w_current(w_current_in),eps_current(eps_in),tree_level(tree_lvl)
  {
    w_lst = lst2set(w_lst_in);
    //update_poles_from_ex();
    gamma_poles = poles_from_ex(full_int_expr_in);
    
  }
  
  // constructor for two MBintegrals concatination
  MBintegral(lst w_lst_in,lst pole_lst_in,ex full_int_expr_in):full_int_expr(full_int_expr_in),tree_level(0)
  {
    w_lst = lst2set(w_lst_in);
    gamma_poles = lst2set(pole_lst_in);
  }
 

  MBintegral(FXmap,lst,numeric, unsigned int displacement = 0); // lst nu is a list of powers of propagators and l is a number of loops
};

typedef std::list<MBintegral> MBlst;
typedef mbtree::tree<MBintegral> MBtree;

MBlst MBcontinue(MBintegral rootint,ex eps0 = 0);
MBtree MBcontinue_tree(MBintegral rootint,ex eps0 = 0);

ex expand_and_integrate(MBintegral& int_in, lst num_subs, int expansion_order = 1); // up to O(eps^1) 

#endif // __MBINTEGRAL_H__
