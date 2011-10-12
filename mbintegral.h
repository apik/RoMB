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



#include <boost/icl/interval.hpp>
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/continuous_interval.hpp>
#include <boost/icl/right_open_interval.hpp>
#include <boost/icl/left_open_interval.hpp>
#include <boost/icl/closed_interval.hpp>
#include <boost/icl/open_interval.hpp>
#include <boost/icl/functors.hpp>

#include <iterator>
#include <functional>

#include <gmpxx.h>

#include <ginacra/ginacra.h>

#include "romb_excompiler.h"
#include "collect_square.h"
//#include "sim.h"
#include <cuba.h>

//#include "utils.h"
#include "uf.h"
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


  lst w_lst;
  lst gamma_poles;
  ex full_int_expr;
  exmap eps_w_current;
  exmap w_current;
  relational eps_current;
  int tree_level;
public:
  typedef std::list<double>::iterator gamma_iterator;
  typedef lst::const_iterator pole_iterator;
  std::pair<pole_iterator,pole_iterator> gamma_args()
  {
    return std::make_pair(gamma_poles.begin(),gamma_poles.end());
  }
  lst gamma_args_with_eps();
 
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
  lst get_pole_lst()
  {

    cout<<"get_pole_lst called and ret "<<gamma_poles<<endl;

    //update_poles_from_ex();
    exset gammaset,psiset,psi2set;

    full_int_expr.find(tgamma(wild()),gammaset);
    full_int_expr.find(psi(wild()),psiset);
    full_int_expr.find(psi(wild(),wild(1)),psi2set);
    cout<<" but must gamma "<<gammaset<<endl;
    cout<<" but must psi(ex) "<<psiset<<endl;
    cout<<" but must psi(int,ex) "<<psi2set<<endl;
    
    return gamma_poles;
  }

  void update_poles_from_ex();
 
  lst get_w_lst()
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

  // add terms to w_lst  
  void add_w_lst(lst w_to_add)
  {
    for(lst::const_iterator wtit  = w_to_add.begin();wtit != w_to_add.end(); ++wtit)
      w_lst.append(*wtit);
  }
  // add terms to gamma_poles_lst  
  void add_poles_lst(lst poles_to_add)
  {
    for(lst::const_iterator wtit  = poles_to_add.begin();wtit != poles_to_add.end(); ++wtit)
      gamma_poles.append(*wtit);
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
  // Constructor for Residue
  MBintegral(lst w_lst_in,lst pole_lst_in,ex full_int_expr_in,exmap w_current_in,relational eps_in):w_lst(w_lst_in),gamma_poles(pole_lst_in),full_int_expr(full_int_expr_in),w_current(w_current_in),eps_current(eps_in),tree_level(0)
  {
 
  }
  // constructor for two MBintegrals concatination
  MBintegral(lst w_lst_in,lst pole_lst_in,ex full_int_expr_in):w_lst(w_lst_in),gamma_poles(pole_lst_in),full_int_expr(full_int_expr_in),tree_level(0)
  {
 
  }
 

  MBintegral(FXmap,lst,numeric, unsigned int displacement = 0); // lst nu is a list of powers of propagators and l is a number of loops
};

typedef std::list<MBintegral> MBlst;

MBlst MBcontinue(MBintegral rootint,ex eps0 = 0);

ex expand_and_integrate(MBintegral& int_in, lst num_subs, int expansion_order = 1); // up to O(eps^1) 

#endif // __MBINTEGRAL_H__
