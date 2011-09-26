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

#include "romb_excompiler.h"
#include "collect_square.h"
#include "fm.h"
#include <cuba.h>
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


   boost::mt19937 rng(static_cast<unsigned char>(std::time(0)));       




/// Factory for symbol creation
const symbol & get_symbol(const string & s)
{
  static std::map<string, symbol> directory;
  std::map<string, symbol>::iterator i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, symbol(s))).first->second;
}
/// Enum for easy numeration of UFx_lst_tuple entries
enum {Uarg,Farg,x_lstarg};
/// Tuple for storing U and F polynoms and list of Feynman parameters x_lst
typedef tuple<ex,ex,lst>  UFx_lst_tuple;
namespace UFX
{
  struct U;
  struct F;
  struct xlst;
}

typedef fusion::map<fusion::pair< UFX::U, ex>, fusion::pair< UFX::F, ex>, fusion::pair< UFX::xlst,lst> > UFXmap;
typedef fusion::map< fusion::pair< UFX::F, ex>, fusion::pair< UFX::xlst,lst> > FXmap;



lst set2lst(const exset & s) {
  lst l;
  exset::const_iterator it = s.begin();
  while (it != s.end()) {
    l.append(*it);
    ++it;
  }
  return l;
}

exset lst2set(const lst & l) {
  exset s;
  lst::const_iterator it = l.begin();
  while (it != l.end()) {
    s.insert(*it);
    ++it;
  }
  return s;
}



// return value for 1st variable
std::pair<ex,double>  hyper_cube(lst pole_list,lst w_list)
{
    using namespace lemon;
    Lp lp;
    exhashmap<LpBase::Col> col_map;
    for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
      {
        col_map[*wi] = lp.addCol();
      }
    for(lst::const_iterator pit = pole_list.begin();pit!= pole_list.end();++pit)
      {
        ex tmp_expr = *pit;
        // cout<<"Expr: "<<*pit<<" subexprs:"<<endl;
        Lp::Expr constr_expr;
        for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
          {
            // cout<<*wi<<"  coeff "<<tmp_expr.coeff(*wi)<<endl;
            ex wi_coeff = tmp_expr.coeff(*wi);
            tmp_expr-=(*wi)*wi_coeff;
            if(is_a<numeric>(wi_coeff))
              {
                constr_expr+=ex_to<numeric>(wi_coeff).to_double()*col_map[*wi];
                //                cout<<ex_to<numeric>(wi_coeff).to_double()<<endl; 
              }
            else throw std::logic_error("Non numeric coefficient in pole term. ");
          }
        //constr_expr+=
        //cout<<"Ostatok "<<tmp_expr<<endl;
        if(is_a<numeric>(tmp_expr))
          lp.addRow(-ex_to<numeric>(tmp_expr).to_double(),constr_expr,Lp::INF);
        else throw std::logic_error("Lower bound is not a numeric");
      }

    double l_bound,r_bound;
    cout<<"Hyper cube"<<endl;
    //   for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
      
        lp.min();
        lp.obj(col_map[*w_list.begin()]);
        // Solve the problem using the underlying LP solver
        lp.solve();
        if (lp.primalType() == Lp::OPTIMAL) 
          {
            cout<<*w_list.begin()<<" = ("<<lp.primal()<<",";
            l_bound = lp.primal();
          }
        else throw std::logic_error("Optimal solution not found.");
        lp.max();
        lp.obj(col_map[*w_list.begin()]);
        // Solve the problem using the underlying LP solver
        lp.solve();
        if (lp.primalType() == Lp::OPTIMAL) 
          {
            cout<<lp.primal()<<")"<<endl;
            r_bound = lp.primal();
          }
        else throw  std::logic_error("Optimal solution not found.");

        //  boost::mt19937 rng;       
        boost::uniform_real<> bounded_distribution(l_bound,r_bound);      // distribution that maps to 1..6
        // see random number distributions
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
          die(rng, bounded_distribution);             // glues randomness with mapping
        //rng.seed(static_cast<unsigned char>(std::time(0)));

        
        double x   = die();                      // simulate rolling a die
        // new edition in center of interval no random
        double x_half = (r_bound + l_bound+2.0*x)/4.0;
        double real_half = (r_bound + l_bound)/2.0;
        cout<<"Point "<<die()<<"   "<<x<<endl;
        return std::make_pair(*w_list.begin(),x_half);
        


}



// return value for 1st variable
std::pair<ex,double>  hyper_cube_den(lst pole_list,lst w_list, ex den)
{
    using namespace lemon;
    Lp lp;
    exhashmap<LpBase::Col> col_map;
    for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
      {
        col_map[*wi] = lp.addCol();
      }
    for(lst::const_iterator pit = pole_list.begin();pit!= pole_list.end();++pit)
      {
        ex tmp_expr = *pit;
        // cout<<"Expr: "<<*pit<<" subexprs:"<<endl;
        Lp::Expr constr_expr;
        for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
          {
            // cout<<*wi<<"  coeff "<<tmp_expr.coeff(*wi)<<endl;
            ex wi_coeff = tmp_expr.coeff(*wi);
            tmp_expr-=(*wi)*wi_coeff;
            if(is_a<numeric>(wi_coeff))
              {
                constr_expr+=ex_to<numeric>(wi_coeff).to_double()*col_map[*wi];
                //                cout<<ex_to<numeric>(wi_coeff).to_double()<<endl; 
              }
            else throw std::logic_error("Non numeric coefficient in pole term. ");
          }
        //constr_expr+=
        //cout<<"Ostatok "<<tmp_expr<<endl;
        if(is_a<numeric>(tmp_expr))
          lp.addRow(-ex_to<numeric>(tmp_expr).to_double(),constr_expr,Lp::INF);
        else throw std::logic_error("Lower bound is not a numeric");
      }

    double l_bound,r_bound;
    cout<<"Hyper cube"<<endl;
    //   for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
      
        lp.min();
        lp.obj(col_map[*w_list.begin()]);
        // Solve the problem using the underlying LP solver
        lp.solve();
        if (lp.primalType() == Lp::OPTIMAL) 
          {
            cout<<*w_list.begin()<<" = ("<<lp.primal()<<",";
            l_bound = lp.primal();
          }
        else throw std::logic_error("Optimal solution not found.");
        lp.max();
        lp.obj(col_map[*w_list.begin()]);
        // Solve the problem using the underlying LP solver
        lp.solve();
        if (lp.primalType() == Lp::OPTIMAL) 
          {
            cout<<lp.primal()<<")"<<endl;
            r_bound = lp.primal();
          }
        else throw  std::logic_error("Optimal solution not found.");

        //  boost::mt19937 rng;       
        boost::uniform_real<> bounded_distribution(l_bound,r_bound);      // distribution that maps to 1..6
        // see random number distributions
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
          die(rng, bounded_distribution);             // glues randomness with mapping
        //rng.seed(static_cast<unsigned char>(std::time(0)));

        
        double x   = die();                      // simulate rolling a die
        // new edition in center of interval no random
        double x_half = l_bound*(1- ex_to<numeric>(pow(den,-1)).to_double()) + r_bound*ex_to<numeric>(pow(den,-1)).to_double();
        double real_half = (r_bound + l_bound)/2.0;
        cout<<"Point "<<die()<<"   "<<x<<endl;
        return std::make_pair(*w_list.begin(),x_half);
        


}





// return value for all variables
exmap  hyper_cube_all(lst pole_list,lst w_list)
{
  exmap whmap;
  using namespace lemon;
  Lp lp;
    exhashmap<LpBase::Col> col_map;
    for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
      {
        col_map[*wi] = lp.addCol();
      }
    for(lst::const_iterator pit = pole_list.begin();pit!= pole_list.end();++pit)
      {
        ex tmp_expr = *pit;
        // cout<<"Expr: "<<*pit<<" subexprs:"<<endl;
        Lp::Expr constr_expr;
        for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
          {
            // cout<<*wi<<"  coeff "<<tmp_expr.coeff(*wi)<<endl;
            ex wi_coeff = tmp_expr.coeff(*wi);
            tmp_expr-=(*wi)*wi_coeff;
            if(is_a<numeric>(wi_coeff))
              {
                constr_expr+=ex_to<numeric>(wi_coeff).to_double()*col_map[*wi];
                //                cout<<ex_to<numeric>(wi_coeff).to_double()<<endl; 
              }
            else throw std::logic_error("Non numeric coefficient in pole term. ");
          }
        //constr_expr+=
        //cout<<"Ostatok "<<tmp_expr<<endl;
        if(is_a<numeric>(tmp_expr))
          lp.addRow(-ex_to<numeric>(tmp_expr).to_double(),constr_expr,Lp::INF);
        else throw std::logic_error("Lower bound is not a numeric");
      }

    for(lst::const_iterator lit = w_list.begin(); lit != w_list.end(); ++lit)
      {
        double l_bound,r_bound;
        cout<<"Hyper cube"<<endl;
        //   for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
      
        lp.min();
        lp.obj(col_map[*lit]);
        // Solve the problem using the underlying LP solver
        lp.solve();
        if (lp.primalType() == Lp::OPTIMAL) 
          {
            cout<<*lit<<" = ("<<lp.primal()<<",";
            l_bound = lp.primal();
          }
        else throw std::logic_error("Optimal solution not found.");
        lp.max();
        lp.obj(col_map[*lit]);
        // Solve the problem using the underlying LP solver
        lp.solve();
        if (lp.primalType() == Lp::OPTIMAL) 
          {
            cout<<lp.primal()<<")"<<endl;
            r_bound = lp.primal();
          }
        else throw  std::logic_error("Optimal solution not found.");

        //  boost::mt19937 rng;       
        boost::uniform_real<> bounded_distribution(l_bound,r_bound);      // distribution that maps to 1..6
        // see random number distributions
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
          die(rng, bounded_distribution);             // glues randomness with mapping
        //rng.seed(static_cast<unsigned char>(std::time(0)));

        
        double x   = die();                      // simulate rolling a die
        // new edition in center of interval no random
        double x_half = (r_bound + l_bound+2.0*x)/4.0;
        double real_half = (r_bound + l_bound)/2.0;
        cout<<"Point "<<die()<<"   "<<x<<endl;
        whmap[*lit] = real_half;
      }
    return whmap;
        


}


exmap start_point(lst pole_list,lst w_list)
{
  exmap subs_map;
    for(lst::const_iterator wi = w_list.begin();wi != w_list.end();++wi) 
    {
      lst tmp_pole_list;
      lst tmp_w_list;
      for(lst::const_iterator pit = pole_list.begin();pit!= pole_list.end();++pit)
        if(!((*pit).subs(subs_map)>0))
          tmp_pole_list.append((*pit).subs(subs_map));
      cout<<"tmp_pole_list"<<tmp_pole_list<<endl;
      for(lst::const_iterator wi2 = wi;wi2 != w_list.end();++wi2) 
        tmp_w_list.append(*wi2);
      std::pair<ex,double> ret_pair = hyper_cube(tmp_pole_list,tmp_w_list);
      subs_map[ret_pair.first] = ret_pair.second;
      cout<<ret_pair.first<<" "<<ret_pair.second<<endl;
    }
  
  cout<<"START POINT SUBS  "<<subs_map<<endl;
  return subs_map;
}
// Start point with different countours
exmap start_point_diff_w(lst pole_list,lst w_list)
{
  exset point_set;
  exmap subs_map;
  ex den = 2;
    for(lst::const_iterator wi = w_list.begin();wi != w_list.end();++wi) 
    {
      lst tmp_pole_list;
      lst tmp_w_list;
      for(lst::const_iterator pit = pole_list.begin();pit!= pole_list.end();++pit)
        if(!((*pit).subs(subs_map)>0))
          tmp_pole_list.append((*pit).subs(subs_map));
      cout<<"tmp_pole_list"<<tmp_pole_list<<endl;
      for(lst::const_iterator wi2 = wi;wi2 != w_list.end();++wi2) 
        tmp_w_list.append(*wi2);

      std::pair<ex,double> ret_pair = hyper_cube_den(tmp_pole_list,tmp_w_list,den);
      if(point_set.count(ret_pair.second) > 0)
	do
	  {
	    den +=1;
	    ret_pair = hyper_cube_den(tmp_pole_list,tmp_w_list,den);
	  }
	while(point_set.count(ret_pair.second) > 0);
      subs_map[ret_pair.first] = ret_pair.second;
      point_set.insert(ret_pair.second);
      cout<<ret_pair.first<<" "<<ret_pair.second<<endl;
    }
  
  cout<<"START POINT SUBS  "<<subs_map<<endl;
  return subs_map;
}



struct compare_xterm_varieties : public std::binary_function<ex,ex,bool>
{
  lst x_lst;
  compare_xterm_varieties(lst x_lst_) : x_lst(x_lst_)
  {
  }
  bool operator()(ex left_ex,ex right_ex)
  {
    unsigned int max_left,max_right;
    max_left = 0;
    max_right = 0;
    for(lst::const_iterator it = x_lst.begin(); it != x_lst.end(); ++it)
      {
        if(left_ex.has(*it))max_left++;
        if(right_ex.has(*it))max_right++;
      }
    return (max_right>max_left);
  }
};

/* *
   Lexicographical ordering of polynoms
   1) overall power 
   2) maximum power
   3) lexi summ
   !change if true!
 */
struct comp_ex_xpow 
{
  lst x_lst;
  comp_ex_xpow(lst x_lst_) : x_lst(x_lst_)
  {
  }
  bool operator()(ex left_ex,ex right_ex)
  {
    unsigned max_left,max_right;
    max_left = 0;
    max_right = 0;
    unsigned op_l = 0,op_r = 0;  // overall power of monomial
    unsigned mp_l,mp_r;// maximum power in monomial
    unsigned li_l = 0,li_r = 0;       // lexi index
    for(lst::const_iterator it = x_lst.begin(); it != x_lst.end(); ++it)
      {
        if(left_ex.has(*it))
          {
            li_l+=std::distance(x_lst.begin(),it);
            max_left++;
          }
        if(right_ex.has(*it))
          {
            max_right++;
            li_r+=std::distance(x_lst.begin(),it);
          }
      }
    if(max_right<max_left)return true;
    else if(max_right ==max_left)return (li_r<li_l);
    else return false;
  }
};

lst bubble_sort_lexi(lst in_lst,lst x_lst)
{
  // instance of sorter
  comp_ex_xpow term_comparator(x_lst);
  bool have_changes = true;
  ex buf;
  while(have_changes==true)
    {
      have_changes = false;
      for(int cntr = 0; cntr<in_lst.nops(); ++cntr)
        {
          if(cntr+1<in_lst.nops())
            {
              if(term_comparator(in_lst.op(cntr),in_lst.op(cntr+1)))
                {
                  // swap(it,boost::next(it));
                  buf = in_lst.op(cntr+1);
                  in_lst.let_op(cntr+1) = in_lst.op(cntr);
                  in_lst.let_op(cntr) = buf;
                  have_changes = true;
                }
            }
        }
      
    }
  return in_lst;
}

class MBintegral
{
  int tree_level;
  ex full_int_expr;
  lst gamma_poles;
  lst w_lst;
  exmap eps_w_current;
  exmap w_current;
  relational eps_current;
public:
  typedef std::list<double>::iterator gamma_iterator;
  typedef lst::const_iterator pole_iterator;
  std::pair<pole_iterator,pole_iterator> gamma_args()
  {
    return std::make_pair(gamma_poles.begin(),gamma_poles.end());
  }
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
    return gamma_poles;
  }
  lst get_w_lst()
  {
    return w_lst;
  }
  exmap new_point()
  {
    lst var_list(w_lst);
    var_list.append(get_symbol("eps"));
    eps_w_current=start_point_diff_w(get_pole_lst(),var_list);
    //eps_w_current=findinstance(get_pole_lst(),var_list);
    for(lst::const_iterator it = w_lst.begin();it!=w_lst.end();++it)
      w_current[*it] = it->subs(eps_w_current);
    eps_current = (get_symbol("eps")==eps_w_current[get_symbol("eps")]);
    return eps_w_current;
  }

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
  MBintegral res(relational w_relation,ex pole,relational new_eps)
  {
    try
      {
        lst new_w_lst;
        for(lst::const_iterator it = w_lst.begin();it!=w_lst.end();++it)
          if(*it != w_relation.lhs()) new_w_lst.append(*it);
        //    std::remove_copy(w_lst.begin(),w_lst.end(),new_w_lst.begin(),w_relation.lhs());
        cout<<"removed "<<w_relation.lhs()<<"  in list "<<new_w_lst<<endl;
        lst new_gamma_pole_list;
        for(lst::const_iterator it = gamma_poles.begin();it!=gamma_poles.end();++it)
          if(*it != pole) new_gamma_pole_list.append(it->subs(w_relation));
        cout<<"removed pole list  "<<pole.subs(w_relation)<<"  in list "<<new_gamma_pole_list<<endl;
        // !!! IMPORTANT!!! no 2*pi*i multiplication and no sign multiplication 
        ex new_no_gamma_part = (full_int_expr.subs(tgamma(pole)==pow(-1,-pole.subs(w_relation))/factorial(-pole.subs(w_relation)))).subs(w_relation);
        // new_no_gamma_part  = pow(-1,pole.subs(w_relation))/factorial(pole.subs(w_relation))*full_int_expr.subs(w_relation);
        cout<< new_no_gamma_part<<endl;
        exmap new_w_current(w_current);
        cout<<" Not modif:  "<<new_w_current<<endl;
        new_w_current.erase(w_relation.lhs());
        cout<<" modif:  "<<new_w_current<<endl;
        MBintegral resINT(new_w_lst,new_gamma_pole_list,new_no_gamma_part,new_w_current,new_eps);
        return resINT;
      }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"MBintegral.res\":\n |___> ")+p.what());
      }
  }

   ex get_expr()
  {
    return full_int_expr; 
  }

  ex get_gamma_expr()
  {
    ex tmpex = full_int_expr;
    for(pole_iterator it = gamma_poles.begin();it!=gamma_poles.end();++it)
      {
        tmpex*=tgamma(*it);
      }
    return tmpex;
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
          cout<<"BARNES1 MATCHES with coeff:  "<<match_map<<endl;
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

  MBintegral(UFXmap ufx_in,lst nu,numeric l):tree_level(0) // lst nu is a list of powers of propagators and l is a number of loops
  {
    try
      {
        eps_current = (get_symbol("eps")==0);
        ex N = accumulate(nu.begin(),nu.end(),ex(0));
        //    cout<<"N(nu)= "<<N<<endl;
        lst  x_lst(fusion::at_key<UFX::xlst>(ufx_in));
        ex U = fusion::at_key<UFX::U>(ufx_in).collect(x_lst,true);
        ex F = fusion::at_key<UFX::F>(ufx_in).collect(x_lst,true);///< distributed polynom factorizing X(i)X(j)*(...)


        // assuming 1/(U^a*F^b)
        ex U_pow = (-N+(l+1)*(2-get_symbol("eps")));
        cout<<"U_pow"<<U_pow<<endl;
        ex F_pow = (N-l*(2-get_symbol("eps")));
    
        cout<<setw(30)<<std::internal<<"+++INTEGRAL PARAMETERS+++"<<endl;
        cout<<setw(30)<<std::left<<"** Number of loops L=   "<<std::internal<<l<<endl;
        cout<<setw(30)<<std::left<<"** Summ of powers N=   "<<N<<endl;
        //cout<<setw(35)<<std::left<<"** Dimension D=   "<<D.subs(D_subs)<<endl;
        cout<<setw(35)<<std::left<<"** U =   "<<U<<endl;
        cout<<setw(35)<<std::left<<"** U power=   "<<U_pow<<endl;
        cout<<setw(35)<<std::left<<"** F =   "<<F<<endl;
        cout<<setw(35)<<std::left<<"** F power=   "<<F_pow<<endl;
        cout<<setw(30)<<std::internal<<"+++++++++++++++++++++++++"<<endl;
        //lst gamma_lst; //Gamma with poles left of contour
        //    lst w_lst;
        //  ex gamma_right; // and right of contour
        exmap x_power_map; // map of X(j) powers
        for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
          {
            x_power_map[*it1] = nu.op(std::distance(x_lst.begin(),it1))-1;
          }

        cout<<"x_map_start "<<x_power_map<<endl;

        //    ex coeff = 1;                             // numerical coeeficient independent of X(j)
        ex coeff = pow(exp(get_symbol("eps")*Euler)*pow(get_symbol("mu"),2*get_symbol("eps")),l)*tgamma(F_pow);
        // important if power = 0????
        for(lst::const_iterator nui = nu.begin();nui!=nu.end();++nui)
          coeff/=tgamma(*nui);
        coeff/=tgamma(F_pow);
        //    ex out_ex = 1;///pow(2*Pi*I,U.nops()-1)/tgamma(al_pow)/pow(2*Pi*I,F.nops()-1)/tgamma(al_pow); //need review
        //working with F-term \Gamma(\nu-L*D/2) contractedx
        ex w_sum = 0;  //F-term generates only integrations in W
        for(const_iterator it = F.begin();it!=F.end();++it)
          {
            int w_index = distance(F.begin(),it);
            ex x_power;
            if(F.end()==boost::next(it)) 
              {
                coeff*=tgamma(F_pow+w_sum);
                x_power = -F_pow-w_sum;
                gamma_poles.append(F_pow+w_sum);
                cout<<"end achieved"<<endl;
              }
            else
              {
                string str = "w_"+boost::lexical_cast<string>(w_index);
                //symbol w(str);
                w_lst.append(get_symbol(str));
                coeff*=tgamma(-get_symbol(str));
                x_power = get_symbol(str);
                gamma_poles.append(-get_symbol(str));
                w_sum+=get_symbol(str);
                cout<<"ok run"<<endl;
              }
            // filling map of X(j) powers
            ex tmp_expr = (*it);
            for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
              {
                if(tmp_expr.has(*it1))
                  {
                    x_power_map[*it1]+=(tmp_expr.degree(*it1)*x_power);
                    tmp_expr = tmp_expr.subs((*it1)==1);
                  }
              }
            // coeff *=pow(tmp_expr,x_power);
          }
    
        //working with U-term, only if (U_pow>0)
        
        if(U_pow.subs(get_symbol("eps")==0) > 0)
          {
            throw std::logic_error("Unsupported topology too small number of propagators versus large number of vertices");
            cout<<">>> Working with MB for U term, U_power = "<< U_pow<<endl;  
            ex v_sum = 0;  //F-term generates only integrations in V
            for(const_iterator it = U.begin();it!=U.end();++it)
              {
                int v_index = distance(U.begin(),it);
                ex x_power;
                if(U.end()==boost::next(it)) 
                  {
                    coeff*=tgamma(U_pow+v_sum);
                    x_power = -U_pow - v_sum;
                    gamma_poles.append(U_pow+v_sum);
                    cout<<"end achieved"<<endl;
                  }
                else
                  {
                    string str = "v_"+boost::lexical_cast<string>(v_index);
                    //symbol v(str);
                    w_lst.append(get_symbol(str));
                    coeff*=tgamma(-get_symbol(str));
                    x_power = get_symbol(str);
                    gamma_poles.append(-get_symbol(str));
                    v_sum+=get_symbol(str);
                    cout<<"ok run"<<endl;
                  }
                // filling map of X(j) powers
                ex tmp_expr = (*it);
                for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
                  {
                    if(tmp_expr.has(*it1))
                      {
                        x_power_map[*it1]+=(tmp_expr.degree(*it1)*x_power);
                        tmp_expr = tmp_expr.subs((*it1)==1);
                      }
                  }
                //coeff *=pow(tmp_expr,x_power);
              }
          }
        else
          {
            /// IMPORTANT add function for arbitrary power of U!!!!
            ex U_term = pow(U,U_pow.subs(get_symbol("eps")==0));
            for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
              if(U_term.has(*it1))
                {
                  x_power_map[*it1]+=((-1)*U_term.degree(*it1));
                  cout<< "U_TERM  "<<U_term<<" deg of "<<*it1<<" is "<<((-1)*U_term.degree(*it1))<<endl;

                }
        
          }
    
        //cout<<x_power_map<<endl;
        cout<<"Gammas after MB: "<<endl<<gamma_poles<<endl;
        cout<<"X powers list:"<<endl<<x_power_map<<endl;
        // applying X integration
        ex gamma_den = 0; // gamma in denominator
        for(exmap::const_iterator mi = x_power_map.begin();mi!=x_power_map.end();++mi)
          {
            cout<<(*mi).first<<" "<<(*mi).second<<endl;
            gamma_poles.append((*mi).second+1);
            coeff*=tgamma((*mi).second+1);
            gamma_den+=((*mi).second+1);
          }
        cout<<"GAMMA_DEN: "<<gamma_den<<endl;
        coeff/=tgamma(gamma_den);
        cout<<"New gamma list:"<<endl<<gamma_poles<<endl;
        full_int_expr = coeff;
        cout<<w_lst<<endl;
      }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"MBintegral)UFX)\":\n |___> ")+p.what());
      }

  }
 


/**
 *
 *  Construction from F, U=1
 *
 */
  MBintegral(FXmap fx_in,lst nu,numeric l, unsigned int displacement = 0):tree_level(0) // lst nu is a list of powers of propagators and l is a number of loops
  {
    try
      {
        eps_current = (get_symbol("eps")==0);
        ex N = accumulate(nu.begin(),nu.end(),ex(0));

        //    cout<<"N(nu)= "<<N<<endl;
        lst  x_lst(fusion::at_key<UFX::xlst>(fx_in));
        // summ of feynman parameters satisfying summ(X_i)=1
        //        ex x_summ = accumulate(x_lst.begin(),x_lst.end(),ex(0));

        //relational delta_subs(x_lst.op(x_lst.nops()-1),lsolve(x_summ==1,x_lst.op(x_lst.nops()-1)));// relation delta function
        // cout<<"xsumm "<<delta_subs.lhs()<<" == "<<delta_subs.rhs()<<endl;
        
        //ex F = fusion::at_key<UFX::F>(fx_in).collect(x_lst,true);///< distributed polynom factorizing X(i)X(j)*(...)
        ex F = (fusion::at_key<UFX::F>(fx_in)).expand();///< distributed polynom factorizing X(i)X(j)*(...)
        //F=F.subs(delta_subs);//applying delta function relation on F-polynom
        //x_lst.remove_last();

        // assuming 1/(U^a*F^b)
        
        ex F_pow = (N-l*(2-get_symbol("eps")));
	exset f_set;
	if(F.find(wild(1)*wild(1)+2*wild(1)*wild(2)+wild(2)*wild(2), f_set))cout<< "HAVE square"<<endl;
        cout<<setw(30)<<std::internal<<"+++INTEGRAL PARAMETERS+++"<<endl;
        cout<<setw(30)<<std::left<<"** Number of loops L=   "<<std::internal<<l<<endl;
        cout<<setw(30)<<std::left<<"** Summ of powers N=   "<<N<<endl;
        //cout<<setw(35)<<std::left<<"** Dimension D=   "<<D.subs(D_subs)<<endl;
        cout<<setw(35)<<std::left<<"** F =   "<<F<<endl;
        cout<<setw(35)<<std::left<<"** F power=   "<<F_pow<<endl;
        cout<<setw(30)<<std::internal<<"+++++++++++++++++++++++++"<<endl;
        
        	lst coe_l,xsq_l;
        	tie(coe_l,xsq_l) = collect_square(F,x_lst);
        	cout<<">>> Found "<<coe_l.nops()<<" full squares in F polynomial"<<endl;
        	cout<< coe_l<<" * "<<xsq_l<<endl;
        	cout<<" F qad " <<F<<endl;
        	
        	F = F.collect(x_lst,true);
		//	assert(false);

        //lst gamma_lst; //Gamma with poles left of contour
        //    lst w_lst;
        //  ex gamma_right; // and right of contour
        exmap x_power_map; // map of X(j) powers
        for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
          {
            x_power_map[*it1] = nu.op(std::distance(x_lst.begin(),it1))-1;
          }

        cout<<"x_map_start "<<x_power_map<<endl;

		
        //    ex coeff = 1;                             // numerical coeeficient independent of X(j)
        ex coeff =tgamma(F_pow)* pow(exp(get_symbol("eps")*Euler),l);///pow(I*pow(Pi,2 - get_symbol("eps")),l);
        // important if power = 0????
        for(lst::const_iterator nui = nu.begin();nui!=nu.end();++nui)
          coeff/=tgamma(*nui);
        coeff/=tgamma(F_pow);
        //    ex out_ex = 1;///pow(2*Pi*I,U.nops()-1)/tgamma(al_pow)/pow(2*Pi*I,F.nops()-1)/tgamma(al_pow); //need review

        //working with F-term \Gamma(\nu-L*D/2) contractedx
        ex w_sum = 0;  //F-term generates only integrations in W

	//--------------------------------
	//      MB for full squares
	//--------------------------------
	size_t z_idx = 0;
	if(F.nops() == 0 && xsq_l.nops() == 1)
	{
	/*
	 F is a full square
	*/
	F_pow *=2;
	F = xsq_l.op(0);
	coeff /= coe_l.op(0);
	}
	else if(F.nops() == 0 && xsq_l.nops() > 1)
	{
	/*
	 F is a summ of full squares
	 */
	 throw std::logic_error(std::string("F = (xs1)^2+ ..(xsn)^2; not realized "));
	}
	else
	  {
	    for(int i = 0; i < xsq_l.nops(); i++)
	      {
		string str = "w_"+boost::lexical_cast<string>(displacement);
		displacement++;
                //symbol w(str);
		ex w_i = get_symbol(str);
                w_lst.append(w_i);
                coeff*=tgamma(-w_i)*pow(coe_l.op(i),w_i)/(2*Pi*I);
                
                cout<<"w_i_power "<<w_i<<endl;
                gamma_poles.append(-w_i); //!!!! review
                w_sum+=w_i;

		// subMB construction
		ex sq_lst(xsq_l.op(i));
		cout<<"SQLST : "<<sq_lst<<endl;
		coeff /= (pow(2*Pi*I,sq_lst.nops()-1)*tgamma(-2*w_i));
		gamma_poles.append(-2*w_i);
		ex z_sum = 0;
		for(const_iterator x_it = sq_lst.begin(); x_it != sq_lst.end(); ++x_it)
		  {
		    ex a_power;
		    if(sq_lst.end() == boost::next(x_it)) // X_k expr
		      {
			coeff *= tgamma(-2*w_i + z_sum);
			a_power = 2*w_i - z_sum;
			gamma_poles.append(-2*w_i + z_sum);
		      }
		    else // ordinary exprs
		      {
			string z_str = "z_"+boost::lexical_cast<string>(z_idx);
			z_idx++;
			ex z_i = get_symbol(z_str);
			cout << z_i<<endl;
			w_lst.append(z_i);
			coeff*=tgamma(-z_i);
			a_power = z_i;
			gamma_poles.append(-z_i); //!!!! review
			z_sum += a_power;
		      } 
		    
		    // add x-part with it's coefficient
		    for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
		      {
			if(x_it->has(*it1))
			  {
			    // simple expression x_i or -x_i, x_it->degree(x) == 1
			    BOOST_ASSERT_MSG(x_it->degree(*it1) == 1,"Not a simple expression in square");
			    x_power_map[*it1]+=(a_power);
			    coeff *=pow(x_it->lcoeff(*it1),a_power);
			  }
		      }
		  }
	      }
	  }
        lst F_to_lst;
        if(is_a<add>(F))
          for(const_iterator it = F.begin();it!=F.end();++it)
            F_to_lst.append(*it);
        else 
          F_to_lst.append(F);
       
       
          //        F_to_lst.sort();
        
        //        comp_ex_xpow F_term_comparator(x_lst);
        //std::sort(F_to_lst.begin(),F_to_lst.end(),F_term_comparator);
        //F_to_lst = bubble_sort_lexi(F_to_lst,x_lst);

	cout<<"displace:  "<<displacement<<endl;
        cout<<"FTLST" <<F_to_lst<<endl;
        for(lst::const_iterator it = F_to_lst.begin();it!=F_to_lst.end();++it)
          {
            cout<<"F_term : "<<(*it)<<endl;
            int w_index = std::distance(F_to_lst.begin(),it);
            ex x_power;
            if(F_to_lst.end()==boost::next(it)) 
              {
                coeff*=tgamma(F_pow+w_sum);
                x_power = -F_pow-w_sum;
                cout<<"x_power last term "<<x_power<<"  w_sum "<<w_sum<<endl;
                gamma_poles.append(F_pow+w_sum);
                cout<<"end achieved"<<endl;
              }
            else
              {
                string str = "w_"+boost::lexical_cast<string>(displacement + w_index);
                //symbol w(str);
                w_lst.append(get_symbol(str));
                coeff*=tgamma(-get_symbol(str));
                x_power = get_symbol(str);
                cout<<"x_power "<<x_power<<endl;
                gamma_poles.append(-x_power); //!!!! review
                w_sum+=x_power;
                cout<<"ok run"<<endl;
              }
            // filling map of X(j) powers
            ex tmp_expr = (*it);
            for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
              {
                if(tmp_expr.has(*it1))
                  {
                    x_power_map[*it1]+=(tmp_expr.degree(*it1)*x_power);
                    cout<<"before subs "<<(*it1)<<"   "<<tmp_expr<<endl;
                    tmp_expr = tmp_expr.subs((*it1)==1);
                    //                    cout<<"after subs "<<tmp_expr<<endl;
                  }
              }
            //     cout<<"after subs "<<tmp_expr<<endl;
             coeff *=pow(tmp_expr,x_power);
          }
    
        //working with U-term, only if (U_pow>0)
        
       
        //cout<<x_power_map<<endl;
        cout<<"Gammas after MB: "<<endl<<gamma_poles<<endl;
        cout<<"X powers list:"<<"  "<<x_power_map<<endl;
        //!!!!!!!!!!!!!!!!!!
        //        assert(false);
        // applying X integration
        ex gamma_den = 0; // gamma in denominator
        for(exmap::const_iterator mi = x_power_map.begin();mi!=x_power_map.end();++mi)
          {
            cout<<(*mi).first<<" "<<(*mi).second<<endl;
            gamma_poles.append((*mi).second+1);
            coeff*=tgamma((*mi).second+1);
            gamma_den+=((*mi).second+1);
          }
        cout<<"GAMMA_DEN: "<<gamma_den<<endl;
        bool gamma_den_has_w = false;
        for(lst::const_iterator wit = w_lst.begin(); wit != w_lst.end(); ++wit)
          if(gamma_den.has(*wit))gamma_den_has_w = true;
        if(gamma_den_has_w) gamma_poles.append(gamma_den);
        coeff/=tgamma(gamma_den);
        coeff/=(pow(2*Pi*I,w_lst.nops()));
        cout<<"New gamma list:"<<endl<<gamma_poles<<endl;
        full_int_expr = coeff;
        cout<<w_lst<<endl;
      }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"MBintegral)UFX)\":\n |___> ")+p.what());
      }

  }
};


lst has_w(ex gamma_arg,lst w_list)
{
  lst out_list;
  for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
    if(gamma_arg.has(*wi))out_list.append(*wi);
  return out_list;
}


void hit_pole(lst pole_list,lst w_list,exmap subs_map)
{
  lst C_list;
  exmap C_subs_map;
  for(lst::const_iterator cit = w_list.begin();cit!=boost::prior(w_list.end()); ++cit)
    {
      C_list.append(*cit);
      C_subs_map[*cit] = subs_map[*cit];
    }

  ex eps0 = *w_list.rbegin();


  cout<<eps0.subs(subs_map)<<"  "<<boost::math::iround(ex_to<numeric>(eps0.subs(subs_map)).to_double())<<endl;
  cout<<boost::math::iround(3.7)<<" 3.7 "<<boost::math::iround(-3.7)<< " -3.7"<<endl;
  cout<<"poles closest to C "<<endl;
  for(lst::const_iterator pit = C_list.begin();pit!=C_list.end();++pit)
    {
      int ipart = boost::math::iround(ex_to<numeric>((*pit).subs(subs_map)).to_double());
      double fpart = ex_to<numeric>((*pit).subs(subs_map)).to_double();
      double delta = std::abs(fpart - ipart);

      cout<< "("<<*pit<<","<<ipart<<")   "<<" float: "<<fpart<<" ipart: "<<ipart<<" delta: "<<delta<<endl;
      

    }
  cout<<endl;
  ex eps_max = 0;

  for(lst::const_iterator ui = pole_list.begin();ui!=pole_list.end();++ui)
    {
      lst U_w_list(has_w(*ui,C_list));
      if((U_w_list.nops()>0)&&(ui->has(eps0)))
        {
          cout<<U_w_list<<endl;
          std::vector<double> delta_vector;
          double min_delta = 1000;
          int U_min;
          for(lst::const_iterator wi=U_w_list.begin();wi!=U_w_list.end();++wi)
            if(std::abs(ex_to<numeric>((*wi).subs(subs_map)).to_double() - boost::math::iround(ex_to<numeric>((*wi).subs(subs_map)).to_double()))<min_delta) 
              {
                min_delta = std::abs(ex_to<numeric>((*wi).subs(subs_map)).to_double() - boost::math::iround(ex_to<numeric>((*wi).subs(subs_map)).to_double()));
                U_min = boost::math::iround(ex_to<numeric>((*wi).subs(subs_map)).to_double());
              }
          cout<<"n(i)="<<U_min<<endl;
          ex tmp_eps = lsolve((*ui).subs(C_subs_map) == U_min,eps0 );
          if(abs(tmp_eps) > abs(eps_max))eps_max = tmp_eps;
          cout<<"EPS1 = "<< tmp_eps<<endl;
          
        }
    }

  cout<<"EPS(1) = "<<eps_max<<endl;
  lst u_n_eqs;
  for(lst::const_iterator ui = pole_list.begin();ui!=pole_list.end();++ui)
    {
      // pole value in new EPS point:
      ex U_pole = (ui->subs(C_subs_map)).subs(eps0 == eps_max);
      cout<<*ui<<"         "<<U_pole<<endl;
      if(U_pole==0) u_n_eqs.append(*ui ==0);
      else if(U_pole==-1) u_n_eqs.append(*ui ==-1);
      else if(U_pole==-2) u_n_eqs.append(*ui ==-2);
      else if(U_pole==-3) u_n_eqs.append(*ui ==-3);
      else if(U_pole==-4) u_n_eqs.append(*ui ==-4);
      else if(U_pole==-5) u_n_eqs.append(*ui ==-5);
    }
  // IF P{...}!=0
  if(u_n_eqs.nops()>0)
    {
      lst w_in_eq;
      for(lst::const_iterator eqi = u_n_eqs.begin();eqi!=u_n_eqs.end();++eqi)
        {
          lst wlst = has_w(eqi->lhs(),C_list);
          cout<<wlst<<endl;
          
          for(lst::const_iterator wi = wlst.begin();wi!=wlst.end();++wi)
          w_in_eq.append(*wi);
        }
      w_in_eq.unique();
      ex w_sub_solve = lsolve(u_n_eqs,w_in_eq);
for(lst::const_iterator ui = pole_list.begin();ui!=pole_list.end();++ui)
      cout<<ui->subs(w_sub_solve).subs(C_subs_map).subs(eps0 == eps_max)<<endl;
    }
}

/*
void gamma_lp(lst pole_list,lst w_list)
{
using namespace lemon;
    Lp lp;
    exhashmap<LpBase::Col> col_map;
    for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
      {
        col_map[*wi] = lp.addCol();
      }
    for(lst::const_iterator pit = pole_list.begin();pit!= pole_list.end();++pit)
      {
        ex tmp_expr = *pit;
        // cout<<"Expr: "<<*pit<<" subexprs:"<<endl;
        Lp::Expr constr_expr;
        for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
          {
            // cout<<*wi<<"  coeff "<<tmp_expr.coeff(*wi)<<endl;
            ex wi_coeff = tmp_expr.coeff(*wi);
            tmp_expr-=(*wi)*wi_coeff;
            if(is_a<numeric>(wi_coeff))
              {
                constr_expr+=ex_to<numeric>(wi_coeff).to_double()*col_map[*wi];
                //                cout<<ex_to<numeric>(wi_coeff).to_double()<<endl; 
              }
            else throw std::logic_error("Non numeric coefficient in pole term. ");
          }
        //constr_expr+=
        //cout<<"Ostatok "<<tmp_expr<<endl;
        if(is_a<numeric>(tmp_expr))
          lp.addRow(-ex_to<numeric>(tmp_expr).to_double(),constr_expr,Lp::INF);
        else throw std::logic_error("Lower bound is not a numeric");
      }

    cout<<"Hyper cube"<<endl;
    for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
      {
        lp.max();
        lp.obj(col_map[*wi]);
        // Solve the problem using the underlying LP solver
        lp.solve();
        if (lp.primalType() == Lp::OPTIMAL) 
          cout<<*wi<<" = ("<<lp.primal()<<",";
        lp.min();
        lp.obj(col_map[*wi]);
        // Solve the problem using the underlying LP solver
        lp.solve();
        if (lp.primalType() == Lp::OPTIMAL) 
          cout<<lp.primal()<<")"<<endl;
      }

    cout<<"Maximal eps value"<<endl;
  // Specify the objective function
  lp.max();
  lp.obj(col_map[*w_list.rbegin()]);
  
  // Solve the problem using the underlying LP solver
  lp.solve();

  // Print the results
  if (lp.primalType() == Lp::OPTIMAL) 
    {
std::cout << "Objective function value: " << lp.primal() << std::endl;
      for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
        {
            std::cout <<*wi<< " = " << lp.primal(col_map[*wi]) << std::endl;
        }
    } 
  else throw std::logic_error( "Optimal solution not found.");

    cout<<"Minimal eps value"<<endl;
  // Specify the objective function
  lp.min();
  lp.obj(col_map[*w_list.rbegin()]);
  
  // Solve the problem using the underlying LP solver
  lp.solve();

  // Print the results
  if (lp.primalType() == Lp::OPTIMAL) 
    {
std::cout << "Objective function value: " << lp.primal() << std::endl;
      for(lst::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
        {
            std::cout <<*wi<< " = " << lp.primal(col_map[*wi]) << std::endl;
        }
    } 
  else throw std::logic_error( "Optimal solution not found.");
  

  
}
*/

UFXmap UF(lst k_lst,lst p_lst,lst subs_lst, unsigned int displacement = 0)
{
  /// constructing expr Xi*Pi
  ex fprop = 0;
  ex sDtmp;
  string str;
  lst xp;
  for(lst::const_iterator it = p_lst.begin();it!=p_lst.end();++it)
    {
      str = "x_" + boost::lexical_cast<string>(displacement+std::distance(p_lst.begin(),it));
      xp.append(get_symbol(str)); // storing list of X identities
      fprop += get_symbol(str) * (*it); 
    }
  fprop = fprop.expand();
  sDtmp = fprop;
  //cout<<fprop<<endl;
  matrix M(k_lst.nops(),k_lst.nops());
  matrix Q(k_lst.nops(),1);//column
  GiNaC::numeric half(1,2);
  for(lst::const_iterator itr = k_lst.begin();itr!=k_lst.end();++itr)
    for(lst::const_iterator itc = itr;itc!=k_lst.end();++itc)
      if(itr == itc)
        {
          M(distance(k_lst.begin(),itr),distance(k_lst.begin(),itc)) = -1*fprop.coeff((*itr)*(*itc),1);
          //cout<<(*itr)*(*itc)<<"M("<<distance(k_lst.begin(),itr)<<","<<distance( k_lst.begin(),itc)<<") coeff "<<fprop.coeff((*itr)*(*itc),1)<<endl;
          sDtmp -= (*itr)*(*itc)*fprop.coeff((*itr)*(*itc),1);
        }
      else
        {
          M(distance(k_lst.begin(),itr),distance( k_lst.begin(),itc)) = -1*half*fprop.coeff((*itr),1).coeff((*itc),1);
          M(distance(k_lst.begin(),itc),distance(k_lst.begin(),itr)) = -1*half*fprop.coeff((*itr),1).coeff((*itc),1);
          //cout<<(*itr)*(*itc)<<"M("<<distance( k_lst.begin(),itr)<<","<<distance( k_lst.begin(),itc)<<") coeff "<<fprop.coeff((*itr),1).coeff((*itc),1)<<endl;
          sDtmp -= (*itr)*(*itc)*fprop.coeff((*itr),1).coeff((*itc),1);
        }
  //cout<<"M: "<<M<<endl;
  sDtmp = sDtmp.expand();
  //cout<<"Expr linear on external momentum: "<<sDtmp.expand()<<endl;

  for(lst::const_iterator itr = k_lst.begin();itr!=k_lst.end();++itr)
    {
      Q(distance( k_lst.begin(),itr),0) = half*sDtmp.coeff((*itr),1);
      sDtmp -= (*itr)*sDtmp.coeff((*itr),1);
    }
  //  cout<<"Q: "<<Q<<endl;
  sDtmp = sDtmp.expand();
  ex minusJ = sDtmp;
  // cout<<"-J: "<<minusJ<<endl;
  ex U = M.determinant();
  ex F = expand(M.determinant()*(minusJ+Q.transpose().mul(M.inverse()).mul(Q)(0,0)));
  lst lp;
  F=F.normal();
  cout<<"U= "<<U<<endl<<"F= "<<F<<endl;
  //  cout<<"pol_list "<<lp<<endl;
  U=U.subs(subs_lst,subs_options::algebraic);
  F=F.subs(subs_lst,subs_options::algebraic);
  //  cout<<"UF_WORK"<<endl;
  return fusion::make_map<UFX::U,UFX::F,UFX::xlst>(U,F,xp);
}


/*
ex D;
//symbol eps("eps");
/// list variant
ex MB_lst(pair<lst,lst> UF_x_lst,lst nu,numeric l,relational D_subs=D==4-2*get_symbol("eps"))
{
  D = D.subs(D_subs);
  cout<<D<<endl;
  ex N = accumulate(nu.begin(),nu.end(),ex(0));
  cout<<"N(nu)= "<<N<<endl;
  lst  x_lst(UF_x_lst.second);
  ex U = UF_x_lst.first.op(0).subs(D_subs);
  ex F = UF_x_lst.first.op(1).subs(D_subs).collect(x_lst,true);///< distributed polynom factorizing X(i)X(j)*(...)

  cout<<setw(30)<<std::internal<<"+++INTEGRAL PARAMETERS+++"<<endl;
  cout<<setw(30)<<std::left<<"** Number of loops L=   "<<std::internal<<l<<endl;
  cout<<setw(30)<<std::left<<"** Summ of powers N=   "<<N<<endl;
  cout<<setw(35)<<std::left<<"** Dimension D=   "<<D.subs(D_subs)<<endl;
  cout<<setw(35)<<std::left<<"** U =   "<<U<<endl;
  cout<<setw(35)<<std::left<<"** F =   "<<F<<endl;

  // assuming 1/(U^a*F^b)
  ex U_pow = (-N+(l+1)*D/2).subs(D_subs);
  cout<<"U_pow"<<U_pow<<endl;
  ex F_pow = (N-l*D/2).subs(D_subs);
  lst gamma_lst; //Gamma with poles left of contour
  lst w_lst;
  //  ex gamma_right; // and right of contour
  exmap x_power_map; // msp of X(j) powers
  ex coeff = 1;                             // numerical coeeficient independent of X(j)
  ex out_ex = 1;///pow(2*Pi*I,U.nops()-1)/tgamma(al_pow)/pow(2*Pi*I,F.nops()-1)/tgamma(al_pow); //need review
  //working with F-term \Gamma(\nu-L*D/2) contracted
  ex w_sum = 0;  //F-term generates only integrations in W
  for(const_iterator it = F.begin();it!=F.end();++it)
    {
      int w_index = distance(F.begin(),it);
      ex x_power;
      if(F.end()==boost::next(it)) 
        {
          //          out_ex*=pow(*it,-al_pow-w_sum)*tgamma(al_pow+w_sum);
          x_power = -F_pow-w_sum;
          gamma_lst.append(F_pow+w_sum);
          cout<<"end achieved"<<endl;
        }
      else
        {
          string str = "w_"+boost::lexical_cast<string>(w_index);
          //symbol w(str);
          w_lst.append(get_symbol(str));
          //  out_ex*=pow(*it,w)*tgamma(-w);
          x_power = get_symbol(str);
          gamma_lst.append(-get_symbol(str));
          w_sum+=get_symbol(str);
          cout<<"ok run"<<endl;
        }
      // filling map of X(j) powers
      ex tmp_expr = (*it);
      for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
        {
          if(tmp_expr.has(*it1))
            {
              x_power_map[*it1]+=(tmp_expr.degree(*it1)*x_power);
              tmp_expr = tmp_expr.subs((*it1)==1);
            }
        }
      coeff *=pow(tmp_expr,x_power);
    }

  //working with U-term, only if (U_pow>0)


  if(U_pow.subs(get_symbol("eps")==0) > 0)
    {
      cout<<">>> Working with MB for U term, U_power = "<< U_pow<<endl;  
      ex v_sum = 0;  //F-term generates only integrations in V
      for(const_iterator it = U.begin();it!=U.end();++it)
        {
          int v_index = distance(U.begin(),it);
          ex x_power;
          if(U.end()==boost::next(it)) 
            {
              //          out_ex*=pow(*it,-al_pow-w_sum)*tgamma(al_pow+w_sum);
              x_power = -U_pow - v_sum;
              gamma_lst.append(U_pow+v_sum);
              cout<<"end achieved"<<endl;
            }
          else
            {
              string str = "v_"+boost::lexical_cast<string>(v_index);
              //symbol v(str);
              w_lst.append(get_symbol(str));
              //  out_ex*=pow(*it,w)*tgamma(-w);
              x_power = get_symbol(str);
              gamma_lst.append(-get_symbol(str));
              w_sum+=get_symbol(str);
              cout<<"ok run"<<endl;
            }
          // filling map of X(j) powers
          ex tmp_expr = (*it);
          for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
            {
              if(tmp_expr.has(*it1))
                {
                  x_power_map[*it1]+=(tmp_expr.degree(*it1)*x_power);
                  tmp_expr = tmp_expr.subs((*it1)==1);
                }
            }
          coeff *=pow(tmp_expr,x_power);
        }
    }


//cout<<x_power_map<<endl;
cout<<"Gammas after MB: "<<endl<<gamma_lst<<endl;
 cout<<"X powers list:"<<endl<<x_power_map<<endl;
// applying X integration
 ex gamma_den = 0; // gamma in denominator
 for(exmap::const_iterator mi = x_power_map.begin();mi!=x_power_map.end();++mi)
   {
     cout<<(*mi).first<<" "<<(*mi).second<<endl;
     gamma_lst.append((*mi).second+1);
     gamma_den+=(*mi).second+1;
   }
 cout<<"New gamma list:"<<endl<<gamma_lst<<endl;
 lst lp_lst(w_lst);
 lp_lst.append(get_symbol("eps"));

 gamma_lp(gamma_lst,lp_lst);
 //cout<<hyper_cube(gamma_lst,lp_lst).rhs();
 exmap w0_map(start_point(gamma_lst,lp_lst));
 cout<<"gamma den:"<<gamma_lst.subs(w0_map)<<endl;
 hit_pole(gamma_lst,lp_lst,w0_map);
return out_ex*coeff;
}

*/

typedef std::list<MBintegral> MBlst;

MBlst MBcontinue(MBintegral rootint,double eps0 = 0)
{
try
  {
  rootint.barnes1();
  rootint.barnes2();

 

  MBlst O(1,rootint);
  MBlst C;
  while(O.size()>0)
    {
      MBlst R;
      for(MBlst::iterator it = O.begin();it!=O.end();++it)
        {
          cout<<std::setw(15+it->get_level())<<std::right<<"shifted on "<<it->get_level()<<endl;
          C.push_back(*it);//need review, multiple entries C=C U I
          MBintegral::pole_iterator pit,pit_end;
          ex eps_i = get_symbol("eps");
          //          cout<<"after barness lemas "<<it->get_eps()<<endl;
          eps_i = eps_i.subs(it->get_eps());

          cout<<"eps_i = "<<eps_i<<endl;
          for(tie(pit,pit_end) = it->gamma_args();pit!=pit_end;++pit) // loop over gamma arguments
            {
              cout<<"F(eps_i) "<<pit->subs(it->get_w()).subs(it->get_eps())<<"F(eps=0) "<<pit->subs(it->get_w()).subs(get_symbol("eps")==eps0)<<"   min  "<<std::min(pit->subs(it->get_w()).subs(it->get_eps()),pit->subs(it->get_w()).subs(get_symbol("eps")==eps0))<<endl;
             
              numeric F_eps0 = ex_to<numeric> ( pit->subs( it->get_w() ).subs( get_symbol("eps") == eps0) );
              numeric F_epsi = ex_to<numeric> ( pit->subs( it->get_w() ).subs( it->get_eps() ) );
              if(F_eps0==F_epsi) 
                cout<<"FFFF"<<std::min(F_eps0,F_epsi)<<endl;

              for(int n = 0;n>std::min(F_eps0,F_epsi);n--)
                {
                  cout<<pit->subs(it->get_w()) <<endl;
                  // test on epsilon existance
                  if(pit->subs(it->get_w()).has(get_symbol("eps")))
                    {
                      ex eps_prime = lsolve(pit->subs(it->get_w()) ==n,get_symbol("eps") );
                      cout<<"solve"<<endl;
                      cout<<"F= "<<*pit<<endl;
                      cout<<"eps_i: "<<lsolve(pit->subs(it->get_w()) ==n,get_symbol("eps") )<<endl;
                      cout<<" Poles of Gamma on eps_i: "<<it->get_pole_lst().subs(it->get_w()).subs(get_symbol("eps")==eps_prime)<<endl;
                      lst w_in_F  = has_w(*pit,it->get_w_lst());
                      if(w_in_F.nops()>0)
                        {
                          cout<<lsolve(*pit==n,w_in_F.op(0))<<endl;
                          cout<<"sign(z) = "<<csgn(pit->coeff(w_in_F.op(0)))<<"     sign(F_i-F_0) = "<<csgn(F_epsi-F_eps0)<<endl;
                      
                          //MBintegral newi(it->get_w_lst(),it->get_pole_lst(),it->get_expr(),it->get_w(),it->get_eps());
                          //      cout<<"NEW INT: "<< newi.get_expr()<<endl;
                          cout<<"debug fepsi"
                              <<"1: "<<lsolve(*pit==n,w_in_F.op(0))<<endl
                              <<"2: "
                              <<endl;
                          MBintegral res_int = it->res(w_in_F.op(0)==lsolve(*pit==n,w_in_F.op(0)),*pit,get_symbol("eps")==eps_prime);
                          res_int.set_level(1+it->get_level());
                          cout<<"Storing RES_INT with eps_prime = " <<eps_prime<<"  "<<res_int.get_eps().rhs()<<endl;
                          res_int*=(2*Pi*I*csgn(pit->coeff(w_in_F.op(0)))*csgn(F_epsi-F_eps0));
                          cout<<"RES EXPR:  "<<res_int.get_expr()<<endl;
                          res_int.barnes1();
                          res_int.barnes2();
                          R.push_back(res_int);
                        }
                    }
                  else BOOST_ASSERT_MSG("EEEEEERRRRRRRROOOORR: no W dependence in pole",false);
                         //cout<<endl<<endl<<"EEEEEERRRRRRRROOOORR: no W dependence in pole"<<endl<<endl;
                }

            }
        }
      O = R;
    }
  cout<<"Continue get "<<C.size()<<" integrals"<<endl;
  cout<< endl<<" Next step?  [Y/n]: ";
  
  char in_ch;
  std::cin>>in_ch;
  if(in_ch=='n')  exit(0);//assert(false);

  for(MBlst::iterator it = C.begin();it!= C.end();++it)
    {
      //cout<<(it->get_pole_lst().subs(it->get_w())).subs(get_symbol("eps")==0)<<endl;
      // cout<<endl <<"series:  " << it->get_gamma_expr().series(get_symbol("eps")==0,0)<<endl<<endl;
    }
  return C;
  }catch(std::exception &p)
  {
    throw std::logic_error(std::string("In function \"MBcontinue\":\n |___> ")+p.what());
  }
}

class RoMB_np
{
public:
  /**
   *
   *  loop momentums,propagator expressions,
   *  invariants substitutions,propagator powers,number of loops
   *
   */
  RoMB_np(lst k_lst,lst p_lst,lst subs_lst,lst nu,numeric l)
  {
    try
      {
        UFXmap inUFmap = UF(k_lst,p_lst,subs_lst);
        //    lst U_lst = ex_to<lst>(fusion::at_key<UFX::U>(inUFmap));
        //  cout<<"U: "<< fusion::at_key<UFX::U>(inUFmap)<<endl;//<<" U_lst: "<<U_lst<<endl;


        ex U_sum =fusion::at_key<UFX::U>(inUFmap).collect(  fusion::at_key<UFX::xlst>(inUFmap),true); 
        cout<<"U_sum: "<<U_sum<<endl;
        if(is_a<add>(U_sum))
          for(const_iterator it = U_sum.begin();it!=U_sum.end();++it)
            {
              MBintegral Uint(fusion::make_map<UFX::U,UFX::F,UFX::xlst>(*it,fusion::at_key<UFX::F>(inUFmap),fusion::at_key<UFX::xlst>(inUFmap)),nu,l);
              //     cout<<"ui9nt eps : "<<Uint.get_eps().lhs()<<endl;
              Uint.new_point();

              MBcontinue(Uint);
            }
        else cout<<" ONE INT"<<endl;


      }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"RoMB\":\n |___> ")+p.what());
      }
  }
};




ex expand_and_integrate(MBintegral& int_in, lst num_subs, int expansion_order = 1) // up to O(eps^1) 
{
  try
    {
      ex out_ex;

      //            cout<<(it->get_pole_lst().subs(it->get_w())).subs(get_symbol("eps")==0)<<endl;
      //         cout<<endl <<"series:  " << it->get_gamma_expr().series(get_symbol("eps")==0,expansion_order)<<endl<<endl;
      lst w_lst = int_in.get_w_lst();
      if(int_in.get_w_lst().nops()>0)// expanding and integrating
        {
          int_in.barnes1();
          int_in.barnes2();
          out_ex = series_to_poly( int_in.get_expr().series(get_symbol("eps"),expansion_order) ).subs(num_subs);
          // loop over W_i, converting integration contour
          for(lst::const_iterator wit = w_lst.begin(); wit != w_lst.end(); ++wit)
            {
              ex c_i = wit->subs(int_in.get_w());
              out_ex = (I*out_ex.subs((*wit)==c_i - I*log( (*wit)/( 1 - (*wit) ) ) ) ) / (*wit)/(1- (*wit));
            }
          cout<<"current : "<<out_ex<<endl;
          ex wo_eps_part = out_ex;
          ex vegas_ex = 0;
          for(int i = out_ex.ldegree( get_symbol("eps") ); i <= out_ex.degree( get_symbol("eps") ); i++)
            {
              cout<<"Ord( "<<i<<" ) coeff : "<< out_ex.coeff(get_symbol("eps"),i)<<endl;
              ex int_expr =  out_ex.coeff(get_symbol("eps"),i);
              RoMB::FUNCP_CUBA2 fp_real,fp_imag;
              RoMB::compile_ex_real(lst(int_expr),int_in.get_w_lst(), fp_real);
              RoMB::compile_ex_imag(lst(int_expr),int_in.get_w_lst(), fp_imag);

              // ----------------------------------- Vegas integration-------------------------
               int  NDIM  = int_in.get_w_lst().nops();
              //#define NCOMP 1
#define USERDATA NULL
#define EPSREL 1e-4
#define EPSABS 1e-12
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 1000000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
//SUAVE
                #define NNEW 1000
                #define FLATNESS 25.

                //#define KEY1 47
                //#define KEY2 1
                //#define KEY3 1
                //#define MAXPASS 5
                //#define BORDER 0.
                //#define MAXCHISQ 10.
                //#define MINDEVIATION .25
                //#define NGIVEN 0
                //#define LDXGIVEN NDIM
                //#define NEXTRA 0

                #define KEY 0
              const int NCOMP = 1;
              int comp, nregions, neval, fail;
              double integral[NCOMP],integral_real[NCOMP],integral_imag[NCOMP], error[NCOMP], prob[NCOMP];
                   
              printf("-------------------- Vegas test --------------------\n");
	      /*              
			      Vegas(NDIM, NCOMP, fp,
			      EPSREL, EPSABS, VERBOSE, 
			      MINEVAL, MAXEVAL,
			      NSTART, NINCREASE,
			      &neval, &fail, integral, error, prob);
	*/     
	     
	      Vegas(NDIM, NCOMP, fp_real, USERDATA,
		    EPSREL, EPSABS, VERBOSE, SEED,
                    MINEVAL, MAXEVAL, 
		    NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE,
		    &neval, &fail, integral_real, error, prob);


          
             // printf("VEGAS RESULT:\tneval %d\tfail %d\n",
              
                for( comp = 0; comp < NCOMP; ++comp )
                printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                       integral_real[comp], error[comp], prob[comp]);
                       
                       	      Vegas(NDIM, NCOMP, fp_imag, USERDATA,
		    EPSREL, EPSABS, VERBOSE, SEED,
                    MINEVAL, MAXEVAL, 
		    NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE,
		    &neval, &fail, integral_imag, error, prob);
		                  for( comp = 0; comp < NCOMP; ++comp )
                printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                       integral_imag[comp], error[comp], prob[comp]);


            
              /*
               printf("\n-------------------- Suave test real--------------------\n");
               
                 Suave(NDIM, NCOMP, fp_real, USERDATA,
                     EPSREL, EPSABS, VERBOSE | LAST, SEED,
                         MINEVAL, MAXEVAL, NNEW, FLATNESS,
                             &nregions, &neval, &fail, integral_real, error, prob);
                             
                               printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
                                   nregions, neval, fail);
                                   
              for( comp = 0; comp < NCOMP; ++comp )
                printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                       integral_real[comp], error[comp], prob[comp]);


               printf("\n-------------------- Suave test imag--------------------\n");
               
                 Suave(NDIM, NCOMP, fp_imag, USERDATA,
                     EPSREL, EPSABS, VERBOSE | LAST, SEED,
                         MINEVAL, MAXEVAL, NNEW, FLATNESS,
                             &nregions, &neval, &fail, integral_imag, error, prob);
                             
                               printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
                                   nregions, neval, fail);*/
                                   
             
              // ----------------------------------- Vegas integration-------------------------              
              vegas_ex+=pow(get_symbol("eps"),i)*(integral_real[0]+I*integral_imag[0]);
            }
          out_ex = vegas_ex;
        }
      else // expanding only
        {
          out_ex = series_to_poly( int_in.get_expr().series(get_symbol("eps"),expansion_order) ).subs(num_subs);
        }
      cout<<endl<<"gp "<<out_ex<<endl;
      return out_ex;
    }catch(std::exception &p)
    {
      throw std::logic_error(std::string("In function \"RoMB\":\n |___> ")+p.what());
    }
}

class RoMB_planar
{
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
  RoMB_planar(
              lst k_lst,
              lst p_lst,
              lst subs_lst,
              lst nu,
              numeric l
              )
  {
    try
      {
        UFXmap inUFmap = UF(k_lst,p_lst,subs_lst);
        //    lst U_lst = ex_to<lst>(fusion::at_key<UFX::U>(inUFmap));
        //  cout<<"U: "<< fusion::at_key<UFX::U>(inUFmap)<<endl;//<<" U_lst: "<<U_lst<<endl;


        ex U_sum =fusion::at_key<UFX::U>(inUFmap).collect(  fusion::at_key<UFX::xlst>(inUFmap),true); 
        cout<<"U_sum: "<<U_sum<<endl;
        MBintegral Uint(
                        fusion::make_map<UFX::F,UFX::xlst>(fusion::at_key<UFX::F>(inUFmap),
                                                           fusion::at_key<UFX::xlst>(inUFmap)
                                                           ),nu,l);
        //     cout<<"ui9nt eps : "<<Uint.get_eps().lhs()<<endl;
        Uint.new_point();
        
        MBlst int_lst = MBcontinue(Uint);
        /*
        ex int_expr_out = 0;
        for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
          int_expr_out+=expand_and_integrate(*it,1);
        
        cout<<" RESULT : "<<endl
            <<"               = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
        */
      }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"RoMB\":\n |___> ")+p.what());
      }
  }
};


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
  RoMB_loop_by_loop(
              lst k_lst,
              lst p_lst,
              lst subs_lst,
              lst nu)
  {
    try
      {
        MBintegral MBlbl_int(lst(),lst(),1);
        exset input_prop_set( p_lst.begin(),p_lst.end());
        lst new_prop_lst;
        exmap prop_pow_map;
        for(lst::const_iterator Pit = p_lst.begin(); Pit != p_lst.end(); ++Pit)
          {
            prop_pow_map[*Pit] = nu.op(std::distance(p_lst.begin(),Pit));
          }
        unsigned int displacement_x = 0;
        unsigned int displacement_w = 0;
        for(lst::const_iterator kit = k_lst.begin(); kit != k_lst.end(); ++kit)
          {
            cout<<"PROP_POW_MAP "<<prop_pow_map<<endl;
            lst tmp_p_lst=set2lst(input_prop_set);
            for(lst::const_iterator Pit = new_prop_lst.begin(); Pit != new_prop_lst.end(); ++Pit)
              tmp_p_lst.append(*Pit);
            new_prop_lst.remove_all();
            lst P_with_k_lst;
            for(lst::const_iterator Pit = tmp_p_lst.begin(); Pit != tmp_p_lst.end(); ++Pit)
              if(Pit->has(*kit))
                {
                  P_with_k_lst.append(*Pit);
                  input_prop_set.erase(*Pit);
                }
            cout<< "Set wo k_i "<<input_prop_set<<endl;
            cout<<" PWKlst "<<P_with_k_lst<<endl;

            // uf and then MB represenatation construction
            // subs only in F for last momentum
            UFXmap inUFmap;
           if(boost::next(kit) == k_lst.end())
            inUFmap = UF(lst(*kit),P_with_k_lst,subs_lst,displacement_x);
           else
            inUFmap = UF(lst(*kit),P_with_k_lst,lst(),displacement_x); // no substitution!!!
            displacement_x +=fusion::at_key<UFX::xlst>(inUFmap).nops(); 

            lst nu_into;
            for(lst::const_iterator nuit = P_with_k_lst.begin(); nuit != P_with_k_lst.end(); ++nuit )
              nu_into.append(prop_pow_map[*nuit]);
            cout<<" Powers list before input: "<<nu_into<<endl;
            MBintegral Uint(
                            fusion::make_map<UFX::F,UFX::xlst>(fusion::at_key<UFX::F>(inUFmap),
                                                               fusion::at_key<UFX::xlst>(inUFmap)
                                                               ),nu_into,1,displacement_w);
            displacement_w+=Uint.get_w_lst().nops();
            cout<<"ui9nt eps : "<<Uint.get_expr()<<endl;

            //expression to mul root integral
            ex expr_k_to_subs_1= Uint.get_expr();

            ex mom_find = Uint.get_expr();
            if(is_a<mul>(mom_find))
              {
                exset found_prop;
                mom_find.find(pow(wild(1),wild(2)),found_prop);
                cout<<" is a mul  "<<found_prop<<endl;
                lst prop_pow_lst;
                if(boost::next(kit)!=k_lst.end())
                  {
                    ex mom_to_find = *(boost::next(kit));
                    BOOST_FOREACH(ex propex,found_prop)
                      {
                        cout<<"kit "<<mom_to_find<<" propex "<<propex<<endl;
                        if(propex.has(mom_to_find))
                          {
                            cout<<"before subs kex : "<<expr_k_to_subs_1<<endl;
                            expr_k_to_subs_1 = expr_k_to_subs_1.subs(propex == 1);
                              cout<<"after subs kex : "<<expr_k_to_subs_1<<endl;
                            new_prop_lst.append(ex_to<power>(propex).op(0));
                            prop_pow_lst.append(ex_to<power>(propex).op(1));
                            prop_pow_map[ex_to<power>(propex).op(0)] = (-1)*ex_to<power>(propex).op(1);
                          }
                      }
                    cout<<"needed props "<<prop_pow_lst<<endl;
                  }
              }
            
            MBlbl_int*=expr_k_to_subs_1;
            MBlbl_int.add_w_lst(Uint.get_w_lst());
            MBlbl_int.add_poles_lst(Uint.get_pole_lst());
            
          }
        cout<<"Constructed integral with:"<<endl;
        cout<<"Poles: "<<MBlbl_int.get_pole_lst()<<endl;
        cout<<"W's : "<<MBlbl_int.get_w_lst()<<endl;
        cout<<"Expr : "<<MBlbl_int.get_expr()<<endl;
        
        cout<< endl<<" Ready for MBcontinue?  [Y/n]: ";
        char in_ch;
        std::cin>>in_ch;
        if(in_ch=='n')  exit(0);//assert(false);

        //                assert(false);
        MBlbl_int.new_point();
	//        MBlst int_lst = MBcontinue(MBlbl_int);
	int_lst = MBcontinue(MBlbl_int);
        /*
	ex int_expr_out = 0;
        for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
         int_expr_out+=expand_and_integrate(*it,1);
        
        cout<<" FRESULT anl : "<<"          = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
        cout<<" FRESULT num: "<<"          = "<<evalf(int_expr_out.expand().collect(get_symbol( "eps" )))<<endl;
	*/
        //    lst U_lst = ex_to<lst>(fusion::at_key<UFX::U>(inUFmap));
        //  cout<<"U: "<< fusion::at_key<UFX::U>(inUFmap)<<endl;//<<" U_lst: "<<U_lst<<endl;


        //        ex U_sum =fusion::at_key<UFX::U>(inUFmap).collect(  fusion::at_key<UFX::xlst>(inUFmap),true); 
        // cout<<"U_sum: "<<U_sum<<endl;
        // MBintegral Uint(
        //                fusion::make_map<UFX::F,UFX::xlst>(fusion::at_key<UFX::F>(inUFmap),
        //                                                   fusion::at_key<UFX::xlst>(inUFmap)
        //                                                  ),nu,l);
        //     cout<<"ui9nt eps : "<<Uint.get_eps().lhs()<<endl;
        // Uint.new_point();
        
        // MBlst int_lst = MBcontinue(Uint);
        //ex int_expr_out = 0;
        //for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
        // int_expr_out+=expand_and_integrate(*it,1);
        
        //cout<<" RESULT : "<<endl
        //    <<"               = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
      }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"RoMB\":\n |___> ")+p.what());
      }
  }
  void integrate(lst number_subs_list, int exp_order = 1)
  {
    ex int_expr_out = 0;
    for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
      int_expr_out+=expand_and_integrate(*it, number_subs_list, exp_order);
    cout<<" FRESULT for parameters: "<<number_subs_list<<endl<<endl;
    cout<<" FRESULT anl : "<<"          = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
    cout<<" FRESULT num: "<<"          = "<<evalf(int_expr_out.expand().collect(get_symbol( "eps" )))<<endl;
  }
};


int main()
{
try
  {
    symbol k("k"),q("q"),p("p"),p1("p1"),p2("p2"),p3("p3"),ms("ms"),l("l"),s("s");
  symbol l1("l1"),l2("l2"),l3("l3"),l4("l4");
  // oneloop box
  //      UFXmap l45 = UF(lst(k),lst(pow(k,2),pow(k+p1,2),pow(k+p1+p2,2),pow(k+p1+p2+p3,2)),lst(pow(p1,2)==0,pow(p2,2)==0));
  // MBintegral root_int(l45,lst(1,1,1,1),1);

   //two loop box bubble
  // UFXmap l45 = UF(lst(k,l),lst(pow(k,2),pow(k+p1,2),pow(k+p1+p2,2),pow(l+p1+p2,2),pow(l+p1+p2+p3,2),pow(l,2),pow(k-l,2)),lst(pow(p1,2)==0,pow(p2,2)==0,pow(p3,2)==0));
  //MBintegral root_int(l45,lst(1,1,1,1,1,1,1),2);

  // B0
  //    UFXmap l45 = UF(lst(k),lst(ms-pow(k,2),ms-pow(-k,2)),lst(ms==1));
  // MBintegral root_int(l45,lst(1,1),1);

  // 2 loop sunrise
  //UFXmap l45 = UF(lst(k,q),lst(ms-pow(k,2),ms-pow(-q-k,2),ms-pow(q,2)),lst(ms==1));
  //MBintegral root_int(l45,lst(1,1,1),2);


  //RoMB_planar box2loop(lst(k,l),lst(pow(k,2),pow(k+p1,2),pow(k+p1+p2,2),pow(l+p1+p2,2),pow(l+p1+p2+p3,2),pow(l,2),pow(k-l,2)),lst(pow(p1,2)==0,pow(p2,2)==0,pow(p3,2)==0),lst(1,1,1,1,1,1,1),2);

  //  RoMB_planar  box1loop(lst(k),lst(pow(k,2),pow(k+p1,2)-ms,pow(k+p1+p2,2),pow(k+p1+p2+p3,2)),lst(pow(p1,2)==0,pow(p2,2)==0,pow(p3,2)==0,p1==0,p2==0,p3==0,ms==1),lst(1,1,1,1),1);

  //  RoMB_planar B0_1loop(lst(k),lst(pow(k,2)-ms,pow(p+k,2)-ms),lst(ms==0,pow(p,2)==1),lst(1,1),1);

  //  RoMB_planar C0_1loop(lst(k),lst(pow(k,2)-ms,pow(p1+k,2)-ms,pow(p1+p2+k,2)),lst(ms==1,pow(p1,2)==0,pow(p2,2)==0,p1*p2==50),lst(1,1,1),1);
//cout<<" new point "<<endl<<root_int.new_point()<<endl;
// cout<<" saved point "<<endl<<root_int.get_point()<<endl;
//  MBcontinue(root_int);
  //cout<<MB_lst(l45,lst(1,1,1,1),1).expand()<<endl;


  // RoMB_loop_by_loop box2loop(lst(k,l),lst(pow(k,2),pow(k+p1,2),pow(k+p1+p2,2),pow(l+p1+p2,2),pow(l+p1+p2+p3,2),pow(l,2),pow(k-l,2)),lst(pow(p1,2)==0,pow(p2,2)==0,pow(p3,2)==0),lst(1,1,1,1,1,1,1));
  //      RoMB_loop_by_loop t2(lst(k,l), lst(pow(k,2),pow(p+k,2),pow(p+k+l,2),pow(l,2),pow(k+l,2)),lst(pow(p,2)==1),lst(1,1,1,1,1));


  // works!!!
//      RoMB_loop_by_loop sunset(lst(k,l), lst(pow(k,2)-1,pow(p-k-l,2)-4,pow(l,2)-5),lst(pow(p,2)==s),lst(1,1,1));
//      RoMB_loop_by_loop sunset(lst(k,l), lst(pow(k,2),pow(p-k-l,2),pow(l,2)),lst(pow(p,2)==s),lst(1,1,1));
//      sunset.integrate(lst(s==1));
      
      RoMB_loop_by_loop t2loop(lst(k,l), lst(pow(k,2),pow(p+k,2),pow(p+k+l,2),pow(k+l,2),pow(l,2)),lst(pow(p,2)==s),lst(1,1,1,1,1));
      t2loop.integrate(lst(s==1));
  
/*        RoMB_loop_by_loop bubble_five_loop(lst(k,l1,l2,l3,l4), 
        lst(pow(k,2)-ms,pow(l1,2)-ms,pow(l2,2)-ms,pow(l3,2)-ms,pow(l4,2)-ms,pow(k+l1,2)-ms,pow(k+l1+l2,2)-ms,pow(k+l1+l2+l3,2)-ms,pow(k+l1+l2+l3+l4,2)-ms,pow(k+l1+l2+l3,2)-ms,pow(k+l1+l2,2)-ms,pow(k+l1,2)-ms),
        lst(ms==1),
        lst(1,1,1,1,1,1,1,1,1,1,1,1));
*/

  // works!!!
//             RoMB_loop_by_loop B0_1loop_lbl(lst(k),lst(pow(k,2)-2-ms,pow(p+k,2)-ms),lst(ms==0,pow(p,2)==1),lst(2,1));

//    RoMB_loop_by_loop B0_1loop_lbl(lst(k),lst(pow(k,2)-ms,pow(p+k,2)-ms),lst(pow(p,2)==s,ms==0),lst(1,1));
//	     B0_1loop_lbl.integrate(lst(s==-1));

  //MB works???
//    RoMB_loop_by_loop C0_1loop_lbl(lst(k),lst(pow(k,2),pow(k+p1,2)-ms,pow(k-p2,2)-ms),lst(ms==1,pow(p1,2)==1,pow(p2,2)==1,p1*p2==-50-1),lst(1,1,1));

  //MB works???
  //RoMB_loop_by_loop box1looplbl(lst(k),lst(pow(k,2)-ms,pow(k+p1,2)-ms,pow(k+p1+p2,2)-ms,pow(k+p1+p2+p3,2)-ms),lst(pow(p1,2)==0,pow(p2,2)==0,pow(p3,2)==0,p1==0,p2==0,p3==0,ms==1),lst(1,1,1,1));
  }
  catch(std::exception &p)
    {
      std::cerr<<"******************************************************************"<<endl;
      std::cerr<<"   >>>ERROR:  "<<p.what()<<endl;
      std::cerr<<"******************************************************************"<<endl;
      return 1;
    }
  return 0;
}





