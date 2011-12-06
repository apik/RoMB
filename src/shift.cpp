#include "shift.h"
#include <iostream>
#include <cmath>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include "compat.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_psi.h>
#include <boost/math/special_functions/digamma.hpp>
#include <limits>

using namespace std;


template<class UserFunc> 
struct  GSLMultiMinFunctionAdapter 
{
  
  static double F( const gsl_vector * x, void * p) 
  { 
    
    UserFunc * function = reinterpret_cast< UserFunc *> (p);
    // get pointer to data from gsl_vector
    return (*function)( x->data ); 
  }
  
  
  static void Df(  const gsl_vector * x, void * p,  gsl_vector * g) 
  { 
    
    UserFunc * function = reinterpret_cast< UserFunc *> (p); 
    (*function).Gradient( x->data, g->data );
    
  }
  
  static void Fdf( const gsl_vector * x, void * p, double *f, gsl_vector * g ) 
  { 
    
    UserFunc * function = reinterpret_cast< UserFunc *> (p);      
    *f = (*function) ( x->data ); 
    (*function).Gradient( x->data,g->data ); 
    
  }
};
/*****************************************
 *         WARNING GLOBAL VARIABLES
 ****************************************/
  ex GLOBAL_F_T_M;
  lst GLOBAL_WL;
  exset GLOBAL_INEQ;
  exmap GLOBAL_F_GRAD;


namespace evalfPsi
{
DECLARE_FUNCTION_1P(digamma)
DECLARE_FUNCTION_2P(polygamma)

  //////////
  // Psi-function (aka digamma-function)
  //////////
  
  static ex digamma_evalf(const ex & x)
  {
    if (is_exactly_a<numeric>(x)) {
    
      return boost::math::digamma(ex_to<numeric>(x).to_double());
    
    }
    
    return digamma(x).hold();
  }
  REGISTER_FUNCTION(digamma,evalf_func(digamma_evalf).
                    latex_name("\\digamma"));
  
  //////////
  // Psi-functions (aka polygamma-functions)  psi(0,x)==psi(x)
  //////////
  
  static ex polygamma_evalf(const ex & n, const ex & x)
  {
    if (is_exactly_a<numeric>(n) && is_exactly_a<numeric>(x)) 
      {
        // if(ex_to<numeric>(n).is_nonneg_integer())
        //  {
        // cout<<"Eval for psi("<<n<<","<<x<<")"<<endl;
            return gsl_sf_psi_n(int(ex_to<numeric>(n).to_double()),ex_to<numeric>(x).to_double());
            // }
            /* else 
          {
            cout<<n<<endl;
            throw std::logic_error(std::string("Non int power of Psi functionNot a 32-bit integer:"));
            }   */    
      }
    
    return polygamma(n,x).hold();
  }
 REGISTER_FUNCTION(polygamma,evalf_func(polygamma_evalf).
                    latex_name("\\polygamma"));
 
 } // namespace evalfPsi

struct minimize_me
{
  double operator()( const double *  t)
  {
    exmap subs_f_val;
    size_t j = 0;
    for(lst::const_iterator wii = GLOBAL_WL.begin(); wii != GLOBAL_WL.end(); ++wii, j++)
      subs_f_val[*wii] = t[j];

    BOOST_FOREACH(ex ineee, GLOBAL_INEQ)
      {
        if(ineee.subs(subs_f_val) < 0) return numeric_limits<double>::infinity();
      }

    //  cout<<ex_to<numeric>(evalf(GLOBAL_F_T_M.subs(subs_f_val))).to_double()<<endl;
    return ex_to<numeric>(
                          evalf(GLOBAL_F_T_M.subs(subs_f_val).subs(psi(wild(1),wild(2))==evalfPsi::polygamma(wild(1),wild(2))).subs(psi(wild(3))==evalfPsi::digamma(wild(3))))).to_double();
  }
  void Gradient( const double *   t, double * g)
  {
    // cout<<"In gradient!!!"<<endl;
    exmap subs_f_val;
    size_t j = 0;
    for(lst::const_iterator wii = GLOBAL_WL.begin(); wii != GLOBAL_WL.end(); ++wii, j++)
      subs_f_val[*wii] = t[j];
    j = 0;
    for(lst::const_iterator wii = GLOBAL_WL.begin(); wii != GLOBAL_WL.end(); ++wii, j++)
      {
        //  cout<<"Jig"<<evalf(GLOBAL_F_GRAD[*wii].subs(subs_f_val).subs(psi(wild(1),wild(2))==evalfPsi::polygamma(wild(1),wild(2))).subs(psi(wild(3))==evalfPsi::digamma(wild(3))))<<endl;       
   
        g[j] = ex_to<numeric>(evalf(GLOBAL_F_GRAD[*wii].subs(subs_f_val).subs(psi(wild(1),wild(2))==evalfPsi::polygamma(wild(1),wild(2))).subs(psi(wild(3))==evalfPsi::digamma(wild(3))))).to_double();
        
      }
  }
};




exmap shift_contours(const exset& fl,const lst& pl, exmap sm)
{
try
  {
  exmap sc;
  ex min_func = 0;
  lst wl;
  exset ineql;
  std::vector<double> inb; 
  // fill list of unknown variables
  typedef std::pair<ex,ex> wpair; 
  BOOST_FOREACH(wpair ep, sm)
    {
      if(pl.has(ep.first))wl.append(ep.first);
    }
  cout<<endl<<wl<<"Original values : "<<sm<<endl;
  
  int n = wl.nops();
  int p = 0;
  for(lst::const_iterator pi = pl.begin(); pi != pl.end(); ++pi)
    {
      double Fb = ex_to<numeric>(pi->subs(sm)).to_double(); 
      if(Fb < 0)
        {
          cout<<Fb<<"  ceil: "<<ceil(Fb)<<"  floor: "<<floor(Fb)<<endl;
          min_func += pow(*pi - numeric(ceil(Fb) + floor(Fb),2),2);
          ineql.insert(*pi - floor(Fb));
          //          inb.push_back(-floor(Fb));
          // cout<<*ineql.rbegin()<<" + ("<<inb.back()<<") >= 0"<<endl;
          ineql.insert(ceil(Fb) - (*pi));
          //inb.push_back(ceil(Fb));
          //cout<<*ineql.rbegin()<<" + ("<<inb.back()<<") >=0"<<endl;
        }
      else
        {
          cout<<Fb<<"  ceil: "<<ceil(Fb)<<"  floor: "<<floor(Fb)<<endl;
          min_func += pow(*pi - ceil(Fb),2);
          ineql.insert(*pi - floor(Fb));
          //inb.push_back(-floor(Fb)); // or 0???
          //cout<<*ineql.rbegin()<<" + ("<<inb.back()<<") >= 0"<<endl;
        }
    }


  ex f_min = 0;
  exmap f_grad;
  BOOST_FOREACH(ex f_term, fl)
    {
      f_min+=pow(f_term,2);
    }
    min_func = min_func.expand();
  cout<<"INEQL:"<<ineql<<endl;
  cout<<"FUNCS: "<< f_min<<endl;
  for(lst::const_iterator wli = wl.begin(); wli != wl.end(); ++wli)
    {
      f_grad[*wli] = 2*f_min.diff(ex_to<symbol>(*wli));
    }
 

  GLOBAL_F_T_M = f_min;
  GLOBAL_WL = wl;
  GLOBAL_INEQ = ineql;
  GLOBAL_F_GRAD = f_grad;
 cout<<"grad: "<<GLOBAL_F_GRAD<<endl;


  size_t iter = 0;
  int status;
  
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  
  /* Position of the minimum (1,2), scale factors 
     10,20, height 30. */
  double par[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 };
     
  gsl_vector *x;
  
  gsl_multimin_function_fdf my_func;
  
  my_func.n = wl.nops();
  my_func.f =  GSLMultiMinFunctionAdapter<minimize_me>::F;
    //reinterpret_cast<double (*)(const gsl_vector*, void*)>(&minimize_me::my_f);

  my_func.df = GSLMultiMinFunctionAdapter<minimize_me>::Df;
  //reinterpret_cast<void (*)(const gsl_vector*, void*, gsl_vector*)>(&minimize_me::my_df);
  my_func.fdf = GSLMultiMinFunctionAdapter<minimize_me>::Fdf;
  //reinterpret_cast<void (*)(const gsl_vector*, void*, double*, gsl_vector*)>(&X::my_fdf);
  my_func.params = par;
   
  /* Starting point, x = (5,7) */
  x = gsl_vector_alloc (wl.nops());
  int h = 0;
  for(lst::const_iterator xi = wl.begin(); xi != wl.end(); ++xi,h++)
    {
      ex w_x = sm[(*xi)];
      // cout<<w_x<<endl;
      double x_in = ex_to<numeric>(w_x).to_double();
      gsl_vector_set (x, h,x_in );
    }
  //  T = gsl_multimin_fdfminimizer_conjugate_fr;
  T = gsl_multimin_fdfminimizer_steepest_descent;
  s = gsl_multimin_fdfminimizer_alloc (T, wl.nops());
  //cout << "Starting minim proc"<<endl;      
  //gsl_multimin_function_fdf my_func;
  //int a = 6;

  // X::get_my_func(3,my_func);
  // my_func.params = par;
  gsl_multimin_fdfminimizer_set (s, &my_func , x, 0.01, 1e-4);
  
  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
     
      if (status)
        break;
     
      status = gsl_multimin_test_gradient (s->gradient, 1e-3);
     
      if (status == GSL_SUCCESS)
        printf ("Minimum found at:\n");
      h = 0;
      printf ("Iter:( %5d )", iter);
      for(lst::const_iterator xi = wl.begin(); xi != wl.end(); ++xi,h++)
        {
          sc[*xi] = gsl_vector_get (s->x, h);
          printf ("%.5f ", gsl_vector_get (s->x, h));
        } 
      printf(" %10.5f\n",s->f);
    }
  while (status == GSL_CONTINUE && iter < 1000);
     
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

  lst newpl;
  lst oldpl;
  
  for(int jk = 0; jk < pl.nops(); jk++)
    {      
      newpl.append(pl.op(jk).subs(sc));
      oldpl.append(pl.op(jk).subs(sm));
    }
  cout<<"OLD_POLES: "<<oldpl<<endl;
  cout<<"NEW_POLES: "<<newpl<<endl;
  for(int jk = 0; jk < pl.nops(); jk++)
    if(csgn(newpl.op(jk)) != csgn(oldpl.op(jk))) throw std::logic_error(std::string("Pole crossed!!!"));
  return sc;
 }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"Shift contours\":\n |___> ")+p.what());
      }
}
