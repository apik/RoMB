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

double npsi(double x_in, int k_in)
{
  double *x = &x_in;
  int *k = &k_in;
static int c__1 = 1;

    /* Initialized data */

    static int nb[6] = { 16,17,17,18,19,20 };
    static struct {
	double e_1[17];
	double fill_2[4];
	double e_3[18];
	double fill_4[3];
	double e_5[18];
	double fill_6[3];
	double e_7[19];
	double fill_8[2];
	double e_9[20];
	double fill_10[1];
	double e_11[21];
	} equiv_31 = { .334838697910949386, -.055187482048730095, 
		.004510190736011502, -3.65705888303721e-4, 2.9434627468223e-5,
		 -2.352776815151e-6, 1.86853176633e-7, -1.4750720184e-8, 
		1.157993337e-9, -9.0439179e-11, 7.029627e-12, -5.43989e-13, 
		4.1925e-14, -3.219e-15, 2.46e-16, -1.9e-17, 1e-18, {0}, 
		-.11259293534547383, .036557001742820941, 
		-.004435942496027282, 4.75475854728926e-4, 
		-4.7471836382632e-5, 4.521815237353e-6, -4.1630007962e-7, 
		3.7338998165e-8, -3.279914474e-9, 2.83211377e-10, 
		-2.4104028e-11, 2.026297e-12, -1.68524e-13, 1.3885e-14, 
		-1.135e-15, 9.2e-17, -7e-18, 1e-18, {0}, .076012604655110384, 
		-.036257186481828739, .005797202338937002, 
		-7.69646513610481e-4, 9.1492082189884e-5, -1.0097131488364e-5,
		 1.055777442831e-6, -1.05929577481e-7, 1.0285494201e-8, 
		-9.7231431e-10, 8.9884635e-11, -8.153171e-12, 7.27572e-13, 
		-6.401e-14, 5.562e-15, -4.78e-16, 4.1e-17, -3e-18, {0}, 
		-.077234724056994793, .047867163451599467, 
		-.009440702186674632, .001489544740103448, 
		-2.0494402334886e-4, 2.5671425065297e-5, -3.001393581584e-6, 
		3.32766437356e-7, -3.5365412111e-8, 3.630622927e-9, 
		-3.62096951e-10, 3.5237509e-11, -3.35744e-12, 3.14068e-13, 
		-2.8908e-14, 2.623e-15, -2.35e-16, 2.1e-17, -2e-18, {0}, 
		.104933034459278632, -.078877901652793557, 
		.018397415112159397, -.003352284159396504, 
		5.22878230918016e-4, -7.317978581474e-5, 9.449729612085e-6, 
		-1.146339856723e-6, 1.32269366108e-7, -1.464666918e-8, 
		1.566940742e-9, -1.62791157e-10, 1.6490345e-11, -1.634028e-12,
		 1.58807e-13, -1.5171e-14, 1.427e-15, -1.32e-16, 1.2e-17, 
		-1e-18, {0}, -.178617622142502753, .155776462200520579, 
		-.041723637673831277, .0085971413032454, -.001496227761073229,
		 2.31089608557137e-4, -3.2632044778436e-5, 4.29609786709e-6, 
		-5.34528790204e-7, 6.3478151644e-8, -7.248699714e-9, 
		8.00521979e-10, -8.5888793e-11, 8.985442e-12, -9.19356e-13, 
		9.2225e-14, -9.09e-15, 8.82e-16, -8.4e-17, 8e-18, -1e-18 };

#define b ((double *)&equiv_31)

    static double c__[42]	/* was [7][6] */ = { .166666666666666667,
	    -.0333333333333333333,.0238095238095238095,-.0333333333333333333,
	    .0757575757575757576,-.253113553113553114,1.16666666666666667,.5,
	    -.166666666666666667,.166666666666666667,-.3,.833333333333333333,
	    -3.29047619047619048,17.5,2.,-1.,1.33333333333333333,-3.,10.,
	    -46.0666666666666667,280.,10.,-7.,12.,-33.,130.,-691.,4760.,60.,
	    -56.,120.,-396.,1820.,-11056.,85680.,420.,-504.,1320.,-5148.,
	    27300.,-187952.,1627920. };
    static double sgn[6] = { -1.,1.,-1.,1.,-1.,1. };
    static double sgf[7] = { 1.,-1.,2.,-6.,24.,-120.,720. };
    static double sgh[6] = { -.5,1.,-3.,12.,-60.,360. };
    static double x0 = 1.46163214496836234;
    static double p1[8] = { 13524.9996677263464,45285.6016995472897,
	    45135.1684697366626,18529.0118185826102,3329.15251494069355,
	    240.680324743572018,5.15778920001390847,.00622835069189847458 };
    static double q1[8] = { 6.93891117537634444e-7,19768.5742630467364,
	    41255.1608353538323,29390.2871199326819,9081.96660748551703,
	    1244.7477785670856,67.4291295163785938,1. };
    static double p2[5] = { -2.72817575131529678e-15,-.64815712376619651,
	    -4.48616543918019358,-7.01677227766758664,-2.12940445131010517 };
    static double q2[5] = { 7.77788548522961604,54.611773810321507,
	    89.292070048186137,32.2703493791143361,1. };

    /* Format strings */
    static char fmt_101[] = "(\002K = \002,i5,\002  (< 0  OR  > 6)\002)";
    static char fmt_102[] = "(\002ARGUMENT EQUALS NON-POSITIVE int =\002"
	    ",1p,e15.6)";

    /* System generated locals */
    int i__1;
    double ret_val, d__1, d__2, d__3, d__4;

    /* Builtin functions */
	    ;

    /* Local variables */
    static double a, h__;
    static int i__, j;
    static double p, r__, s, v, b0, b1, b2;
    static int k1;
    static double ap, aq;
    static int ix;
    static double alfa;
    static char errtxt[80];

    /* Fortran I/O blocks */




/* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $ */

/* $Log: imp64.inc,v $ */
/* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni */
/* Mathlib gen */


/* imp64.inc */



    a = abs(*x);
    v = a;
    ix = (int) (*x - 1e-13);
    if (*k < 0 || *k > 6) {
	h__ = 0.;
	//s_wsfi(&io___17);
	//do_fio(&c__1, (char *)&(*k), (ftnlen)sizeof(int));
	//e_wsfi();
	//mtlprt_("RPSIPG/DPSIPG", "C316.1", errtxt, (ftnlen)13, (ftnlen)6, (
	//	ftnlen)80);
        throw std::logic_error("npsi large K, not implemented");
    } else if ((d__1 = ix - *x, abs(d__1)) <= 1e-13) {
	h__ = 0.;
	//s_wsfi(&io___18);
	//do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(double));
	//e_wsfi();
        throw std::logic_error("npsi arg near numerical treshold");
	//mtlprt_("RPSIPG/DPSIPG", "C316.2", errtxt, (ftnlen)13, (ftnlen)6, (
	//	ftnlen)80);
    } else if (*k == 0) {
	if (a <= 3.) {
	    s = 0.;
	    if (a < .5) {
		s = 1 / v;
		v += 1;
	    }
	    ap = p1[7];
	    aq = q1[7];
	    for (i__ = 6; i__ >= 0; --i__) {
		ap = p1[i__] + v * ap;
/* L11: */
		aq = q1[i__] + v * aq;
	    }
	    h__ = (v - x0) * ap / aq - s;
	} else {
/* Computing 2nd power */
	    d__1 = v;
	    r__ = 1 / (d__1 * d__1);
	    ap = p2[4];
	    aq = q2[4];
	    for (i__ = 3; i__ >= 0; --i__) {
		ap = p2[i__] + r__ * ap;
/* L12: */
		aq = q2[i__] + r__ * aq;
	    }
	    h__ = log(v) - .5 / v + ap / aq;
	}
	if (*x < 0.) {
	    h__ = h__ + 1 / a + 3.14159265358979324 / tan(a * 
		    3.14159265358979324);
	}
    } else {
	k1 = *k + 1;
	if (a <= 10.) {
	    if (a < 3.) {
		s = -1 / pow(v, k1);
		i__1 = 2 - (int) a;
		for (j = 1; j <= i__1; ++j) {
		    v += 1;
/* L1: */
		    s -= 1 / pow(v, k1);
		}
		v += 1;
	    } else if (a <= 4.) {
		s = 0.;
	    } else {
		v += -1;
		s = 1 / pow(v, k1);
		i__1 = (int) a - 4;
		for (j = 1; j <= i__1; ++j) {
		    v += -1;
/* L5: */
		    s += 1 / pow(v, k1);
		}
	    }
	    h__ = v * 2 - 7;
	    alfa = h__ + h__;
	    b1 = 0.;
	    b2 = 0.;
	    for (j = nb[*k - 1]; j >= 0; --j) {
		b0 = b[j + *k * 21 - 21] + alfa * b1 - b2;
		b2 = b1;
/* L2: */
		b1 = b0;
	    }
	    h__ = b0 - h__ * b2 + sgf[*k] * s;
	} else {
	    s = 0.;
	    if (a < 15.) {
		s = 1 / pow(v, k1);
		i__1 = 14 - (int) a;
		for (j = 1; j <= i__1; ++j) {
		    v += 1;
/* L3: */
		    s += 1 / pow(v, k1);
		}
		v += 1;
	    }
/* Computing 2nd power */
	    d__1 = v;
	    r__ = 1 / (d__1 * d__1);
	    p = r__ * c__[*k * 7 - 1];
	    for (j = 6; j >= 1; --j) {
/* L4: */
		p = r__ * (c__[j + *k * 7 - 8] + p);
	    }
	    h__ = ((sgf[*k - 1] - sgn[*k - 1] * p) * v - sgh[*k - 1]) / 
		    pow(v, k1) - sgf[*k] * s;
	}
	if (*x < 0.) {
	    p = a * 3.14159265358979324;
	    if (*k == 1) {
/* Computing 2nd power */
		d__1 = sin(p);
		v = -9.869604401089358 / (d__1 * d__1);
	    } else if (*k == 2) {
/* Computing 3rd power */
		d__1 = sin(p);
		v = cos(p) * 62.012553360599632 / (d__1 * (d__1 * d__1));
	    } else if (*k == 3) {
/* Computing 2nd power */
		d__1 = sin(p);
		s = d__1 * d__1;
/* Computing 2nd power */
		d__1 = s;
		v = (s * 2 - 3) * 194.81818206800486 / (d__1 * d__1);
	    } else if (*k == 4) {
		s = sin(p);
/* Computing 2nd power */
		d__1 = s;
/* Computing 5th power */
		d__2 = s, d__3 = d__2, d__2 *= d__2;
		v = cos(p) * -2448.1574782822513 * (d__1 * d__1 - 3) / (d__3 *
			 (d__2 * d__2));
	    } else if (*k == 5) {
/* Computing 2nd power */
		d__1 = sin(p);
		s = d__1 * d__1;
/* Computing 2nd power */
		d__1 = s;
/* Computing 3rd power */
		d__2 = s;
		v = (15 - s * 15 + d__1 * d__1 * 2) * -7691.1135486024341 / (
			d__2 * (d__2 * d__2));
	    } else if (*k == 6) {
		s = sin(p);
/* Computing 2nd power */
		d__1 = s;
/* Computing 4th power */
		d__2 = s, d__2 *= d__2;
/* Computing 7th power */
		d__3 = s, d__4 = d__3, d__3 *= d__3, d__4 *= d__3;
		v = cos(p) * 48324.691644428662 * (45 - d__1 * d__1 * 30 + 
			d__2 * d__2 * 2) / (d__4 * (d__3 * d__3));
	    }
	    h__ = sgn[*k - 1] * (h__ + v + sgf[*k] / pow(a, k1));
	}
    }
    ret_val = h__;
#undef b
    return ret_val;
} /* dpsipg_ */



    

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
lst GLOBAL_PL;

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
        //            return gsl_sf_psi_n(int(ex_to<numeric>(n).to_double()),ex_to<numeric>(x).to_double());

        return npsi(ex_to<numeric>(x).to_double(),int(ex_to<numeric>(n).to_double()));
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

//          minimize function approximated with hyperbolas

struct minimize_hyper
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
  GLOBAL_PL = pl;
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
  gsl_multimin_fdfminimizer_set (s, &my_func , x, 0.01, 0.1);
  const size_t ITER_MAX = 1000;
  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
     
      if (status)
        break;
     
      status = gsl_multimin_test_gradient (s->gradient, 1e-3);
     
      if (status == GSL_SUCCESS || iter == ITER_MAX)
        {
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
    }
  while (status == GSL_CONTINUE && iter < ITER_MAX);
     
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
