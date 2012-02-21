#ifndef __UTILS_H__
#define __UTILS_H__

#include <boost/tuple/tuple.hpp>
#include <boost/fusion/container/map.hpp>
#include <ginac/ginac.h>


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


/*
#include <boost/icl/interval.hpp>
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/continuous_interval.hpp>
#include <boost/icl/right_open_interval.hpp>
#include <boost/icl/left_open_interval.hpp>
#include <boost/icl/closed_interval.hpp>
#include <boost/icl/open_interval.hpp>
#include <boost/icl/functors.hpp>
  */
#include <iterator>
#include <functional>

#include <gmpxx.h>


#include "romb_excompiler.h"
#include "collect_square.h"

//#include "sim.h"
#include <cuba.h>


using namespace GiNaC;
using namespace boost;

  typedef std::list<ex> exlist;



lst set2lst(const exset & );

exset lst2set(const lst & );



/***********************************************************************
 *  NAME
 *
 *  fp2rat - convert floating-point number to rational number
 *
 *  SYNOPSIS
 *
 *  #include "glplib.h"
 *  int fp2rat(double x, double eps, double *p, double *q);
 *
 *  DESCRIPTION
 *
 *  Given a floating-point number 0 <= x < 1 the routine fp2rat finds
 *  its "best" rational approximation p / q, where p >= 0 and q > 0 are
 *  integer numbers, such that |x - p / q| <= eps.
 *
 *  RETURNS
 *
 *  The routine fp2rat returns the number of iterations used to achieve
 *  the specified precision eps.
 *
 *  EXAMPLES
 *
 *  For x = sqrt(2) - 1 = 0.414213562373095 and eps = 1e-6 the routine
 *  gives p = 408 and q = 985, where 408 / 985 = 0.414213197969543.
 *
 *  BACKGROUND
 *
 *  It is well known that every positive real number x can be expressed
 *  as the following continued fraction:
 *
 *     x = b[0] + a[1]
 *                ------------------------
 *                b[1] + a[2]
 *                       -----------------
 *                       b[2] + a[3]
 *                              ----------
 *                              b[3] + ...
 *
 *  where:
 *
 *     a[k] = 1,                  k = 0, 1, 2, ...
 *
 *     b[k] = floor(x[k]),        k = 0, 1, 2, ...
 *
 *     x[0] = x,
 *
 *     x[k] = 1 / frac(x[k-1]),   k = 1, 2, 3, ...
 *
 *  To find the "best" rational approximation of x the routine computes
 *  partial fractions f[k] by dropping after k terms as follows:
 *
 *     f[k] = A[k] / B[k],
 *
 *  where:
 *
 *     A[-1] = 1,   A[0] = b[0],   B[-1] = 0,   B[0] = 1,
 *
 *     A[k] = b[k] * A[k-1] + a[k] * A[k-2],
 *
 *     B[k] = b[k] * B[k-1] + a[k] * B[k-2].
 *
 *  Once the condition
 *
 *     |x - f[k]| <= eps
 *
 *  has been satisfied, the routine reports p = A[k] and q = B[k] as the
 *  final answer.
 *
 *  In the table below here is some statistics obtained for one million
 *  random numbers uniformly distributed in the range [0, 1).
 *
 *      eps      max p   mean p      max q    mean q  max k   mean k
 *     -------------------------------------------------------------
 *     1e-1          8      1.6          9       3.2    3      1.4
 *     1e-2         98      6.2         99      12.4    5      2.4
 *     1e-3        997     20.7        998      41.5    8      3.4
 *     1e-4       9959     66.6       9960     133.5   10      4.4
 *     1e-5      97403    211.7      97404     424.2   13      5.3
 *     1e-6     479669    669.9     479670    1342.9   15      6.3
 *     1e-7    1579030   2127.3    3962146    4257.8   16      7.3
 *     1e-8   26188823   6749.4   26188824   13503.4   19      8.2
 *
 *  REFERENCES
 *
 *  W. B. Jones and W. J. Thron, "Continued Fractions: Analytic Theory
 *  and Applications," Encyclopedia on Mathematics and Its Applications,
 *  Addison-Wesley, 1980. */

int fp2rat(double , double, double *, double *);

                                   

mpq_class set_d_eps(double);



ex ginac_set_d(double);


// return value for 1st variable
std::pair<ex,ex>  hyper_cube_den(lst,lst, ex);



// Start point with different countours
exmap start_point_diff_w(lst,lst);
// test if system solvable
bool zero_volume(const lst&,const lst&);

/*
  OPERATORS
 */
template <typename T>
std::ostream & operator<<(std::ostream & os, const std::list<T>& e)
 {
     typename std::list<T>::const_iterator i = e.begin();
     typename std::list<T>::const_iterator vend = e.end();
 
     if (i==vend) {
         os << "()";
         return os;
     }
 
     os << "(";
     while (true) {
       os<<*(i);
         ++i;
         if (i==vend)
             break;
         os << ",";
     }
     os << ")";
 
     return os;
 }

template <typename T1, typename T2, typename C>
    std::ostream & operator<<(std::ostream & os, const std::multimap<T1,T2,C>& e)
 {
     typename std::multimap<T1,T2,C>::const_iterator i = e.begin();
     typename std::multimap<T1,T2,C>::const_iterator vend = e.end();
 
     if (i==vend) {
         os << "()";
         return os;
     }
 
     os << "{";
     while (true) {
       os << i->first << "->" << i->second;
         ++i;
         if (i==vend)
             break;
         os << ",";
     }
     os << "}";
 
     return os;
 }



#endif
