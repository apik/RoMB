#ifndef __COMPLEX_C_PRINT_H__
#define __COMPLEX_C_PRINT_H__
#include <limits>
#include "ginac/ginac.h"
using namespace GiNaC;
#include <cln/output.h>
#include <cln/integer_io.h>
#include <cln/integer_ring.h>
#include <cln/rational_io.h>
#include <cln/rational_ring.h>
#include <cln/lfloat_class.h>
#include <cln/lfloat_io.h>
#include <cln/real_io.h>
#include <cln/real_ring.h>
#include <cln/complex_io.h>
#include <cln/complex_ring.h>
#include <cln/numtheory.h>

static void print_integer_csrc(const print_context & c, const cln::cl_I & x)
{
  // Print small numbers in compact float format, but larger numbers in
  // scientific format
  const int max_cln_int = 536870911; // 2^29-1
  if (x >= cln::cl_I(-max_cln_int) && x <= cln::cl_I(max_cln_int))
    c.s << cln::cl_I_to_int(x) << ".0";
  else
    c.s << cln::double_approx(x);
}


static void print_real_csrc(const print_context & c, const cln::cl_R & x)
{
  if (cln::instanceof(x, cln::cl_I_ring)) {

    // Integer number
    print_integer_csrc(c, cln::the<cln::cl_I>(x));

  } else if (cln::instanceof(x, cln::cl_RA_ring)) {

    // Rational number
    const cln::cl_I numer = cln::numerator(cln::the<cln::cl_RA>(x));
    const cln::cl_I denom = cln::denominator(cln::the<cln::cl_RA>(x));
    if (cln::plusp(x) > 0) {
      c.s << "(";
      print_integer_csrc(c, numer);
    } else {
      c.s << "-(";
      print_integer_csrc(c, -numer);
    }
    c.s << "/";
    print_integer_csrc(c, denom);
    c.s << ")";

  } else {

    // Anything else
    c.s << cln::double_approx(x);
  }
}


void complex_fort_print(const numeric & value,const print_csrc & c, unsigned level) 
{
  std::ios::fmtflags oldflags = c.s.flags();
  c.s.setf(std::ios::scientific);
  int oldprec = c.s.precision();
  
  // Set precision
  if (is_a<print_csrc_double>(c))
    c.s.precision(std::numeric_limits<double>::digits10 + 1);
  else
    c.s.precision(std::numeric_limits<float>::digits10 + 1);
  
  if (value.is_real()) {
    
    // Real number
    print_real_csrc(c, cln::the<cln::cl_R>(value.to_cl_N()));
    
  } else {
    // Complex number
    /*    c.s << "std::cooooooooooomplex<";
          if (is_a<print_csrc_double>(c))
          c.s << "double>(";
          else
          c.s << "float>(";*/
    /*  
  if( zerop(cln::realpart(value.to_cl_N())))
      {
        print_real_csrc(c, cln::imagpart(value.to_cl_N()));
        c.s << "*I";
      }
    else if(zerop(cln::imagpart(value.to_cl_N())))
      print_real_csrc(c, cln::realpart(value.to_cl_N()));
    else
      {
        print_real_csrc(c, cln::realpart(value.to_cl_N()));
        c.s << "+";
        print_real_csrc(c, cln::imagpart(value.to_cl_N()));
        c.s << "*I";
      }
*/
    c.s<<"(";
    print_real_csrc(c, cln::realpart(value.to_cl_N()));
    c.s << ",";
    print_real_csrc(c, cln::imagpart(value.to_cl_N()));
    c.s << ")";
  }
  
  c.s.flags(oldflags);
  c.s.precision(oldprec);
}


static void print_sym_pow(const print_context & c, const symbol &x, int exp)
{
        // Optimal output of integer powers of symbols to aid compiler CSE.
        // C.f. ISO/IEC 14882:1998, section 1.9 [intro execution], paragraph 15
        // to learn why such a parenthesation is really necessary.
        if (exp == 1) {
                x.print(c);
        } else if (exp == 2) {
                x.print(c);
                c.s << "*";
                x.print(c);
        } else if (exp & 1) {
                x.print(c);
                c.s << "*";
                print_sym_pow(c, x, exp-1);
        } else {
                c.s << "(";
                print_sym_pow(c, x, exp >> 1);
                c.s << ")*(";
                print_sym_pow(c, x, exp >> 1);
                c.s << ")";
        }
}


// POWERS PRINTIG AS **
void power_fort_print(const power & value,const print_csrc & c, unsigned level) 
{
  ex basis = value.op(0);
  ex exponent = value.op(1);
  // Integer powers of symbols are printed in a special, optimized way
  if (exponent.info(info_flags::integer)
      && (is_a<symbol>(basis) || is_a<constant>(basis))) {
    int exp = ex_to<numeric>(exponent).to_int();
    if (exp > 0)
      c.s << '(';
    else {
      exp = -exp;
      c.s << "1.0/(";
    }
    print_sym_pow(c, ex_to<symbol>(basis), exp);
    c.s << ')';
                                                                                                                                                                                                 
    // <expr>^-1 is printed as "1.0/<expr>" or with the recip() function of CLN
  } else if (exponent.is_equal(-1)) {
    c.s << "1.0/(";
    basis.print(c);
    c.s << ')';
                                                                                                                                                                                                                                                                 
    // Otherwise, use the pow() function
  } else {
    c.s << "(";
    basis.print(c);
    c.s << ")**(";
    exponent.print(c);
    c.s << ')';
  }
}
                                                                                                                                                                                                                                                                                                                                                                         
DECLARE_FUNCTION_1P(wgamma)
REGISTER_FUNCTION(wgamma, dummy())

DECLARE_FUNCTION_2P(wpsipg)
REGISTER_FUNCTION(wpsipg, dummy())
#endif
