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


void complex_c_print(const numeric & value,const print_csrc & c, unsigned level) 
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
  }
  
  c.s.flags(oldflags);
  c.s.precision(oldprec);
}

DECLARE_FUNCTION_1P(dummy_psi1)
REGISTER_FUNCTION(dummy_psi1, dummy())

DECLARE_FUNCTION_2P(dummy_psi2)
REGISTER_FUNCTION(dummy_psi2, dummy())
#endif
