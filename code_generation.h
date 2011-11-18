
/** @file code_generation.h
 *
 * Interface to code_generation
 *
 */

#ifndef __GINAC_CODE_GENERATION_H__
#define __GINAC_CODE_GENERATION_H__

#include <iostream>
#include <string>

#include "ginac/ex.h"
#include "ginac/lst.h"

namespace GiNaC {

  /**
   *
   * An enumeration class, which can be used as flags to determine the behaviour of
   * assign_lst::replace_common_subex(unsigned flag).
   *
   */
  class subex_type {
  public:
    enum {
      add        = 0x00000001,
      mul        = 0x00000002,
      power      = 0x00000004,
      function   = 0x00000008,
      all        = 0xFFFFFFFF
    };
  };

  /**
   *
   * A class which transforms a list of assignments into an equivalent list of assignments.
   * The intention is that when these assignments are printed out as C code, the transformed ones
   * will give "better" C code.
   *
   * Methods to manipulate the assignment list are:
   *  - replace_common_subex(unsigned flag), which replaces subexpressions, which occur syntactically 
   *    more than once by temporary variables. The flag determines which classes from the
   *    set (add, mul, power, function) are considered for this substitution.
   *    The default is all of them.
   * 
   *  - split_large_subex(size_t nmax), splits subexpressions with more than nmax operands into smaller
   *    pieces.
   *
   * Data members:
   * A list of assignments 
   * \verbatim
      t0 = expr0;
      t1 = expr1;
       ...
      tn = exprn;
     \endverbatim
   * is stored in two lists
   * \verbatim
      lhs = lst(t0,t1,...,tn);
      rhs = lst(expr0,expr1,...,exprn);
     \endverbatim
   * T is a string giving the data type in C of the expressions, usually something like "double" or "float".
   * When printing out the assignment list, T is printed in front of every assignment:
   * \verbatim
      T t0 = expr0;
      T t1 = expr1;
       ...
     \endverbatim
   * T can also be the empty string.
   *
   * Intermediate variables obtain their names from temp_var_base_name, to which a serial number is attached, e.g.
   * if temp_var_base_name=string("t"), the intermediate variables are labelled t0, t1, t2, etc..
   *
   * temp_var is used only internally and is a vector holding the intermediate variables.
   * A list of the additional intermediate variables can be obtained with the method get_temporary_variables().
   *
   */
  class assign_lst {

  public:
    assign_lst(const std::string & T, const std::string & temp_var_base_name);
    assign_lst(const std::string & T, const std::string & temp_var_base_name, const lst & lhs, const lst & rhs);

  public:
    void append(const ex & a, const ex & b);
    void prepend(const ex & a, const ex & b);
    void replace_common_subex(unsigned flag_type = subex_type::all);
    void split_large_subex(size_t nmax);
    void print(std::ostream & os) const;
    ex get_temporary_variables() const;

  private:
    std::string itos(int arg) const;
    bool check_substitution_type(const ex & expr, unsigned flag_type) const;
    bool check_large_subex(const ex & expr, size_t max_nops, lst & subexpr_lst) const;
    void split_large_line(const ex & a, const ex & b, lst & new_lhs_lst, lst & new_rhs_lst, size_t max_nops);

  public:
    // I/O operators
    friend std::ostream & operator<< (std::ostream & os, const assign_lst & arg);


  private:
    /// a list containing the l.h.s. of the assignment list
    ex lhs; 
    /// a list containing the r.h.s. of the assignment list
    ex rhs;

    /// string giving the data type in C of the expressions, usually something like "double" or "float".
    std::string T; 
    /// string used as a base for labelling temporary variables.
    std::string temp_var_base_name;

    /// vector holding intermediate variables, used internally.
    std::vector<symbol> temp_var;
  };

  std::ostream & operator<< (std::ostream & os, const assign_lst & arg);

} // namespace GiNaC

#endif // ndef __GINAC_CODE_GENERATION_H__

