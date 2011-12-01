

/** @file code_generation.cc
 *
 * Implementation of code_generation 
 *
 */

#include <cstddef>
#include <vector>
#include <sstream>

#include "ginac/ginac.h"
#include "code_generation.h"
#include "complex_c_print.h"
namespace GiNaC {

  /**
   *
   * Constructor.
   *
   * T is initialised with arg_T, temp_var_base_name is initialised with arg_temp_var_base_name.
   * lhs and rhs are empty.
   *
   */
  assign_lst::assign_lst(const std::string & arg_T, const std::string & arg_temp_var_base_name) :
    lhs(lst()), rhs(lst()), T(arg_T), temp_var_base_name(arg_temp_var_base_name), temp_var()
  {}

  /**
   *
   * Constructor.
   *
   * T is initialised with arg_T, temp_var_base_name is initialised with arg_temp_var_base_name.
   * lhs and rhs are initialised with arg_lhs and arg_rhs, respectively.
   *
   */
  assign_lst::assign_lst(const std::string & arg_T, const std::string & arg_temp_var_base_name, const lst & arg_lhs, const lst & arg_rhs) :
    lhs(arg_lhs), rhs(arg_rhs), T(arg_T), temp_var_base_name(arg_temp_var_base_name), temp_var()
  {}

  /**
   *
   * Appends the assignment
   * a = b
   * to the assignment list.
   *
   */
  void assign_lst::append(const ex & a, const ex & b)
  {
    lst lhs_lst = ex_to<lst>(lhs);
    lst rhs_lst = ex_to<lst>(rhs);

    lhs = lhs_lst.append(a);
    rhs = rhs_lst.append(b);
  }

  /**
   *
   * Prepends the assignment
   * a = b
   * to the assignment list.
   *
   */
  void assign_lst::prepend(const ex & a, const ex & b)
  {
    lst lhs_lst = ex_to<lst>(lhs);
    lst rhs_lst = ex_to<lst>(rhs);

    lhs = lhs_lst.prepend(a);
    rhs = rhs_lst.prepend(b);
  }

  /**
   *
   * This function substitutes all subexpressions which occur syntactially more than once 
   * (in the sense of "subs") in the assignment list.
   *
   * The function proceeds iteratively and substitutes first subexpressions, which are
   * "closer" to the root of the expression tree.
   * If these subexpressions in turn contain subsubexpressions, which occur multiple times,
   * the latter are substituted in subsequent iterations.
   *
   */
  void assign_lst::replace_common_subex(unsigned flag_type)
  {
    bool flag_try_substitution = true;

    lst replies;

    // initialize
    lst lhs_lst = ex_to<lst>(lhs);
    lst rhs_lst = ex_to<lst>(rhs);

    while ( flag_try_substitution )
      {
	// expand tree into a long list
	lst subexpr_lst;
	for (const_preorder_iterator it = rhs.preorder_begin(); it != rhs.preorder_end(); it++)
	  {
	    if ( check_substitution_type(*it, flag_type) ) subexpr_lst.append(*it);
	  }
	subexpr_lst.sort();

	// find dublicate entries
	replies.remove_all();
	// 08.05.2007: subexpr_lst might be empty
	if ( subexpr_lst.nops() > 0 )
	  {
	    ex subexpr = subexpr_lst.op(0);
	    int counter = 1;

	    lst::const_iterator it_start = subexpr_lst.begin();
	    it_start++;
	    for ( lst::const_iterator it = it_start; it != subexpr_lst.end(); it++ )
	      {
		if ( subexpr.is_equal(*it) ) 
		  {
		    counter++;

		    if ( counter == 2 )
		      {
			bool flag_subs = true;

			size_t k = 0;
			for ( lst::const_iterator k_it = replies.begin(); k_it != replies.end(); k_it++ )
			  {
			    // subexpression of something already in replies
			    if ( (*k_it).has(subexpr) ) 
			      {
				flag_subs = false;
			      }
			    // replies contains a subsubexpression of subexpr
			    // note that subexpr can be the mother of more than one child
			    else if ( subexpr.has(*k_it) )
			      {
				replies.let_op(k) = subexpr;
				flag_subs = false;
			      }

			    // update k
			    k++;
			  } // const_iterator k_it

			// new substitution
			if ( flag_subs )
			  {
			    replies.append(subexpr);
			  }
		      } // counter == 2

		  } // equal
		else // not equal
		  {
		    // new subexpression
		    subexpr = *it;
		    counter = 1;
		  } // subexpr not equal
	      } // lst::const_iterator it
	  } // subexpr_lst.nops() > 0

	// sort and remove double entries
	replies.sort();
	replies.unique();

	// substitute
	if ( replies.nops() > 0 )
	  {
	    exmap subs_map;

	    for ( lst::const_iterator k = replies.begin(); k != replies.end(); k++ )
	      {
		int i = temp_var.size();
		temp_var.push_back( symbol( temp_var_base_name+itos(i) ) );
		subs_map[*k] = temp_var[i];
	      }
	    rhs = rhs.subs( subs_map, subs_options::no_pattern );

	    rhs_lst = ex_to<lst>(rhs);
	    for ( lst::const_iterator k = replies.begin(); k != replies.end(); k++ )
	      {
		lhs_lst.prepend(subs_map[*k]);
		rhs_lst.prepend(*k);
	      }
	    lhs = lhs_lst;
	    rhs = rhs_lst;
	  }
	else
	  {
	    flag_try_substitution = false;
	  }

      } // flag_try_substitution
  }
   
  /**
   *
   *  This method splits subexpressions with more than max_nops operands into smaller
   *  pieces.
   *
   */
  void assign_lst::split_large_subex(size_t max_nops)
  {
    lst new_lhs_lst;
    lst new_rhs_lst;

    lst lhs_lst = ex_to<lst>(lhs);
    lst rhs_lst = ex_to<lst>(rhs);

    lst::const_iterator it_lhs, it_rhs;
    for ( it_lhs = lhs_lst.begin(), it_rhs = rhs_lst.begin(); it_lhs != lhs_lst.end(); it_lhs++, it_rhs++)
      {
	split_large_line( *it_lhs, *it_rhs, new_lhs_lst, new_rhs_lst, max_nops);
      } 

    lhs = new_lhs_lst;
    rhs = new_rhs_lst;
  }

  /**
   *
   * Helper function, which converts an integer into a string
   *
   */
  std::string assign_lst::itos(int arg) const
  {
    std::ostringstream buffer;
    buffer << arg; 
    return buffer.str(); 
  }

  /**
   *
   * Checks if an expression can potentially be substituted.
   *
   * The flag "flag_type" determines which classes from the
   * set (add, mul, power, function) are considered for this substitution.
   *
   */
  bool assign_lst::check_substitution_type(const ex & expr, unsigned flag_type) const
  {
    if ( ( is_a<add>(expr) && (flag_type & subex_type::add) ) 
	 || ( is_a<mul>(expr) && (flag_type & subex_type::mul) ) 
	 || ( is_a<power>(expr) && (flag_type & subex_type::power) ) 
	 || ( is_a<function>(expr) && (flag_type & subex_type::function) ) 
	 ) return true;

    return false;
  }

  /**
   *
   * Checks if expr contains large subexpressions of type add or mul with more
   * than max_nops terms.
   *
   * These large subexpressions are appended to subexpr_lst.
   *
   */
  bool assign_lst::check_large_subex(const ex & expr, size_t max_nops, lst & subexpr_lst) const
  {
    bool flag_long_expression = false;

    for (const_preorder_iterator it = expr.preorder_begin(); it != expr.preorder_end(); it++)
      {
	if ( is_a<add>(*it) || is_a<mul>(*it) ) 
	  {
	    if ( (*it).nops() > max_nops )
	      {
		bool flag_subs = true;

		size_t k = 0;
		for ( lst::const_iterator k_it = subexpr_lst.begin(); k_it != subexpr_lst.end(); k_it++ )
		  {
		    // subexpression of something already in replies
		    if ( k_it->has( *it ) ) 
		      {
			flag_subs = false;
		      }
		    // replies contains a subsubexpression of subexpr
		    // note that subexpr can be the mother of more than one child
		    else if ( it->has(*k_it) )
		      {
			subexpr_lst.let_op(k) = *it;
			flag_subs = false;
		      }

		    // update k
		    k++;
		  } // const_iterator k_it

		    // new substitution
		if ( flag_subs )
		  {
		    subexpr_lst.append(*it);
		    flag_long_expression = true;
		  }

	      } // (*it).nops() > max_nops
	  } // is_a<add>(*it) || is_a<mul>(*it)
      } // const_preorder_iterator it

    // sort and remove double entries
    subexpr_lst.sort();
    subexpr_lst.unique();

    return flag_long_expression;
  }

  /**
   *
   * Splits a single assignment into smaller pieces, if the right-hand-side contains a subexpression
   * with more than max_nops operators.
   *
   * The algorithm proceeds recursively, by introducing temporary variables for the first max_nops operands
   * of the large subexpression and the remaining operands of the subexpression.
   *
   * The assignments are appended to new_lhs_lst and new_rhs_lst.
   * In the trivial case, where b does not contain any large subexpression, the method only appends
   * a to new_lhs_lst and b to new_rhs_lst.
   *
   */
  void assign_lst::split_large_line(const ex & a, const ex & b, lst & new_lhs_lst, lst & new_rhs_lst, size_t max_nops)
  {
    lst subexpr_lst;
    if ( check_large_subex( b, max_nops, subexpr_lst) )
      {
	exmap subs_map;

	for ( lst::const_iterator k = subexpr_lst.begin(); k != subexpr_lst.end(); k++ )
	  {
	    // two new symbols
	    size_t i1 = temp_var.size();
	    temp_var.push_back( symbol( temp_var_base_name+itos(i1) ) );
	    size_t i2 = temp_var.size();
	    temp_var.push_back( symbol( temp_var_base_name+itos(i2) ) );

	    ex expr_1, expr_2;

	    if ( is_a<add>(*k) )
	      {
		expr_1 = 0;
		const_iterator j_it = k->begin();
		for (size_t j=0; j<max_nops; j++)
		  {
		    expr_1 += *j_it;
		    j_it++;
		  }
		expr_2 = (*k) - expr_1;
 
		subs_map[*k] = temp_var[i1]+temp_var[i2];
	      }
	    else if ( is_a<mul>(*k) )
	      {
		expr_1 = 1;
		const_iterator j_it = k->begin();
		for (size_t j=0; j<max_nops; j++)
		  {
		    expr_1 *= *j_it;
		    j_it++;
		  }
		expr_2 = (*k) / expr_1;
 
		subs_map[*k] = temp_var[i1]*temp_var[i2];
	      }

	    // print expr_2 and expr_1
	    split_large_line( temp_var[i2], expr_2, new_lhs_lst, new_rhs_lst, max_nops);
	    split_large_line( temp_var[i1], expr_1, new_lhs_lst, new_rhs_lst, max_nops);

	  } // lst::const_iterator k

	new_lhs_lst.append(a);
	new_rhs_lst.append(b.subs( subs_map, subs_options::no_pattern ));
      }
    else // no large subex
      {
	new_lhs_lst.append(a);
	new_rhs_lst.append(b);
      }
  }

  /**
   *
   * Returns a list of the temporary variables.
   *
   */
  ex assign_lst::get_temporary_variables() const
  {
    lst symbol_lst;

    for (size_t j=0; j<temp_var.size(); j++) symbol_lst.append(temp_var[j]);

    return symbol_lst;
  }


  /**
   *    
   * Print method: Prints the assignment list as C code.
   *
   */
  void assign_lst::print(std::ostream & os) const
  {
    lst lhs_lst = ex_to<lst>(lhs);
    lst rhs_lst = ex_to<lst>(rhs.subs(psi(wild(1),wild(2)) == wpsipg(wild(2),wild(1))).subs(psi(wild(1)) 
                                                                                            == wpsipg(wild(1),0)).subs(tgamma(wild(3))==wgamma(wild(3))));
set_print_func<numeric, print_csrc_double>(complex_fort_print);
set_print_func<power, print_csrc_double>(power_fort_print);

    lst::const_iterator it_lhs, it_rhs;
    for ( it_lhs = lhs_lst.begin(), it_rhs = rhs_lst.begin(); it_lhs != lhs_lst.end(); it_lhs++, it_rhs++)
      os<<"     double complex "<<*it_lhs<<std::endl;
    for ( it_lhs = lhs_lst.begin(), it_rhs = rhs_lst.begin(); it_lhs != lhs_lst.end(); it_lhs++, it_rhs++)
      {
	os << " " 
          //<< T << " " 
           << *it_lhs << " = ";

	if ( T == std::string("double") )
	  {
	    it_rhs->print(print_csrc_double(os));
	  }
	else if ( T == std::string("float") )
	  {
	    it_rhs->print(print_csrc_double(os));
	  }
	else if ( T == std::string("cl_N") )
	  {
	    it_rhs->print(print_csrc_cl_N(os));
	  }
	else // default
	  {
	    it_rhs->print(print_csrc_double(os));
	  }

	os << ";" << std::endl;
      }
  }

  /**
   *
   * Output operator. Prints the assignment list as C code.
   *
   */
  std::ostream & operator<< (std::ostream & os, const assign_lst & arg)
  {
    arg.print(os);

    return os;
  }

} // namespace GiNaC

  
