#include <ginac/ginac.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <functional>
#include <exception>

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"

#include <boost/numeric/ublas/matrix_proxy.hpp> /* matrix_row */
#include <boost/numeric/ublas/vector.hpp> /*vector*/
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/icl/continuous_interval.hpp>
#include  <boost/icl/open_interval.hpp>



using namespace GiNaC;
using namespace std;
using GiNaC::symbol;
using GiNaC::lst;
using GiNaC::exvector;
using std::cout;
using std::endl;
using namespace boost::numeric::ublas;

//typedef boost::numeric::ublas::vector_of_vector<GiNaC::ex> rmatrix;
typedef boost::numeric::ublas::matrix<GiNaC::ex> rmatrix;



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


// all inequalities in for A*x>0, fcn c*x>0

 rmatrix standart_form(lst ineq_lst, lst x_lst, ex fcn)
{
  lst new_x_lst;
  exmap wpm_map;
  for(lst::const_iterator it = x_lst.begin(); it != x_lst.end(); ++it)
    {
      if(is_a<symbol>(*it))
	{
	  cout<<"dfgdfg"<<endl;
	  string str(ex_to<symbol>(*it).get_name());
          new_x_lst.append(get_symbol(str + "p"));
          new_x_lst.append(get_symbol(str + "m"));
	  wpm_map[*it] = get_symbol(str + "p") - get_symbol(str + "m");
	  cout<<str<<endl;
	}
    }
  cout<<wpm_map<<endl;
  lst new_ineq_lst =  ex_to<lst>(ineq_lst.subs(wpm_map));
                   //ineq_lst = ineq_lst.subs(wpm_map);
  cout<<new_ineq_lst<<endl;
  lst eq_lst;
  size_t sj = 0;
  for(lst::const_iterator it = new_ineq_lst.begin(); it != new_ineq_lst.end(); ++it)
    {
      ex b = *it;
     
      for(lst::const_iterator xit = new_x_lst.begin(); xit != new_x_lst.end(); ++xit)
	{
	  //T(j,i) = -1*(it->coeff(*xit,1));  /* (row,col) */
          b = b.coeff(*xit,0);
	}
      //      ex to_eq_lst = 
      
      if(b > 0)
        {
          ex s_i = get_symbol("s" + boost::lexical_cast<string>(sj));
          new_x_lst.append(s_i);
          eq_lst.append( *it-s_i);
          sj++;
        }
      else if(b < 0)
        {
          ex s_i = get_symbol("s" + boost::lexical_cast<string>(sj));
          new_x_lst.append(s_i);
          eq_lst.append((-1)*(*it-s_i)); 
          sj++;
        }
      else
        eq_lst.append(*it); 
      
      //Q(j) = b;
      //j++;
    }
  cout<<"EQ_LST: "<<endl<<eq_lst<<endl;  

  rmatrix T (1 + eq_lst.nops(), 1 + new_x_lst.nops() + 1);
  size_t r = T.size1();
  size_t s = T.size2();
  cout<< r<< "   "<<s<<endl;
  int i = 1,j = 1;
  T(0,0) = 1;
  for(lst::const_iterator it = eq_lst.begin(); it != eq_lst.end(); ++it)
    {
      ex b = *it;
      i = 1;
      for(lst::const_iterator xit = new_x_lst.begin(); xit != new_x_lst.end(); ++xit)
	{
	  T(j,i) = it->coeff(*xit,1);  /* (row,col) */
	  i++;
	  b = b.coeff(*xit,0);
	}
       T(j,r) = b;
      j++;
    }
  i = 1;
  // add function
  ex fcn_subs = fcn.subs(wpm_map);
  for(lst::const_iterator xit = new_x_lst.begin(); xit != new_x_lst.end(); ++xit)
    {
      T(0,i) = fcn_subs.coeff(*xit,1);  /* (row,col) */
      i++;
    }
  return T;
}



int main () 
{
  symbol w1("w1"),w2("w2"),w3("w3"),w4("w4"),w5("w5"),ep("ep");
  //bound_vec bv;
  //  fm(lst(1-w1+2*w2,2-7*w3-w2,3-w1+w2-3*w3,4-w1+5*w3,5-w2,6-2*w3,w1+30), lst(w1,w2,w3), bv);
  //  fm(lst(w1,3-w1), lst(w1), bv);
  //  fm(lst(-w1,-w2,-w3,2+ep+w1+w2+w3,1+w3,1+w1+w2+w3,-1-ep-w1-w3,-1-ep-w2-w3), lst(w1,w2,w3,ep), bv);
  //  standart_form(lst(-w1,-w2,-w3,2+ep+w1+w2+w3,1+w3,1+w1+w2+w3,-1-ep-w1-w3,-1-ep-w2-w3), lst(w1,w2,w3,ep),ep);
  cout<< standart_form(lst(w1-5,-w2-2*w3+3,-w4+3*w5-2), lst(w1,w2,w3,w4,w5),w5)<<endl;

  
}






