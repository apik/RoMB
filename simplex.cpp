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

void standart_form(lst ineq_lst, lst x_lst, ex fcn)
{
  exmap wpm_map;
  for(lst::const_iterator it = x_lst.begin(); it != x_lst.end(); ++it)
    {
      if(is_a<symbol>(*it))
	{
	  cout<<"dfgdfg"<<endl;
	  string str(ex_to<symbol>(*it).get_name());
	  wpm_map[*it] = get_symbol(str + "p") - get_symbol(str + "m");
	  cout<<str<<endl;
	}
    }
  cout<<wpm_map<<endl;
  cout<<ineq_lst.subs(wpm_map)<<endl;
}


int main () 
{
  symbol w1("w1"),w2("w2"),w3("w3"),ep("ep");
  //bound_vec bv;
  //  fm(lst(1-w1+2*w2,2-7*w3-w2,3-w1+w2-3*w3,4-w1+5*w3,5-w2,6-2*w3,w1+30), lst(w1,w2,w3), bv);
  //  fm(lst(w1,3-w1), lst(w1), bv);
  //  fm(lst(-w1,-w2,-w3,2+ep+w1+w2+w3,1+w3,1+w1+w2+w3,-1-ep-w1-w3,-1-ep-w2-w3), lst(w1,w2,w3,ep), bv);
  standart_form(lst(-w1,-w2,-w3,2+ep+w1+w2+w3,1+w3,1+w1+w2+w3,-1-ep-w1-w3,-1-ep-w2-w3), lst(w1,w2,w3,ep),ep);

  
}
