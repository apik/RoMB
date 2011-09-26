/*
 * Fourier-Motzkin elimination algoritmAlgorithm described in 
 * Ke.ler, Christoph W.. "Parallel Fourier.Motzkin Elimination". Universit.t Trier.
 * http://en.scientificcommons.org/43069420
 */
#include "fm.h"
#include <iostream>
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



// We can change the library default for the interval types by defining 
//#define BOOST_ICL_USE_STATIC_BOUNDED_INTERVALS
// prior to other inluces from the icl.
// The interval type that is automatically used with interval
// containers then is the statically bounded right_open_interval.

#include <boost/icl/interval_set.hpp>
//#include <boost/icl/split_interval_set.hpp>
// The statically bounded interval type 'right_open_interval'
// is indirectly included via interval containers.


using GiNaC::ex;
using GiNaC::symbol;
using GiNaC::lst;
using GiNaC::exvector;
using std::cout;
using std::endl;
using namespace boost::numeric::ublas;

//typedef boost::numeric::ublas::vector_of_vector<GiNaC::ex> rmatrix;
typedef boost::numeric::ublas::matrix<GiNaC::ex> rmatrix;

/* true (swap) >0,<0,==0 */
bool slc_comp(GiNaC::ex a,GiNaC::ex b)
{
  return (b <0 && a == 0) || (b > 0 && a <= 0);
}

/* true (swap) >0,<0,==0 */
bool not_slc_comp(GiNaC::ex a,GiNaC::ex b)
{
  return !((b <0 && a == 0) || (b > 0 && a <= 0));
}

struct f_not_slc_comp
{
  bool operator()(GiNaC::ex a,GiNaC::ex b)
  {
    return !((b <0 && a == 0) || (b > 0 && a <= 0));
  }
};

/* true (swap) >0,<0,==0 */
bool map_not_slc_comp(std::pair<GiNaC::ex,size_t> a,std::pair<GiNaC::ex,size_t> b)
{
  return !((b.first <0 && a.first == 0) || (b.first > 0 && a.first <= 0));
}


//std::pair<size_t,size_t> sort_last_column(boost::numeric::ublas::matrix<GiNaC::ex>& A, boost::numeric::ublas::vector<GiNaC::ex>& B)
std::pair<size_t,size_t> sort_last_column(rmatrix& A, boost::numeric::ublas::vector<GiNaC::ex>& B)
{
  using namespace boost::numeric::ublas;
  int r = A.size2() - 1;
  bool changed = true;
  while(changed)
    {
      changed = false;
      for(int j = 0; j < A.size1() - 1; j++)
	{
	  if(slc_comp(A(j,r),A(j + 1,r)))
	    {
	      matrix_row<rmatrix>  rowk(A, j);
	      matrix_row<rmatrix>  rowl(A, j + 1);
	      rowk.swap(rowl);
	      std::swap(B(j),B(j + 1));
	      changed = true;
	      //	      break;
	    }
	  //else continue;
	  // break;
	}
    }
  matrix_column<rmatrix > Acol(A, r);
  cout<<Acol<<endl;
  size_t n1 = std::count_if (Acol.begin(), Acol.end(), std::bind1st(std::less<GiNaC::ex>(),0) );
  size_t n2 = n1 + std::count_if (Acol.begin(), Acol.end(), std::bind1st(std::greater<GiNaC::ex>(),0) );

  cout<<"Zero terms "<<std::count_if (Acol.begin(), Acol.end(), std::bind1st(std::equal_to<GiNaC::ex>(),0) )<<endl;
  return std::make_pair(n1,n2);
}



std::pair<size_t,size_t> fast_sort_last_column(rmatrix& A, boost::numeric::ublas::vector<GiNaC::ex>& B)
{
  using namespace boost::numeric::ublas;
  //cout<<"before sort: "<<A<<endl;
  rmatrix A_out(A.size1(), A.size2());
  boost::numeric::ublas::vector<GiNaC::ex> B_out(B.size());
  int r = A.size2() - 1;
  /*
  bool changed = true;
  while(changed)
    {
      changed = false;
      for(int j = 0; j < A.size1() - 1; j++)
	{
	  if(slc_comp(A(j,r),A(j + 1,r)))
	    {
	      matrix_row<rmatrix > rowk(A, j);
	      matrix_row<rmatrix > rowl(A, j + 1);
	      rowk.swap(rowl);
	      std::swap(B(j),B(j + 1));
	      changed = true;
	      //	      break;
	    }
	    else continue;
	    break;
	}
    }
  */

  matrix_column<rmatrix > Acol(A, r);
  std::multimap<GiNaC::ex,size_t,f_not_slc_comp> last_row_sort_map;
  for(matrix_column<rmatrix >::iterator it = Acol.begin(); it != Acol.end(); ++it)
    {
      last_row_sort_map.insert(std::pair<GiNaC::ex,size_t>(*it,it.index()));
      //      cout<< it.index()<<endl;
    }
  size_t cntr = 0;
  for(std::multimap<GiNaC::ex,size_t>::iterator it = last_row_sort_map.begin();it != last_row_sort_map.end(); ++it)
    {
      size_t idx = size_t(it->second);
      //      matrix_row<rmatrix > A_out_row(A_out, cntr);
      matrix_row<rmatrix > rowk(A, idx);
      matrix_row<rmatrix > rowl(A_out, cntr);
      //   cout<<rowk<<endl;
      // cout<<rowl<<endl;

      rowl = rowk;
      //std::swap(B(cntr), B(idx));
      B_out(cntr) = B(idx);
      cntr++;
      //  cout<<"MAMAMAMPPP "<<it->first<<" sec: "<<it->second<<endl;
    }

//stable_sort(Acol.begin(), Acol.end(),not_slc_comp);
//cout<<Acol<<endl;
  size_t n1 = std::count_if (Acol.begin(), Acol.end(), std::bind1st(std::less<GiNaC::ex>(),0) );
  size_t n2 = n1 + std::count_if (Acol.begin(), Acol.end(), std::bind1st(std::greater<GiNaC::ex>(),0) );

  //  cout<<"Zero terms "<<std::count_if (Acol.begin(), Acol.end(), std::bind1st(std::equal_to<GiNaC::ex>(),0) )<<endl;
  A = A_out;
  B = B_out;


  return std::make_pair(n1,n2);
}

typedef boost::tuples::tuple<GiNaC::ex, GiNaC::exvector, GiNaC::exvector> bound;
typedef std::vector<bound> bound_vec;
/*
  Input is a list ofinequalities a_i_j*x_j+b >= 0
  and list of independent variables x_i
*/
void fm_instance(lst ineq_lst, lst x_lst, bound_vec&  xb)
{
  using namespace boost::numeric::ublas;
  rmatrix T (ineq_lst.nops(),x_lst.nops());  /* (rows,cols) */
  vector<GiNaC::ex> Q (ineq_lst.nops());
  int i = 0,j = 0;
  for(lst::const_iterator it = ineq_lst.begin(); it != ineq_lst.end(); ++it)
    {
      ex b = *it;
      i = 0;
      for(lst::const_iterator xit = x_lst.begin(); xit != x_lst.end(); ++xit)
	{
	  T(j,i) = -1*(it->coeff(*xit,1));  /* (row,col) */
	  i++;
	  b = b.coeff(*xit,0);
	}
      Q(j) = b;
      j++;
    }
  size_t n1, n2;

  int r = T.size2() - 1;
  int s = T.size1() - 1;
  cout<<endl<<" START SORT for r = "<<r<<endl<<endl;
  boost::tie(n1,n2) = fast_sort_last_column(T,Q);
    cout<<endl<<" END SORT LAST"<<endl<<endl;

  // Normalizing inequalities with x_r
  for(int i = 0; i < n2; i++)
    {
      matrix_row<rmatrix > rowi(T, i);
      Q(i) /= T(i,r);
      rowi /= T(i,r);
    }
  vector<GiNaC::ex> xvec(x_lst.nops() - 1);
  for( int i = 0; i != x_lst.nops() - 1;i++)
    {
      xvec(i) = x_lst.op(i);
    }    
  //cout<<"xvec "<<xvec<<endl;
  //std::cout<< T<<std::endl;
  //std::cout<< Q<<endl;

  // U nad L bounds on x_r
  exvector Lbound,Ubound;
  for(int i = 0; i < n1; i++)
    {
      matrix_vector_slice<rmatrix > mvs (T, slice (i, 0, r), slice (0, 1, r));
      Ubound.push_back(Q(i)-inner_prod(mvs,xvec));
      //cout<<inner_prod(mvs,xvec)<<endl;
      //  std::cout <<"slice "<< mvs << std::endl;
    }
  for(int i = n1; i < n2; i++)
    {
      matrix_vector_slice<rmatrix > mvs (T, slice (i, 0, r), slice (0, 1, r));
      Lbound.push_back(Q(i)-inner_prod(mvs,xvec));
      // cout<<inner_prod(mvs,xvec)<<endl;
      //std::cout <<"slice "<< mvs << std::endl;
    }
  xb.push_back(boost::make_tuple(*(x_lst.rbegin()), Lbound, Ubound));
  if (s+1-n2+n1*(n2-n1) > 0)
    {
      GiNaC::lst in_lst;
      for(int i = n2; i <= s; i++)
	{
	  matrix_vector_slice<rmatrix > mvs (T, slice (i, 0, r), slice (0, 1, r));
	  in_lst.append(Q(i)-inner_prod(mvs,xvec));
	  //  cout<<inner_prod(mvs,xvec)<<endl;
	  // std::cout <<"slice "<< mvs << std::endl;
	}
      for(int i = 0; i < n1; i++)
	for(int ip = n1; ip <n2; ip++)
	  {
	    matrix_vector_slice<rmatrix > mvsi (T, slice (i, 0, r), slice (0, 1, r));
	    matrix_vector_slice<rmatrix > mvsip (T, slice (ip, 0, r), slice (0, 1, r));
	    in_lst.append(Q(i)-Q(ip) + inner_prod(mvsip,xvec) - inner_prod(mvsi,xvec));
	  }
      x_lst.remove_last();
      
      fm_instance(in_lst, x_lst,xb);
    }

}

void fm(lst ineq_lst, lst x_lst, bound_vec&  xb)
{
  using namespace boost::icl;
  fm_instance(ineq_lst, x_lst, xb);

  GiNaC::exhashmap<interval<GiNaC::ex>::type > interval_map;
  lst known_x;
  BOOST_REVERSE_FOREACH(bound bx, xb)
    {
      cout<<bx<<endl;
      exvector Lbound, Ubound;
      GiNaC::exset Ub_min,Lb_max;
      ex x;
      boost::tie(x, Lbound, Ubound) = bx;

      BOOST_FOREACH(ex Ub_ex, Ubound)
	{
	  ex tmp_ex = Ub_ex;
	  for(GiNaC::lst::const_iterator xit = known_x.begin(); xit != known_x.end(); ++xit)
	    {
	      if(Ub_ex.lcoeff(*xit) > 0)
		tmp_ex = tmp_ex.subs((*xit) == interval_map[*xit].lower());
	      else
		tmp_ex = tmp_ex.subs((*xit) == interval_map[*xit].upper());
	      // cout<<"tmpex"<<tmp_ex<<endl;
	    }
	  Ub_min.insert(tmp_ex);
	}

      BOOST_FOREACH(ex Lb_ex, Lbound)
	{
	  ex tmp_ex = Lb_ex;
	  for(GiNaC::lst::const_iterator xit = known_x.begin(); xit != known_x.end(); ++xit)
	    {
	      if(Lb_ex.lcoeff(*xit) > 0)
		tmp_ex = tmp_ex.subs((*xit) == interval_map[*xit].upper());
	      else
		tmp_ex = tmp_ex.subs((*xit) == interval_map[*xit].lower());
	      //      cout<<"tmpex"<<tmp_ex<<endl;
	    }
	  Lb_max.insert(tmp_ex);
	}


      cout<<"ubmin "<<Ub_min<<endl;
      //cout<<"L_bound: "<< *max_element(Lbound.begin(),Lbound.end())<<endl;
      cout<<"U_bound: "<< *min_element(Ub_min.begin(),Ub_min.end())<<endl;

      known_x.append(x);// append x to set with known intervals 
      interval_map[x] = interval<GiNaC::ex>::open(*max_element(Lb_max.begin(),Lb_max.end()), *min_element(Ub_min.begin(),Ub_min.end()));

      cout<<" inter "<<      interval_map[x].upper()<<endl;
      cout<<"Lbound: "<< *max_element(Lbound.begin(),Lbound.end())<<endl;
      cout<<"Ubound: "<< *min_element(Ubound.begin(),Ubound.end())<<endl;
    }

  for(GiNaC::exhashmap<interval<GiNaC::ex>::type >::iterator it = interval_map.begin(); it != interval_map.end(); ++it)
    {
      cout<<it->second.lower()<<" < "<<it->first<<" < "<<it->second.upper()<<endl;
    }

}

GiNaC::exmap findinstance(lst ineq_lst, lst x_lst)
{
  bound_vec bv;
  fm_instance(ineq_lst, x_lst, bv);
  GiNaC::exmap x_subs_map;
  ex den = 2;
  GiNaC::exset point_set;
  BOOST_REVERSE_FOREACH(bound bx, bv)
    {
      exvector Lbound, Ubound;
      GiNaC::exset Ub_min,Lb_max;
      ex x;
      boost::tie(x, Lbound, Ubound) = bx;
      cout<<"X: "<<x<<endl;
      cout<<"Ubo "<<Ubound<<endl;
      cout<<"Lbo "<<Lbound<<endl;

      BOOST_FOREACH(ex Ub_ex, Ubound)
	{
	  Ub_min.insert(Ub_ex.subs(x_subs_map));
	}
      BOOST_FOREACH(ex Lb_ex, Lbound)
	{
	  Lb_max.insert(Lb_ex.subs(x_subs_map));
	}
      ex lb_xi = *max_element(Lb_max.begin(),Lb_max.end());
      ex ub_xi = *min_element(Ub_min.begin(),Ub_min.end());
      if(!(lb_xi < ub_xi)) throw std::logic_error( std::string("Lower bound is not less then Upper bound. solution cn not be found.") );

      ex fraction = pow(den,-1);
      cout<<"point "<< (1-fraction)*lb_xi + fraction*ub_xi<<endl;
      ex point = (1-fraction)*lb_xi + fraction*ub_xi;
      if(point_set.count(point) > 0) // contours do not cross
	do
	  {
	    den += 1;
	    point = (1-pow(den,-1))*lb_xi + pow(den,-1)*ub_xi;
	  }
	while(point_set.count(point) > 0);
      point_set.insert(point);
      x_subs_map[x] = point;
    }
  cout<<"XSUUM   : "<<x_subs_map<<endl;
  return x_subs_map;
}

/*
int main () {
  using namespace boost::numeric::ublas;

  symbol w1("w1"),w2("w2"),w3("w3"),ep("ep");
  bound_vec bv;
  //  fm(lst(1-w1+2*w2,2-7*w3-w2,3-w1+w2-3*w3,4-w1+5*w3,5-w2,6-2*w3,w1+30), lst(w1,w2,w3), bv);
  //  fm(lst(w1,3-w1), lst(w1), bv);
  //  fm(lst(-w1,-w2,-w3,2+ep+w1+w2+w3,1+w3,1+w1+w2+w3,-1-ep-w1-w3,-1-ep-w2-w3), lst(w1,w2,w3,ep), bv);
  findinstance(lst(-w1,-w2,-w3,2+ep+w1+w2+w3,1+w3,1+w1+w2+w3,-1-ep-w1-w3,-1-ep-w2-w3), lst(w1,w2,w3,ep));

  
}

*/
