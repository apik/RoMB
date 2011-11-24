#include "collect_square.h"
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>
#include <cassert>
std::pair<lst,lst> collect_square(ex& F,lst x)
{
  //symbol l("l"),u("u"),g("g"),x1("x1"),x2("x2"),x3("x3"),x4("x4");
  //x = lst(x1,x2,x3,x4);
  //  ex F_out =(3*u*F).expand();
  //  ex F_out = (l*(x1*x1+x2*x2+x3*x3+x4*x4 -2*x1*x4 +2*x2*x4 -2*x3*x4 -2*x1*x2 -2*x2*x3 ) -g*x1*x3 +u*x2*x3).expand();
  ex F_out = F;

  lst x_sq_lst;
  lst x_coeff_lst;
  exset coeff_set;
  std::multimap<ex,ex,ex_is_less> square_map;
  for(lst::const_iterator xit = x.begin(); xit != x.end(); ++xit)
    {
     // cout<<(*xit)<<" -- "<<F.coeff(*xit,2)<<endl;
      if(F.coeff(*xit,2) != 0 )
	{
	  square_map.insert(pair<ex,ex>(F.coeff(*xit,2) ,*xit));
	  coeff_set.insert(F.coeff(*xit,2));
	}
    }
  /*
    iterate over coefficients infront of x_i^2
  */
  BOOST_FOREACH(ex x_coeff_key,coeff_set)
    {
     // cout<<"xinlst: "<<x_coeff_key<<endl;
      if(square_map.count( x_coeff_key ) > 1) // count more then one element with the same coefficient
	{
	  //vector for mapping between X_i and matrix index
	  exvector x_to_index;
	  size_t x=0,y=0;
	  // adjacency matrix for coefficient 
	  matrix xm(square_map.count( x_coeff_key ),square_map.count( x_coeff_key ));
	  // list for zero columns and rows
	  std::deque<size_t> zero_cr;  // deque due to push_fron() operation
	  cout<<"M("<<square_map.count( x_coeff_key )<<"x"<<square_map.count( x_coeff_key )<<")"<<endl;
	  std::multimap<ex,ex,ex_is_less>::iterator xi1,xi2,xi_end;;
	  boost::tie(xi1,xi_end) = square_map.equal_range(x_coeff_key);
	  for(xi1; xi1 != xi_end; ++ xi1){y=0;                 //iterate over all pairs of x_i
	    x_to_index.push_back(xi1->second); // filling vector of X_i
	    for(xi2 = square_map.equal_range(x_coeff_key).first; xi2 != boost::next(xi1); ++ xi2)         //may be next(xi1)???
	      {
		if(F.coeff(xi1->second,1).coeff(xi2->second,1).has(2*(x_coeff_key)))
		  {
		    //cout<<"("<<x<<","<<y<<")"<<endl;
		    xm(x,y) = (xi1->second)*(xi2->second);
		    xm(y,x) = (xi1->second)*(xi2->second);
		  }
		else if(F.coeff(xi1->second,1).coeff(xi2->second,1).has(-2*(x_coeff_key)))
		  {
		    //cout<<"("<<x<<","<<y<<")"<<endl;
		    xm(x,y) = -(xi1->second)*(xi2->second);
		    xm(y,x) = -(xi1->second)*(xi2->second);
		  }
		else if(xi1 == xi2)
		  {
		    xm(x,y) = (xi1->second)*(xi2->second);
		    xm(y,x) = (xi1->second)*(xi2->second);
		  }
		else 
		  {
		    xm(x,y) = 0;
		    xm(y,x) = 0;
		    zero_cr.push_front(x);//removing by X
		  }
		y++;
	      }
	    x++;
	  }
	  /*
	    removing columns and rows with zeros
	  */
	  BOOST_FOREACH(size_t cr_todel, zero_cr)
	    {
	      ex m_ex = reduced_matrix(xm, cr_todel, cr_todel);
	      xm = ex_to<matrix>(m_ex);
	      x_to_index.erase(x_to_index.begin()+cr_todel);
	    }
	  //cout<<"xm after no term : "<<xm<<endl;
	  /*
	    signum matrix s23*s31==s21, if not remove r3,c3
	  */
	  ex rm = xm;
	  bool matrix_reduced = true;
	  while(matrix_reduced)
	    {
	      matrix_reduced = false;
	      //xm = ex_to<matrix>(rm);
	      cout<<"xm "<<xm<<endl;
	      int i=0,j=0,k=0;
	      while(i < xm.rows() && !matrix_reduced)
		{
		  while(k < xm.rows() && !matrix_reduced)
		    {
		      j = i;
		      while(j < xm.rows() && !matrix_reduced)
			{
			  /*
			    Sign conventions
			    (xi*xj).nops()==2 =>sij==(-1)^2
			    (-xi*xj).nops()==3 => sij==(-1)^3
		          
			  */
			  ex s_ij = pow(-1, xm(i,j).nops());
			  ex s_ik = pow(-1, xm(i,k).nops());
			  ex s_kj = pow(-1, xm(k,j).nops());
			  if(s_ij != s_ik*s_kj)
			    {
			      //cout<<k<<" "<<k<<" to be reduced"<<endl;
			      ex rm =  reduced_matrix(xm,k,k);
			      x_to_index.erase(x_to_index.begin()+k);
			      //cout<<"rm: "<<ex_to<matrix>(rm)<<endl;
			      xm = ex_to<matrix>(rm);
			      //assert(false);
			      matrix_reduced = true;
			    }
			  j++;
			}
		      k++;
		    }
		  i++;
		}
	    }
          if(!(xm.rows() < 2))
            {
         
	  cout<<xm<<endl;
	  /*
	    Signum calculation
	    S_1 = 1
	    S_(i+1) = S_(i+1,i)*S_i
	   */
	  ex x_sum_sq = 0;
	  ex Si = 1;
	  for(int j = 0; j < x_to_index.size(); j++)
	    {
	      cout<<x_to_index[j]<<endl;
	      x_sum_sq += Si*x_to_index[j];
	      if((j+1) < x_to_index.size())
		Si *= pow(-1, xm(j+1,j).nops());
	    }
	  x_sq_lst.append(x_sum_sq);
	  x_coeff_lst.append(x_coeff_key);  
	  // Result test part
	  //cout<<endl<<endl<<"The square : "<<pow(x_sum_sq,2)<<endl<<endl;
	  ex out_sq = pow(x_sum_sq,2);
	  exvector unit_vec(xm.rows(),1); // unit vector to fill matrix with units
	  lst unit_lst(unit_vec.begin(),unit_vec.end());
	  lst matrix_lst;
	  for(int rcnt = 0; rcnt < unit_lst.nops(); rcnt++)
	    matrix_lst.append(unit_lst);
	  ex ltm = xm.mul(ex_to<matrix>(lst_to_matrix(matrix_lst))).trace();
	  //cout<<" LTM: "<<ltm<<endl;
	  BOOST_ASSERT_MSG((ltm - out_sq).expand() == 0,"FULL SQUARE EXPANDED INEQUAL");
	  F_out -= (x_coeff_key*out_sq).expand();
	}
    }
    F=F_out;
    }
  return std::make_pair(x_coeff_lst,x_sq_lst);
   
}
