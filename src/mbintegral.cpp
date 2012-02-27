#define GMPRATIONAL
#include <boost/assert.hpp>
#include "mbintegral.h"
#include "utils.h"

#include "constracc.h"
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <setoper.h>

//#include <ppl.hh>
//using namespace Parma_Polyhedra_Library;
//using namespace Parma_Polyhedra_Library::IO_Operators;
#include "ppl_interface.h"
extern "C" {
#include <cdd.h>
}



//#include "tree_util.hh"
/**
 *
 *  Construction from F, U=1
 *
 */

MBintegral::MBintegral(
                       UFXmap fx_in,
                       lst nu,
                       numeric l,
                       bool subs_U, 
                       unsigned int displacement) :tree_level(0),res_pole(0),opt_flag(false) // lst nu is a list of powers of propagators and l is a number of loops
{
  try
    {
      eps_current = (get_symbol("eps")==0);
      ex N = accumulate(nu.begin(),nu.end(),ex(0));


      lst  x_lst(fusion::at_key<UFX::xlst>(fx_in));
     
      ex F = (fusion::at_key<UFX::F>(fx_in)).expand();///< distributed polynom factorizing X(i)X(j)*(...)
      ex U = (fusion::at_key<UFX::U>(fx_in)).expand();
      ex F_pow = (N-l*(2-get_symbol("eps")));
      ex U_pow = N - (l+1)*(2-get_symbol("eps"));
      exset f_set;
      if(F.find(wild(1)*wild(1)+2*wild(1)*wild(2)+wild(2)*wild(2), f_set))cout<< "HAVE square"<<endl;
      cout<<setw(30)<<std::internal<<"+++INTEGRAL PARAMETERS+++"<<endl;
      cout<<setw(30)<<std::left<<"** Number of loops L=   "<<std::internal<<l<<endl;
      cout<<setw(30)<<std::left<<"** Summ of powers N=   "<<N<<endl;
      //cout<<setw(35)<<std::left<<"** Dimension D=   "<<D.subs(D_subs)<<endl;
      cout<<setw(35)<<std::left<<"** F =   "<<F<<endl;
      cout<<setw(35)<<std::left<<"** F power=   "<<F_pow<<endl;
      cout<<setw(30)<<std::internal<<"+++++++++++++++++++++++++"<<endl;

      //	assert(false);
        
      exmap x_power_map; // map of X(j) powers
      for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
        {
          x_power_map[*it1] = nu.op(std::distance(x_lst.begin(),it1))-1;
        }

      cout<<"x_map_start "<<x_power_map<<endl;

		
      ex coeff = 1;                             // numerical coeeficient independent of X(j)
      // ex coeff = pow(exp(get_symbol("eps")*Euler),l)*pow(1,U_pow);///pow(I*pow(Pi,2 - get_symbol("eps")),l);
     // ex Laporta_factor =I/ tgamma(1+get_symbol("eps"));
    //    ex coeff = Laporta_factor;
      // important if power = 0????
      for(lst::const_iterator nui = nu.begin();nui!=nu.end();++nui)
        coeff/=tgamma(*nui);
     
      /************************************************************ 
       |  collecting squares only if it's efficient,              |
       |  number of terms in expanded F="n_0", in collected "n_c" |
       |  and length of square "s"                                |
       |  n_c > n_0 - C(s,2), where C is a binomial coefficient   |
       |  for arbitrary number of squares                         |
       |  n_c > n_o - Summ_j( C(s_j,2) )                          |
       ***********************************************************/        
      size_t n_0 = F.expand().nops();
      size_t n_c = F.collect(x_lst,true).nops();
      // F_col_sq modified, for decision
      ex F_col_sq = F;
      lst coe_l,xsq_l;
      cout<<"test collect squares"<<endl;
      tie(coe_l,xsq_l) = collect_square(F_col_sq,x_lst);
      cout<<">>> Found "<<coe_l.nops()<<" full squares in F polynomial"<<endl;
      cout<< coe_l<<" * "<<xsq_l<<endl;
      cout<<" F qad " <<F_col_sq<<endl;
      n_0  = n_0 - F_col_sq.expand().nops() + F_col_sq.collect(x_lst,true).nops();
      //      ex nFprime = ex_to<add>(F_col_sq.collect(x_lst,true)).nops();
      size_t nFprime = is_a<add>(F_col_sq.collect(x_lst,true)) ? F_col_sq.collect(x_lst,true).nops() : 1;
      //working with F-term \Gamma(\nu-L*D/2) contractedx
      ex w_sum = 0;  //F-term generates only integrations in W


      // decide to collect or not
      // sumCj = sum_j(C(j,2))
      ex sumCj(0);
      ex sumlj(0);
      for(lst::const_iterator cit = xsq_l.begin(); cit != xsq_l.end(); ++cit)
        {
          sumCj += binomial(cit->nops(),2);
          sumlj += std::min((*(cit)+1-U).nops(),(*(cit)-1+U).nops());
        }

      cout<<" n_c = "<<n_c<<" n_0("<<nFprime<<") + sumlj("<<sumlj<<") =  "<< nFprime + sumlj<<endl;
      /*********************************
         COLLECT FULL SQARE PART
      **********************************/
      //      if(sumCj >0 && n_c >= n_0 - sumCj) // collecting squares
      if(sumCj >0 && n_c >= nFprime + sumlj) // collecting squares
        {
          cout<< ">>>  COLLECTING SQUARES!!!!!!!"<<endl;
          F = F_col_sq;
          F = F.collect(x_lst,true);
            
            
            
          //--------------------------------
          //      MB for full squares
          //--------------------------------
          size_t z_idx = displacement;
          if(F.nops() == 0 && xsq_l.nops() == 1)
            {
              /*
                F is a full square
              */
              F_pow *=2;
              F = xsq_l.op(0);
              coeff /= coe_l.op(0);
            }
          else if(F.nops() == 0 && xsq_l.nops() > 1)
            {
              /*
                F is a summ of full squares
              */
              throw std::logic_error(std::string("F = (xs1)^2+ ..(xsn)^2; not realized "));
            }
          else
            {
              for(size_t i = 0; i < xsq_l.nops(); i++)
                {
                  ex sq_lst(xsq_l.op(i)); // summ for MB
                  // substituting U=Summ(x_j)=1
                  if(subs_U && (sq_lst.nops() > std::min((sq_lst + U -1).nops(), (sq_lst - U + 1).nops())))
                    {
                      sq_lst = (sq_lst + U -1).nops() < (sq_lst - U + 1).nops() ? (sq_lst+U-1) : (sq_lst+1-U);
                    }
                  
                  string str = "w"+boost::lexical_cast<string>(displacement);
                  displacement++;
                  //symbol w(str);
                  ex w_i = get_symbol(str);
                  //     w_lst.append(w_i);
                  insert_w(w_i);
                  coeff*=tgamma(-w_i);
                  cout<<"w_i_power "<<w_i<<endl;
                  //                gamma_poles.append(-w_i); //!!!! review
                  insert_pole(-w_i); //!!!! review
                  w_sum+=w_i;
                  cout<<"SQ_LST: "<<sq_lst<<"  .nops()"<<sq_lst.nops()<<endl;                  

                  if(sq_lst.nops() > 1)
                    {
                    //edited 2Pi*I
                      coeff *= pow(coe_l.op(i),w_i);///(2*Pi*I);
                      cout<<"MEGA COEFF: "<<pow(coe_l.op(i),w_i)<<endl;
                      // subMB construction
                      // ordering in full square expr
                      // x_lst -> to list<symbol> for comparision
                      sq_lst = sq_lst.collect(x_lst,true);
                      std::list<symbol> tmp_x_sym_list;
                      for(lst::const_iterator simit = x_lst.begin(); simit != x_lst.end(); ++simit)
                        if(is_a<symbol>(*simit)) tmp_x_sym_list.push_back(ex_to<symbol>(*simit));
                      ex_is_lesseq_degrevlex subF_comp(tmp_x_sym_list);
                      exlist subFl(sq_lst.begin(),sq_lst.end());
                      // sorting lexicographicaly
                      subFl.sort(subF_comp);

                  
                  
                      cout<<"SQLST : "<<subFl<<endl;
                      //edited 2Pi*I
                      coeff /= tgamma(-2*w_i);//(pow(2*Pi*I,sq_lst.nops()-1)*tgamma(-2*w_i));
                      insert_pole(-2*w_i);
                      ex z_sum = 0;
                      for(exlist::iterator x_it = subFl.begin(); x_it != subFl.end(); ++x_it)
                        {
                          ex a_power;
                          if(subFl.end() == boost::next(x_it)) // X_k expr
                            {
                              coeff *= tgamma(-2*w_i + z_sum);
                              a_power = 2*w_i - z_sum;
                              insert_pole(-2*w_i + z_sum);
                            }
                          else // ordinary exprs
                            {
                              string z_str = "z_"+boost::lexical_cast<string>(z_idx);
                              z_idx++;
                              ex z_i = get_symbol(z_str);
                              cout << z_i<<endl;
                              insert_w(z_i);
                              //			w_lst.append(z_i);
                              coeff*=tgamma(-z_i);
                              a_power = z_i;
                              insert_pole(-z_i); //!!!! review
                              z_sum += a_power;
                            } 
		    
                          // add x-part with it's coefficient
                          for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
                            {
                              if(x_it->has(*it1))
                                {
                                  // simple expression x_i or -x_i, x_it->degree(x) == 1
                                  BOOST_ASSERT_MSG(x_it->degree(*it1) == 1,"Not a simple expression in square");
                                  x_power_map[*it1]+=(a_power);
                                  coeff *=pow(x_it->lcoeff(*it1),a_power);
                                }
                            }//it1
                        }//xit
                    }// more then one term in x_sq_lst
                  else // only one term in square
                    {
                      ex tmsq = sq_lst;
                      for(lst::const_iterator it11=x_lst.begin();it11!=x_lst.end();++it11)
                        if(sq_lst.has(*it11))
                          {
                            x_power_map[*it11] += w_i*sq_lst.degree(*it11);
                            tmsq = tmsq.lcoeff(*it11);
                          }
                      coeff *= pow(tmsq,w_i);
                      coeff *= pow(coe_l.op(i),w_i);///(2*Pi*I);
                      cout<<"MEGA COEFF: "<<pow(coe_l.op(i),w_i)<<endl;
                      cout<<"one term constructed"<<endl; 
                    }

                }
            }

    
    
          // x_lst -> to list<symbol> for comparision
          std::list<symbol> x_sym_list;
          for(lst::const_iterator sit = x_lst.begin(); sit != x_lst.end(); ++sit)
            if(is_a<symbol>(*sit)) x_sym_list.push_back(ex_to<symbol>(*sit));
          f_lesseq_revlex F_comp(x_sym_list);
          exlist Fl;
          if(is_a<add>(F))
            Fl.assign(F.begin(),F.end());
          else 
            Fl.push_back(F);
          //   cout<<"FEX "<<Fex<<endl;

          // sorting lexicographicaly
          Fl.sort(F_comp);
          Fl.reverse();
          cout<<"FEX "<<Fl<<endl;



          /*        if(is_a<add>(F))
                    for(const_iterator it = F.begin();it!=F.end();++it)
                    F_to_lst.append(*it);
                    else 
                    F_to_lst.append(F);
       
        
                    //        F_to_lst.sort();
        
                    //        comp_ex_xpow F_term_comparator(x_lst);
                    //std::sort(F_to_lst.begin(),F_to_lst.end(),F_term_comparator);
                    F_to_lst = bubble_sort_lexi(F_to_lst,x_lst);
          */
          cout<<"displace:  "<<displacement<<endl;


          //        for(lst::const_iterator it = F_to_lst.begin();it!=F_to_lst.end();++it)

          // for_each term F_term in sorted Fl
          size_t w_index = 0;
          BOOST_FOREACH(ex F_term, Fl)
            {
              cout<<"F_term : "<<F_term<<endl;
              //            int w_index = std::distance(F_to_lst.begin(),it);
            
              ex x_power;
              //            if(Fl.end()==boost::next(it)) 
              if(Fl.back() == F_term)
                {
                  coeff*=tgamma(F_pow+w_sum);
                  x_power = -F_pow-w_sum;
                  cout<<"x_power last term "<<x_power<<"  w_sum "<<w_sum<<endl;
                  insert_pole(F_pow+w_sum);
                  cout<<"end achieved"<<endl;
                }
              else
                {
                  string str = "w"+boost::lexical_cast<string>(displacement + w_index);
                  //symbol w(str);
                  //                w_lst.append(get_symbol(str));
                  insert_w(get_symbol(str));
                  coeff*=tgamma(-get_symbol(str));
                  x_power = get_symbol(str);
                  cout<<"x_power "<<x_power<<endl;
                  insert_pole(-x_power); //!!!! review
                  w_sum+=x_power;
                  cout<<"ok run"<<endl;
                }
              // filling map of X(j) powers
              ex tmp_expr = F_term;
              for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
                {
                  if(tmp_expr.has(*it1))
                    {
                      x_power_map[*it1]+=(tmp_expr.degree(*it1)*x_power);
                      cout<<"before subs "<<(*it1)<<"   "<<tmp_expr<<endl;
                      tmp_expr = tmp_expr.subs((*it1)==1);
                      //                    cout<<"after subs "<<tmp_expr<<endl;
                    }
                }
              //     cout<<"after subs "<<tmp_expr<<endl;
              coeff *=pow(tmp_expr,x_power);
              //increment W index
              w_index++;
            }
    
          //working with U-term, only if (U_pow>0)
        
       
          //cout<<x_power_map<<endl;
          //        cout<<"Gammas after MB: "<<endl<<get_gamma_poles()<<endl;
          cout<<"X powers list:"<<"  "<<x_power_map<<endl;
          //!!!!!!!!!!!!!!!!!!
          //        assert(false);
          // applying X integration
          ex gamma_den = 0; // gamma in denominator
          for(exmap::const_iterator mi = x_power_map.begin();mi!=x_power_map.end();++mi)
            {
              cout<<(*mi).first<<" "<<(*mi).second<<endl;
              insert_pole((*mi).second+1);
              coeff*=tgamma((*mi).second+1);
              gamma_den+=((*mi).second+1);
            }
          cout<<"GAMMA_DEN: "<<gamma_den<<endl;
          /*
            bool gamma_den_has_w = false;

            for(lst::const_iterator wit = w_lst.begin(); wit != w_lst.end(); ++wit)
            if(gamma_den.has(*wit))gamma_den_has_w = true;
            if(gamma_den_has_w) 
          */  
          //      insert_pole(gamma_den);

          coeff/=tgamma(gamma_den);
          coeff/=(pow(2*Pi*I,w_lst.size()));
          //        cout<<"New gamma list:"<<endl<<get_gamma_poles()<<endl;
          full_int_expr = coeff;
          cout<<w_lst<<endl;
        }

      /*********************************
         NO COLLECT FULL SQARE PART
      **********************************/

      else // no collect squares
        {
          F = F.collect(x_lst,true);
          // x_lst -> to list<symbol> for comparision
          std::list<symbol> x_sym_list;
          for(lst::const_iterator sit = x_lst.begin(); sit != x_lst.end(); ++sit)
            if(is_a<symbol>(*sit)) x_sym_list.push_back(ex_to<symbol>(*sit));
       
          //          ex_is_lesseq_degrevlex F_comp(x_sym_list);
          f_lesseq_revlex F_comp(x_sym_list);
          exlist Fl;
          if(is_a<add>(F))
            Fl.assign(F.begin(),F.end());
          else 
            Fl.push_back(F);

          //   cout<<"FEX "<<Fex<<endl;

          // sorting lexicographicaly
          Fl.sort(F_comp);
          Fl.reverse();
          cout<<"FEX "<<Fl<<endl;



          /*        if(is_a<add>(F))
                    for(const_iterator it = F.begin();it!=F.end();++it)
                    F_to_lst.append(*it);
                    else 
                    F_to_lst.append(F);
       
        
                    //        F_to_lst.sort();
        
                    //        comp_ex_xpow F_term_comparator(x_lst);
                    //std::sort(F_to_lst.begin(),F_to_lst.end(),F_term_comparator);
                    F_to_lst = bubble_sort_lexi(F_to_lst,x_lst);
          */
          cout<<"displace:  "<<displacement<<endl;


          //        for(lst::const_iterator it = F_to_lst.begin();it!=F_to_lst.end();++it)

          // for_each term F_term in sorted Fl
          size_t w_index = 0;
          BOOST_FOREACH(ex F_term, Fl)
            {
              cout<<"F_term : "<<F_term<<endl;
              //            int w_index = std::distance(F_to_lst.begin(),it);
            
              ex x_power;
              //            if(Fl.end()==boost::next(it)) 
              if(Fl.back() == F_term)
                {
                  coeff*=tgamma(F_pow+w_sum);
                  x_power = -F_pow-w_sum;
                  cout<<"x_power last term "<<x_power<<"  w_sum "<<w_sum<<endl;
                  insert_pole(F_pow+w_sum);
                  cout<<"end achieved"<<endl;
                }
              else
                {
                  string str = "w"+boost::lexical_cast<string>(displacement + w_index);
                  //symbol w(str);
                  //                w_lst.append(get_symbol(str));
                  insert_w(get_symbol(str));
                  coeff*=tgamma(-get_symbol(str));
                  x_power = get_symbol(str);
                  cout<<"x_power "<<x_power<<endl;
                  insert_pole(-x_power); //!!!! review
                  w_sum+=x_power;
                  cout<<"ok run"<<endl;
                }
              // filling map of X(j) powers
              ex tmp_expr = F_term;
              for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
                {
                  if(tmp_expr.has(*it1))
                    {
                      x_power_map[*it1]+=(tmp_expr.degree(*it1)*x_power);
                      cout<<"before subs "<<(*it1)<<"   "<<tmp_expr<<endl;
                      tmp_expr = tmp_expr.subs((*it1)==1);
                      //                    cout<<"after subs "<<tmp_expr<<endl;
                    }
                }
              //     cout<<"after subs "<<tmp_expr<<endl;
              coeff *=pow(tmp_expr,x_power);
              //increment W index
              w_index++;
            }
    
          //working with U-term, only if (U_pow>0)
        
       
          //cout<<x_power_map<<endl;
          //        cout<<"Gammas after MB: "<<endl<<get_gamma_poles()<<endl;
          cout<<"X powers list:"<<"  "<<x_power_map<<endl;
          //!!!!!!!!!!!!!!!!!!
          //        assert(false);
          // applying X integration
          ex gamma_den = 0; // gamma in denominator
          for(exmap::const_iterator mi = x_power_map.begin();mi!=x_power_map.end();++mi)
            {
              cout<<(*mi).first<<" "<<(*mi).second<<endl;
              insert_pole((*mi).second+1);
              coeff*=tgamma((*mi).second+1);
              gamma_den+=((*mi).second+1);
            }
          cout<<"GAMMA_DEN: "<<gamma_den<<endl;
          /*
            bool gamma_den_has_w = false;

            for(lst::const_iterator wit = w_lst.begin(); wit != w_lst.end(); ++wit)
            if(gamma_den.has(*wit))gamma_den_has_w = true;
            if(gamma_den_has_w) 
          */  
          //      insert_pole(gamma_den);

          coeff/=tgamma(gamma_den);
          coeff/=(pow(2*Pi*I,w_lst.size()));
          //        cout<<"New gamma list:"<<endl<<get_gamma_poles()<<endl;
          full_int_expr = coeff;
          cout<<w_lst<<endl;

        }// no collect squares part           
      
      /*
       *
       *  checking class invariants
       *  
       * 1) w_lst 
       * 2) poles, sorted lexi
       *
       */
      p_lst_type poles_before_inv(gamma_poles);
      gamma_poles.clear();
      ex expr_numer = full_int_expr.numer();
      BOOST_FOREACH(ex epe, poles_before_inv)
        {
          if(expr_numer.has(tgamma(epe))) gamma_poles.push_back(epe);
          else if(expr_numer.has(psi(epe))) gamma_poles.push_back(epe);
          else if(expr_numer.has(psi(wild(),epe))) gamma_poles.push_back(epe);
          else cout<<" UNUSED POLE:"<<epe<<endl<<endl;
        }
      w_lst_type w_before_inv(w_lst);
      w_lst.clear();
      // cerating comparator
      std::list<symbol> w_lst_for_comparision;
      BOOST_FOREACH(ex ewe, w_before_inv)
        {
          if(full_int_expr.has(ewe))
            {
              w_lst.push_back(ewe);
              if(is_a<symbol>(ewe)) w_lst_for_comparision.push_back(ex_to<symbol>(ewe));
            }
          else cout<<" UNUSED W:"<<ewe<<endl<<endl;
        }
      w_lst_for_comparision.push_back(get_symbol("eps"));
      comparator = new ex_is_lesseq_degrevlex(w_lst_for_comparision);
      w_lst.sort(*comparator);
      gamma_poles.sort(*comparator);
      //gamma_poles = poles_from_ex(full_int_expr);      
   

    }catch(std::exception &p)
    {
      throw std::logic_error(std::string("In function \"MBintegral-UFXversion\":\n |___> ")+p.what());
    }

}
void MBintegral::fix_inv()
{
    p_lst_type poles_before_inv(gamma_poles);
    gamma_poles.clear();
    ex expr_numer = full_int_expr.numer();
    BOOST_FOREACH(ex epe, poles_before_inv)
    {
        if(expr_numer.has(tgamma(epe))) gamma_poles.push_back(epe);
        else if(expr_numer.has(psi(epe))) gamma_poles.push_back(epe);
        else if(expr_numer.has(psi(wild(),epe))) gamma_poles.push_back(epe);
        else cout<<" UNUSED POLE:"<<epe<<endl<<endl;
    }
    w_lst_type w_before_inv(w_lst);
    w_lst.clear();
    // cerating comparator
    std::list<symbol> w_lst_for_comparision;
    BOOST_FOREACH(ex ewe, w_before_inv)
    {
        if(full_int_expr.has(ewe))
            {
              w_lst.push_back(ewe);
              if(is_a<symbol>(ewe)) w_lst_for_comparision.push_back(ex_to<symbol>(ewe));
            }
          else cout<<" UNUSED W:"<<ewe<<endl<<endl;
        }
      w_lst_for_comparision.push_back(get_symbol("eps"));
      comparator = new ex_is_lesseq_degrevlex(w_lst_for_comparision);
      w_lst.unique();
      w_lst.sort(*comparator);
      gamma_poles.unique();
      gamma_poles.sort(*comparator);
}
/*
  lst MBintegral::get_pole_lst()
  {

  cout<<"get_pole_lst called and ret "<<get_gamma_poles()<<endl;

  //update_poles_from_ex();
  exset gammaset,psiset,psi2set;

  full_int_expr.find(tgamma(wild()),gammaset);
  full_int_expr.find(psi(wild()),psiset);
  full_int_expr.find(psi(wild(),wild(1)),psi2set);
  cout<<" but must gamma "<<gammaset<<endl;
  cout<<" but must psi(ex) "<<psiset<<endl;
  cout<<" but must psi(int,ex) "<<psi2set<<endl;
    
  return get_gamma_poles();
  }
*/


MBintegral MBintegral::res(relational w_relation,ex pole,relational new_eps)
{
  try
    {
      // remove W-residue from w-list
      exvector cut_w_vec(w_lst.size()-1);
      std::remove_copy(w_lst.begin(),w_lst.end(),cut_w_vec.begin(),w_relation.lhs());
      /*
        REsidue by Laurent series !!!!!
      */
      ex res_loran = full_int_expr.series(w_relation,0).coeff(w_relation.lhs(),-1);
      /*
        ex new_no_gamma_part = (full_int_expr.subs(tgamma(pole)==pow(-1,-pole.subs(w_relation))/factorial(-pole.subs(w_relation)))).subs(w_relation);
        // new_no_gamma_part  = pow(-1,pole.subs(w_relation))/factorial(pole.subs(w_relation))*full_int_expr.subs(w_relation)
        cout<< new_no_gamma_part<<endl;
      */
      exmap new_w_current(w_current);
      //   cout<<" Not modif:  "<<new_w_current<<endl;

      new_w_current.erase(w_relation.lhs());
      MBintegral resINT(lst(cut_w_vec.begin(),cut_w_vec.end()),res_loran,new_w_current,new_eps,tree_level+1);
      resINT.set_respole(pole); // save gamma argument, generated residue
      return resINT;
    }catch(std::exception &p)
    {
      throw std::logic_error(std::string("In function \"MBintegral.res\":\n |___> ")+p.what());
    }
}

lst MBintegral::poles_with_eps()
{
  lst lst_with_eps;
  for(w_lst_type::const_iterator it = gamma_poles.begin(); it != gamma_poles.end(); ++it)
    if(it->has(get_symbol("eps"))) lst_with_eps.append(*it);
  lst_with_eps.unique();
  return lst_with_eps;
}
/*
  lst MBintegral::get_gamma_poles()
  {
  update_poles_from_ex();
  return set2lst(gamma_poles);
  }
*/

lst MBintegral::has_w(const ex& gamma_arg,w_lst_type w_list)
{
  lst out_list;
  for(w_lst_type::const_iterator wi = w_list.begin();wi!=w_list.end();++wi)
    if(gamma_arg.has(*wi))out_list.append(*wi);
  return out_list;
}

lst MBintegral::has_w(const ex& gamma_arg)
{
  lst out_list;
  for(w_lst_type::const_iterator wi = w_lst.begin();wi!=w_lst.end();++wi)
    if(gamma_arg.has(*wi))out_list.append(*wi);
  return out_list;
}


/*exmap sp_cdd(MBintegral::w_lst_type pole_list,MBintegral::p_lst_type w_list)
{
}
*/
dd_MatrixPtr dd_ConstrList2Matrix (MBintegral::w_lst_type pole_list,MBintegral::p_lst_type w_list)//, dd_ErrorType *Error)
{
  dd_MatrixPtr M=NULL;
  dd_rowrange m_input,i;
  dd_colrange d_input,j;
  dd_RepresentationType rep=dd_Inequality;
  mytype value;
  dd_boolean found=dd_FALSE, newformat=dd_FALSE, successful=dd_FALSE, linearity=dd_FALSE;
  //  char command[dd_linelenmax], comsave[dd_linelenmax], numbtype[dd_wordlenmax];
  dd_NumberType NT;
#if !defined(GMPRATIONAL)
  double rvalue;
#endif

  dd_init(value);
  // (*Error)=dd_NoError;

  //  M=dd_CreateMatrix(m_input, d_input);
  M=dd_CreateMatrix(pole_list.size(), w_list.size()+1);
  M->representation = dd_Inequality;
  M->numbtype = dd_Rational;

  i = 0;
  //j = 1;
  signed long a_num;
  unsigned long b_den;
     
  for(MBintegral::p_lst_type::const_iterator pit = pole_list.begin();pit!= pole_list.end();++pit) 
    {
      ex b = *pit;
      j = 1;
      for(MBintegral::w_lst_type::const_iterator wi = w_list.begin();wi != w_list.end();++wi) 
        {
          a_num = numer(ex_to<numeric>(pit->coeff(*wi,1))).to_long();
          b_den = denom(ex_to<numeric>(pit->coeff(*wi,1))).to_long();
          //cout<<"VAL = "<<a_num<<"/"<<b_den<<endl;

          dd_set_si2(value,a_num,b_den);
          // dd_WriteNumber(stderr, value);
          dd_set(M->matrix[i][j],value);
          b = b.coeff(*wi,0);
          j++;
        }      
      a_num = numer(ex_to<numeric>(b)).to_long(); 
      b_den = denom(ex_to<numeric>(b)).to_long();
      dd_set_si2(value,a_num,b_den);
      dd_set(M->matrix[i][0],value);
      i++;
    }

  dd_clear(value);
  /* if (f!=NULL) fclose(f); */
    return M;
}


exmap MBintegral::start_point_diff_w(MBintegral::w_lst_type pole_list,MBintegral::p_lst_type w_list)
{
  try{

  /* 
     test part x_0= (A'*A)^(-1)*A'*b 
  
  matrix A(pole_list.nops(),w_list.nops());
  matrix B(pole_list.nops(),1);
  size_t i,j;
  i = 0;
  j = 0;
for(lst::const_iterator it = pole_list.begin(); it != pole_list.end(); ++it)
    {
      ex b = *it;
      i = 0;
      for(lst::const_iterator xit = w_list.begin(); xit != w_list.end(); ++xit)
	{
	  A(j,i) = (it->coeff(*xit,1));  // (row,col) 
	  i++;
	  b = b.coeff(*xit,0);
	}
      B(j,0) = b;
      j++;
    }
 matrix x0 =  ex_to<matrix>( (pow(A.transpose()*A,-1)*A.transpose()*B).evalm());
 cout<<"AB: "<<(pow(A.transpose()*A,-1)*A.transpose()*B).evalm()<<endl;
 exmap x0_w;
 for(size_t ct = 0; ct < x0.rows();ct++)
   x0_w[w_list.op(ct)] = x0(ct,0);
 cout<<x0_w<<"   int point  "<<interior_point(pole_list,x0_w)<<endl;
  // end test
*/
    lst constraints, constraints_wo_eps;
    lst w_con,w_con_eps0;
    BOOST_FOREACH(ex ce,pole_list)
      {
        constraints.append(ce);
        constraints_wo_eps.append(ce.subs(get_symbol("eps") == 0));
      }
    BOOST_FOREACH(ex we,w_list)
      {
        w_con.append(we);
        if(we != get_symbol("eps"))
          w_con_eps0.append(we);
      }
    cout<<constraints<<endl;
    // BOOST_ASSERT_MSG(!zero_volume(constraints,w_con)," ZERO VOLUME AT START");
    exmap subs_map;

    /********************************
     *               EPS = 0 part                  *
     *******************************/
    if(false && !zero_volume(constraints_wo_eps,w_con_eps0)) 
      {
        cout<<"!!! NO CONTINUATION eps=0"<<endl;
        subs_map[get_symbol("eps")] = 0;
        constraints = constraints_wo_eps;
        exset point_set;
        ex den = 2.75;
        for(MBintegral::w_lst_type::const_reverse_iterator wi = w_list.rbegin();wi != w_list.rend();++wi) 
          {
            lst tmp_pole_list;
            lst tmp_w_list;
            for(MBintegral::p_lst_type::const_iterator pit = pole_list.begin();pit!= pole_list.end();++pit)
              if(!((*pit).subs(subs_map)>0))
                tmp_pole_list.append((*pit).subs(subs_map));
            cout<<"tmp_pole_list"<<tmp_pole_list<<endl;
            for(MBintegral::w_lst_type::const_reverse_iterator wi2 = wi;wi2 != w_list.rend();++wi2) 
              if((*wi2 != get_symbol("eps")) && (subs_map.count(*wi2) ==0) )tmp_w_list.append(*wi2);
            if(subs_map.count(*wi) ==0)
              {
            //      std::pair<ex,double> ret_pair = 
            //  cout<<" Hyper test   "<<hyper_cube_den(tmp_pole_list,tmp_w_list,den).second;
            cout<<"Simplex called for:"<<endl;
            cout<<"w set: "<<tmp_w_list<<endl;
            cout<<"Poles set: "<<tmp_pole_list<<endl;
            //simplex_zero(tmp_pole_list,tmp_w_list);
            // assert(false);
            std::pair<ex,ex> ret_pair;
              ret_pair = hyper_cube_den(tmp_pole_list,tmp_w_list,3.5);

            subs_map[ret_pair.first] = ret_pair.second;
            point_set.insert(ex_to<numeric>(ret_pair.second).to_double());
            cout<<ret_pair.first<<" "<<ret_pair.second<<endl;
              }
          }
        
      }

    /********************************
     *               EPS != 0 part                  *
     *******************************/
    else
      {
        exset point_set;
        ex den = 2.75;
        for(MBintegral::w_lst_type::const_reverse_iterator wi = w_list.rbegin();wi != w_list.rend();++wi) 
          {
            lst tmp_pole_list;
            lst tmp_w_list;
            for(MBintegral::p_lst_type::const_iterator pit = pole_list.begin();pit!= pole_list.end();++pit)
              if(!((*pit).subs(subs_map)>0))
                tmp_pole_list.append((*pit).subs(subs_map));
            cout<<"tmp_pole_list"<<tmp_pole_list<<endl;
            for(MBintegral::w_lst_type::const_reverse_iterator wi2 = wi;wi2 != w_list.rend();++wi2) 
              tmp_w_list.append(*wi2);
            
            //      std::pair<ex,double> ret_pair = 
            //  cout<<" Hyper test   "<<hyper_cube_den(tmp_pole_list,tmp_w_list,den).second;
            cout<<"Simplex called for:"<<endl;
            cout<<"w set: "<<tmp_w_list<<endl;
            cout<<"Poles set: "<<tmp_pole_list<<endl;
            //simplex_zero(tmp_pole_list,tmp_w_list);
            // assert(false);
            std::pair<ex,ex> ret_pair;
            if(wi == w_list.rbegin())  //  epsilon minimum
              ret_pair = hyper_cube_den(tmp_pole_list,tmp_w_list,den);
            else  
              ret_pair = hyper_cube_den(tmp_pole_list,tmp_w_list,den);

            if(point_set.count(ret_pair.second) > 0)
              do
                {
                  den +=1;
                  ret_pair = hyper_cube_den(tmp_pole_list,tmp_w_list,den);
                }
              while(point_set.count(ret_pair.second) > 0);
            subs_map[ret_pair.first] = ret_pair.second;
            point_set.insert(ret_pair.second);
            cout<<ret_pair.first<<" "<<ret_pair.second<<endl;
          }
      }
    cout<<"START POINT SUBS  "<<subs_map<<endl;
    return subs_map;
  }catch(std::exception &p)
    {
      throw std::logic_error(std::string("In function \"start_point_diff_w\":\n |___> ")+p.what());
    }
}


exset MBintegral::poles_from_ex(ex ie_)
{
  ex ie = normal(ie_).numer();
  exset gammaset,psiset,psi2set;
  ie.find(tgamma(wild()),gammaset);
  ie.find(psi(wild()),psiset);
  ie.find(psi(wild(),wild(1)),psi2set);

  exset poles_set;
  BOOST_FOREACH(ex gex,gammaset)
    {
      exmap repls;
      if(gex.match(tgamma(wild()), repls))
        {
          ex gpole = repls[wild()];
          if(has_w(gpole,w_lst).nops() > 0) poles_set.insert(gpole);
        }
    }
  BOOST_FOREACH(ex psiex,psiset)
    {
      exmap repls;
      if(psiex.match(psi(wild()), repls))
        {
          ex psipole = repls[wild()];
          if(has_w(psipole,w_lst).nops() > 0) poles_set.insert(psipole);
        }
    }
  BOOST_FOREACH(ex psi2ex,psi2set)
    {
      exmap repls;
      if(psi2ex.match(psi(wild(),wild(1)), repls))
        {
          ex psipole = repls[wild(1)];
          if(has_w(psipole,w_lst).nops() > 0) poles_set.insert(psipole);
        }
    }    
  return poles_set;
}
/*
  void MBintegral::update_poles_from_ex()
  {
  exset gammaset,psiset,psi2set;
  full_int_expr.find(tgamma(wild()),gammaset);
  full_int_expr.find(psi(wild()),psiset);
  full_int_expr.find(psi(wild(),wild(1)),psi2set);

  exset poles_set;
  BOOST_FOREACH(ex gex,gammaset)
  {
  exmap repls;
  if(gex.match(tgamma(wild()), repls))
  {
  ex gpole = repls[wild()];
  if(has_w(gpole,w_lst).nops() > 0) poles_set.insert(gpole);
  }
  }
  BOOST_FOREACH(ex psiex,psiset)
  {
  exmap repls;
  if(psiex.match(psi(wild()), repls))
  {
  ex psipole = repls[wild()];
  if(has_w(psipole,w_lst).nops() > 0) poles_set.insert(psipole);
  }
  }
  BOOST_FOREACH(ex psi2ex,psi2set)
  {
  exmap repls;
  if(psi2ex.match(psi(wild(),wild(1)), repls))
  {
  ex psipole = repls[wild(1)];
  if(has_w(psipole,w_lst).nops() > 0) poles_set.insert(psipole);
  }
  }    
  gamma_poles = poles_set;
  cout<<"Update poles "<<poles_set<<endl;
  }
*/
/*
exmap chebyshevSphere(MBintegral::w_lst_type wIn, MBintegral::p_lst_type pIn)
{

    std::vector<Variable> varVector;
    size_t varN = 0;
    BOOST_FOREACH(ex w, wIn)
    {
        Variable vNew(varN);
        varVector.push_back(vNew);
        varN++;
    }

    
// Sphere radius

    Variable r(varN);

    Constraint_System cs;


// Add constraints to the system    

    BOOST_FOREACH(ex p, pIn)
    {
        Linear_Expression l;
        size_t dimN = 0;
        ex bI = p;
        ex aNormSquared;
        BOOST_FOREACH(ex w, wIn)
        {
            ex linCoeff = p.coeff(w,1);
            bI = bI.coeff(w,0);
            aNormSquared += pow(linCoeff,2);
            if(linCoeff.info(info_flags::integer))
            {
                mpz_class gmpInt(ex_to<numeric>(linCoeff).to_long());
                l = sub_mul_assign(l, gmpInt, varVector[dimN]);
                dimN++;
            }
            else throw std::logic_error(std::string( "Not an integer coefficient in Linear_Expression" ));
        }
// Add B_i        
        if(bI.info(info_flags::integer))
        {
            mpz_class gmpInt(ex_to<numeric>(bI).to_long());
            l -= gmpInt;
        }
        else throw std::logic_error(std::string( "Not an integer B_i in Linear_Expression" ));

// Add radius
        
        GMP_Integer aNormSqrt;
        if(aNormSquared.info(info_flags::integer))
        {

            mpz_class gmpInt(ex_to<numeric>(aNormSquared).to_long());

// Assign integer value of sqare root: 
//
//      int(sqrt) < sqrt
//
            sqrt_assign(aNormSqrt, gmpInt);
            aNormSqrt += 1;
            cout << "SQRT: " << aNormSqrt << endl;
        }
        else throw std::logic_error(std::string( "Not an integer B_i in Linear_Expression" ));

        l = add_mul_assign(l, aNormSqrt, r);
        
        cout << l << endl;
        cs.insert(l <= 0);
    }

    NNC_Polyhedron ph(cs);
    bool maxVal;
    Coefficient supN,supD;
    Generator g  = closure_point();
    if( ph.maximize(r,supN,supD,maxVal,g) )
    {
        std::cout<< maxVal<<" then " << supN<<"/"<<supD<< "Generator: " << g << std::endl;
       
    }

    ex rEx = numeric(supN.get_si(),supD.get_si());
    Generator_System gs;
    gs.insert(g);
    NNC_Polyhedron ph2(gs);
    if(ph.strictly_contains(ph2)) cout<< " Contains" <<endl;
    
    
// 168 first prime numbers
    unsigned int primeNumbers[168] = {
// no first 28
//        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
 109, 113, 127, 131, 137, 
        139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 
        283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 
        457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 
        631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 
        821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997 
    };
      
    unsigned int* primeNumbersEnd = primeNumbers + 167 - 28;


    cout<< " dimension " << g.space_dimension() << endl;
    
// no r variable
    long int largestNumerator = 0;
    for(dimension_type i = 0; i < g.space_dimension() - 2; i++)
        if(largestNumerator < abs(g.coefficient(varVector[i]).get_si()))
            largestNumerator = abs(g.coefficient(varVector[i]).get_si());
    
    cout << largestNumerator << endl;        
    unsigned int* primeToDivide = upper_bound (primeNumbers, primeNumbersEnd, largestNumerator);

    if (primeToDivide == primeNumbersEnd)
        throw std::logic_error(string("Unavalable prime number to divide"));

    exmap intPtest;
    size_t cncn = 0;
    BOOST_FOREACH(ex w, wIn)
    {
        intPtest[w] = ex_to<numeric>(numeric(g.coefficient(varVector[cncn]).get_si(),g.divisor().get_si()) + rEx/(*primeToDivide)).to_double();
        primeToDivide++;
        cncn++;
    }
        
//       if (binary_search (primeNumbers, primeNumbersEnd, g.coefficient(varVector[i]).get_si()))

//        cout << "found!\n"; else cout << "not found.\n";
    
    cout<< intPtest <<endl;
    
     BOOST_ASSERT_MSG(interior_point(pIn,intPtest),"Not a convex polyhedron interior point");


     return intPtest;
}

*/
bool compExmapSecond(exmap::value_type a, exmap::value_type b)
{
    return (a.second < b.second);
}




class CompExmapFirst : public std::binary_function<numeric,numeric,bool>
{
public:
    bool operator()(const numeric& lh, const numeric& rh)
        {
//            if(is_a<numeric>(lh) && is_a<numeric>(rh))
                return (lh < rh);
                //           else
                //  return lh.compare(rh);
        }
};

template <typename T>
bool alwaysTrue (T i) { return true; }


exmap MBintegral::newPoint()
{
    w_lst_type var_list(w_lst);
    var_list.push_back(get_symbol("eps"));
    exmap eps_w_current = chebyshevSphere(var_list,gamma_poles);

    BOOST_ASSERT_MSG(interior_point(gamma_poles,eps_w_current),"Not a convex polyhedron interior point");


    for(w_lst_type::const_iterator it = w_lst.begin();it!=w_lst.end();++it)
        w_current[*it] = it->subs(eps_w_current);
    eps_current = (get_symbol("eps")==eps_w_current[get_symbol("eps")]);

    return eps_w_current;
    
}




//exmap MBintegral::new_point()

// exmap npppp()   // Absolet function
// {
//     w_lst_type var_list(w_lst);
//     p_lst_type localPoles(gamma_poles);
//     var_list.push_back(get_symbol("eps"));
    
//     // eps_w_current=start_point_diff_w(gamma_poles,var_list);
    
//     exmap wch = chebyshevSphere(var_list,localPoles);
//     eps_w_current = wch;

//     //                    start point from CDDLIB package
//     //  sp_cdd(localPoles,var_list);
//     dd_set_global_constants();  /* First, this must be called. */
//     dd_MatrixPtr mpmp =   dd_ConstrList2Matrix(localPoles,var_list);
//     dd_WriteMatrix(stderr, mpmp);
//     dd_ErrorType err;
//     /*
//       dd_PolyhedraPtr phpp = dd_DDMatrix2Poly(mpmp,&err);
//       dd_MatrixPtr mp34 =  dd_CopyGenerators(phpp);
//       dd_WriteMatrix(stderr, mp34);
//       phpp = dd_DDMatrix2Poly(mp34,&err);
//       mpmp =  dd_CopyInequalities(phpp);
//       dd_WriteMatrix(stderr, mpmp);*/
//     cout<<endl<<"CAnonicalized"<<endl;
//     dd_MatrixPtr dn_cmp;
//     dd_rowset dn_rs;
//     dd_rowset dn_rs2;
//     dd_rowindex dn_ri;
//     // dd_ErrorType *
//     dd_boolean boo  = dd_MatrixCanonicalize(&mpmp, &dn_rs,&dn_rs2, &dn_ri, &err);
//     dd_WriteMatrix(stderr, mpmp);
//     dd_LPSolutionPtr lps;
//     dd_boolean boo2 =  dd_FindRelativeInterior(mpmp, &dn_rs,&dn_rs2, &lps, &err);

//     exmap wCdd;
    

//     typedef multimap<numeric,ex,CompExmapFirst> exmultimap;  
//     exmultimap wInverted;
    
//     exmap wFinal;
//     exmap wFinalInverted;

//     int cn_so = 1;

// // First iteration of start point finding
//     for(w_lst_type::const_iterator it = var_list.begin();it!=var_list.end();++it)
//     {
        
//         if(*it == get_symbol("eps"))
//             wFinal[get_symbol("eps")] = numeric(mpz_get_si(mpq_numref(lps->sol[cn_so])), mpz_get_si(mpq_denref(lps->sol[cn_so])));

//         else
//         {
//             dd_WriteNumber(stderr, lps->sol[cn_so]);
            
//             wCdd[*it] = numeric(mpz_get_si(mpq_numref(lps->sol[cn_so])), mpz_get_si(mpq_denref(lps->sol[cn_so])));
            
//             wInverted.insert(pair<numeric,ex>(ex_to<numeric>(wCdd[*it]), (*it)));
//         }
        
//         cn_so++;
//     }

//     cout<< "\nINVERTED: " << wInverted << endl;
// //    cout<< *wInverted.lower_bound() << endl;

// // inverted map, minimal element is the first
// //    exmap::iterator minElemIter = min_element(wCdd.begin(),wCdd.end(),compExmapSecond);
//     exmultimap::iterator minElemIter = wInverted.begin();

// //    cout << "\nMin: " << minElemIter->first << endl;
// //    wFinalInverted[minElemIter->first] = minElemIter->second;
//     //wCdd.erase(minElemIter);
// //   minElemIter = min_element(wCdd.begin(),wCdd.end(),compExmapSecond);
    
//     while (wInverted.size() > 0)
//     {

// //        minElemIter = min_element(wCdd.begin(),wCdd.end(),compExmapSecond);
//         minElemIter = wInverted.begin();
//         cout<<"tmp\n";        

//         if (count_if (wInverted.equal_range(minElemIter->first).first,wInverted.equal_range(minElemIter->first).second,alwaysTrue<exmultimap::value_type>) == 1)
//         {
//             wFinalInverted[wInverted.begin()->first] = wInverted.begin()->second;
//             wFinal[wInverted.begin()->second] = wInverted.begin()->first;
// /*
//             BOOST_FOREACH(ex lp, localPoles)
//             {
//                 lp.subs(wFinal);
//             }
// */
//             wInverted.erase(minElemIter);
//             cout << " W_F: " << wFinalInverted << endl;
//         }
//         else // Contour duplication
//         {
//         cout<<"Dupl tmp\n";        
// //     Append inequalities : (w - w' >= 0)
//             exmultimap::const_iterator wIter,wIter_end;

//             for(tie(wIter,wIter_end) = wInverted.equal_range(wInverted.begin()->first); wIter != wIter_end; ++wIter)
//             {
// // Add first element 
//                 if(wIter != wInverted.begin())
//                 {
//                     localPoles.push_back(wIter->second - wIter->first);
//                     cout<< " Contour: " << *localPoles.rbegin() << endl;
//                 }                
//             }
//             mpmp =   dd_ConstrList2Matrix(localPoles,var_list);
//             // dd_WriteMatrix(stderr, mpmp);

//             boo  = dd_MatrixCanonicalize(&mpmp, &dn_rs,&dn_rs2, &dn_ri, &err);
// //          dd_WriteMatrix(stderr, mpmp);
            
//             boo2 =  dd_FindRelativeInterior(mpmp, &dn_rs,&dn_rs2, &lps, &err);

//             BOOST_ASSERT_MSG(err == dd_NoError,"Interior not feasible");
//             wFinalInverted[wInverted.begin()->first] = wInverted.begin()->second;
//             wFinal[wInverted.begin()->second] = wInverted.begin()->first;

// /*
//             BOOST_FOREACH(ex lp, localPoles)
//             {
//                 lp.subs(wFinal);
//             }
// */
//             wInverted.clear();
            
//             cn_so = 1;
//             for(w_lst_type::const_iterator it = w_lst.begin();it!=w_lst.end();++it)
//             {
//                 //   dd_WriteNumber(stderr, lps->sol[cn_so]);

//                 if ( wFinal.count(*it) == 0 )
//                 {
//                     wCdd[*it] = numeric(mpz_get_si(mpq_numref(lps->sol[cn_so])), mpz_get_si(mpq_denref(lps->sol[cn_so])));
//                     wInverted.insert(pair<numeric,ex>(ex_to<numeric>(wCdd[*it]), (*it)));
//                 }

//                 cn_so++;
//             }
            
//             cout << " NEW CDD: " << wCdd << endl;
//             cout << " NEW INV: " << wInverted << endl;

// /*            minElemIter = min_element(wCdd.begin(),wCdd.end(),compExmapSecond);
            
//             if (wFinal.count(minElemIter->second) == 0)
//             {
//                 wFinal[minElemIter->second] = minElemIter->first;
//                 wCdd.erase(minElemIter);
//             }
// */
//         }
        

//     }
    
//     cout<< "Final w: " << localPoles << endl;
// // cout<<"Int: "<<  interior_point(set2lst(get_poles_set()),eps_w_current)<<endl;
//     cout<< "poles  " << gamma_poles <<endl;
//       BOOST_ASSERT_MSG(interior_point(gamma_poles,wch),"Not a convex polyhedron interior point");

//     cout << "\nCDD point" << endl;
//     // cout << wCdd << endl;
// /*    for(w_lst_type::const_iterator it = w_lst.begin();it!=w_lst.end();++it)
//         w_current[*it] = it->subs(wCdd);
//     eps_current = (get_symbol("eps")==wFinal[get_symbol("eps")]);
//     return eps_w_current;
// */
//     for(w_lst_type::const_iterator it = w_lst.begin();it!=w_lst.end();++it)
//         w_current[*it] = it->subs(wch);
//     eps_current = (get_symbol("eps")==wch[get_symbol("eps")]);

//     return wch;

// }








MBlst MBcontinue(MBintegral rootint,ex eps0)
{
  try
    {
      rootint.barnes1();
      rootint.barnes2();

  

      MBlst O(1,rootint);
      MBlst C;
      while(O.size()>0)
        {
          MBlst R;
          for(MBlst::iterator it = O.begin();it!=O.end();++it)
            {
              //   cout<<std::setw(15+it->get_level())<<std::right<<"shifted on "<<it->get_level()<<endl;
              C.push_back(*it);//need review, multiple entries C=C U I
              MBintegral::pole_iterator pit,pit_end;
              ex eps_i = get_symbol("eps");
              //          cout<<"after barness lemas "<<it->get_eps()<<endl;
              eps_i = eps_i.subs(it->get_eps());



              // Iterate over gamma arguments with eps dependence only!!!!!!!
              lst with_eps_lst(it->poles_with_eps());
              for(lst::const_iterator pit  = with_eps_lst.begin(); pit != with_eps_lst.end(); ++pit)
                {
                  cout<<*pit<<endl;
                  cout<<"F(eps_i) "<<pit->subs(it->get_w()).subs(it->get_eps())<<"F(eps=0) "<<pit->subs(it->get_w()).subs(get_symbol("eps")==eps0)<<"   min  "<<std::min(pit->subs(it->get_w()).subs(it->get_eps()),pit->subs(it->get_w()).subs(get_symbol("eps")==eps0))<<endl;
                 
             
                  ex F_eps0 = pit->subs( it->get_w() ).subs( get_symbol("eps") == eps0) ;
                  ex F_epsi =  pit->subs( it->get_w() ).subs( it->get_eps() ) ;

                  if(F_eps0==F_epsi) 
                    cout<<"Terminating eps=0 achieved  "<<std::min(F_eps0,F_epsi)<<endl;

                  ex dir__ = csgn(F_eps0 - F_epsi);
                  int dir = int( ex_to<numeric>(dir__).to_double());

                  //              for(int n =0;n>std::min(F_eps0,F_epsi);n--)

                  int pole;
              
                  if(dir > 0)
                    pole = int(ceil(ex_to<numeric>(F_epsi).to_double())); 
                  else if(F_epsi < 0)
                    pole = int(floor(ex_to<numeric>(F_epsi).to_double())); 
                  else
                    pole = 0;

                  for(int n = pole; dir*(F_eps0 - n) >= 0 && n <= 0; n += dir)
                    {
                      //        if( n < std::max(F_eps0,F_epsi))
                      {
                  
                        // cout<<pit->subs(it->get_w()) <<endl;
                        // test on epsilon existance
                        if(pit->subs(it->get_w()).has(get_symbol("eps")))
                          {
                            ex eps_prime = lsolve(pit->subs(it->get_w()) ==n,get_symbol("eps") );
                            //                      BOOST_ASSERT_MSG( F_eps0 < F_epsi," Wrong direction");
                            //       cout<<"solve"<<endl;
                            // cout<<"F= "<<*pit<<endl;
                            // cout<<"eps_i: "<<lsolve(pit->subs(it->get_w()) ==n,get_symbol("eps") )<<endl;
                            // cout<<" Poles of Gamma on eps_i: "<<it->get_pole_lst().subs(it->get_w()).subs(get_symbol("eps")==eps_prime)<<endl;
                            lst w_in_F  = it->has_w(*pit);
                            if(w_in_F.nops()>0)
                              {
                                cout<<endl<<"LEVEL "<<it->get_level()<<" Epsilon continue from eps_i = "<<eps_i<<" to "<<eps_prime<<endl<<endl;
                                BOOST_ASSERT_MSG(abs(ex_to<numeric>(eps_i).to_double())>abs(ex_to<numeric>(eps_prime).to_double()), "Bad continuation");
                                // cout<<lsolve(*pit==n,w_in_F.op(0))<<endl;
                                //   cout<<"sign(z) = "<<csgn(pit->coeff(w_in_F.op(0)))<<"     sign(F_i-F_0) = "<<csgn(F_epsi-F_eps0)<<endl;
                      
                                //MBintegral newi(it->get_w_lst(),it->get_pole_lst(),it->get_expr(),it->get_w(),it->get_eps());
                                //      cout<<"NEW INT: "<< newi.get_expr()<<endl;
                                // cout<<"debug fepsi"
                                //   <<"1: "<<lsolve(*pit==n,w_in_F.op(0))<<endl
                                //   <<"2: "
                                //   <<endl;
                                // cout<<"BEFORE RESIDUE!: "<<it->get_w_lst()<<endl
                                //    <<it->get_expr()<<endl;
                                MBintegral res_int = it->res(w_in_F.op(w_in_F.nops()-1)==lsolve(*pit==n,w_in_F.op(w_in_F.nops()-1)),*pit,get_symbol("eps")==eps_prime);
                                res_int.set_level(1+it->get_level());
                                //        cout<<"after RESIDUE!: "<<res_int.get_w_lst()<<endl
                                //   <<res_int.get_expr()<<endl;
                                //res_int.update_poles_from_ex();

                                //  cout<<"Storing RES_INT with eps_prime = " <<eps_prime<<"  "<<res_int.get_eps().rhs()<<endl;
                                res_int*=(2*Pi*I*csgn(pit->coeff(w_in_F.op(w_in_F.nops()-1)))*csgn(F_epsi-F_eps0));
                                // cout<<"RES EXPR:  "<<res_int.get_expr()<<endl;
                                res_int.barnes1();
                                res_int.barnes2();
                                R.push_back(res_int);
                              }
                            // else BOOST_ASSERT_MSG(false,"EEEEEERRRRRRRROOOORR: no W dependence in pole");
                          }
                        else BOOST_ASSERT_MSG(false,"EEEEEERRRRRRRROOOORR: no eps dependence in pole");
                        //cout<<endl<<endl<<"EEEEEERRRRRRRROOOORR: no W dependence in pole"<<endl<<endl;
                      }// if n>max
                    }

                }
            }
          O = R;
        }
      cout<<"Continue get "<<C.size()<<" integrals"<<endl;
  

      for(MBlst::iterator it = C.begin();it!= C.end();++it)
        {
          //cout<<(it->get_pole_lst().subs(it->get_w())).subs(get_symbol("eps")==0)<<endl;
          // cout<<endl <<"series:  " << it->get_gamma_expr().series(get_symbol("eps")==0,0)<<endl<<endl;
        }
      return C;
    }catch(std::exception &p)
    {
      throw std::logic_error(std::string("In function \"MBcontinue\":\n |___> ")+p.what());
    }
}


MBtree MBcontinue_tree(MBintegral rootint,ex eps0)
{
    using namespace mbtree;
    try 
    {
        rootint.barnes1();
        rootint.barnes2();
        
        //           accumulator for constraints
        //  cout<<"start constrainta : "<<rootint.get_poles_set()<< " "<<rootint.get_w_eps_set()<<endl;
        constr_acc ca(rootint.get_poles(),rootint.get_w_eps());
        //            tree root creation 
        MBtree C;
        MBtree::iterator lastChildIt;//,root_it;
        lastChildIt= C.insert(C.begin(), rootint);
        size_t nChildrenAdded = 0;  
        do 
        {
            nChildrenAdded = 0;           // no children added in new level
            //      for(;it!=same_depth_start_iter;next_at_same_depth (it))
            /* MBtree::fixed_depth_iterator it,it_end;
               it = C.begin_fixed (root_it, C.max_depth ());
    
               it_end = C.end_fixed (root_it, C.max_depth ());*/
            MBtree::leaf_iterator it,it_end;
            it = C.begin_leaf();
            it_end = C.end_leaf();
            cout<<"work"<<endl;
            for(it = C.begin_leaf();it != it_end; ++it ) 
            {
                if((C.depth(it) == C.max_depth() && nChildrenAdded == 0) || (C.depth(it) == C.max_depth()-1 && nChildrenAdded > 0) ) 
                {
                    cout<<"Leaf depth "<<C.depth(it)<<" Max depth "<<C.max_depth()<<endl;
                    //   cout<<std::setw(15+it->get_level())<<std::right<<"shifted on "<<it->get_level()<<endl;
                    //C.push_back(*it);//need review, multiple entries C=C U I
                    MBintegral::pole_iterator pit,pit_end;
                    ex eps_i = get_symbol("eps");
                    //          cout<<"after barness lemas "<<it->get_eps()<<endl;
                    eps_i = eps_i.subs(it->get_eps());

                    // finding nearest EPS

                    typedef multimap<ex,ex> EpsPolesMap;
                    typedef EpsPolesMap::iterator mapIter;
                    //                  while(no_continuat

                    EpsPolesMap epsPoles;

                    //             Iterate over gamma arguments with eps dependence only!!!!!!!
                    lst polesWithEps(it->poles_with_eps());
                    lst polesMinEps;
                    for(lst::const_iterator pit  = polesWithEps.begin(); pit != polesWithEps.end(); ++pit) 
                    {
                        ex fEps0 = pit->subs( it->get_w() ).subs( get_symbol("eps") == eps0) ;
                        ex fEpsI =  pit->subs( it->get_w() ).subs( it->get_eps() ) ;
                        ex poleValue;
                        if(fEpsI > fEps0) poleValue = int(floor(ex_to<numeric>(fEpsI).to_double()));
                        else if(fEpsI < fEps0) poleValue = int(ceil(ex_to<numeric>(fEpsI).to_double()));
                        if (poleValue <= 0) 
                        {
                            cout << setw(5) << right << poleValue <<  " --------------->  " 
                                 << setw(25) << left << *pit << endl;
                            ex eps_pole_sol = lsolve(pit->subs(it->get_w()) == poleValue,get_symbol("eps") );
                            polesMinEps.append(eps_pole_sol);
                            epsPoles.insert ( pair<ex,ex>(eps_pole_sol,*pit) );
                        }
                    }

                    cout << " ______________POLES NEAREST_______________  " << endl;
                    cout << "|                                          | " << endl;
                    cout << "| " << setw(40) << std::left <<it->get_w() << endl;
                    cout << "|        " << polesMinEps << endl;
                    cout << "|__________________________________________| " << endl;
                    cout << endl;//"*** min elem: "<<*min_element(polesMinEps.begin(),polesMinEps.end())<<endl;
                    
                    mapIter m_it, s_it;
                    
                    for (m_it = epsPoles.begin();  m_it != epsPoles.end();  m_it = s_it) 
                    {
                        ex theKey = (*m_it).first;
                        
                        cout << endl;
                        cout << "  key = '" << theKey << "'" << endl;

                        // working with key here:
                      
                      
                        pair<mapIter, mapIter> keyRange = epsPoles.equal_range(theKey);
                      
                        // Iterate over all map elements with key == theKey
                      
                        for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it) 
                        {
                            cout << "    value = " << (*s_it).second << endl;
                          
                            lst wInPole = it->has_w((*s_it).second);
                            wInPole.sort();
                            
                            if (epsPoles.count(theKey) > 1)                             
                            {
                                cout << endl;
                                cout << endl;
                                cout << "(((((((((((       w in pole " << wInPole << endl;
                                cout << endl;
                                cout << endl;
                            }
                        }
                    }



                    for(lst::const_iterator pit  = polesWithEps.begin(); pit != polesWithEps.end(); ++pit) 
                    {
                        cout<<"POLE EXPR:   "<< *pit <<endl;
                        
                        cout<<"F(eps_i) "<<pit->subs(it->get_w()).subs(it->get_eps())
                            <<"F(eps=0) "<<pit->subs(it->get_w()).subs(get_symbol("eps")==eps0)
                            <<"   min  "<<std::min(pit->subs(it->get_w()).subs(it->get_eps()),pit->subs(it->get_w()).subs(get_symbol("eps")==eps0))<<endl;
                    
                        cout<<"eps_i = "<<eps_i<<"  w current: "<<it->get_w()<<endl;
                    
                        ex fEps0 = pit->subs( it->get_w() ).subs( get_symbol("eps") == eps0) ;
                        ex fEpsI =  pit->subs( it->get_w() ).subs( it->get_eps() ) ;
                    
                        if(fEps0==fEpsI) 
                        {
                            cout<<"Terminating, eps=0 achieved  "<<std::min(fEps0,fEpsI)<<endl;
                            throw std::logic_error(string("Contour hit the pole"));
                            assert(false);
                        }

                      
                        else 
                        {
                      
                            ex dir__ = csgn(fEps0 - fEpsI);
                            int dir = int( ex_to<numeric>(dir__).to_double());
                      
                            //              for(int n =0;n>std::min(fEps0,fEpsI);n--)
                      
                            int pole;
                      
                            if(dir > 0)
                                pole = int(ceil(ex_to<numeric>(fEpsI).to_double())); 
                            else if(fEpsI < 0)
                                pole = int(floor(ex_to<numeric>(fEpsI).to_double())); 
                            else
                                pole = 0;

                            for(int n = pole; dir*(fEps0 - n) >= 0 && n <= 0; n += dir) 
                            {
                                //        if( n < std::max(fEps0,fEpsI))
                            
			    
                                //             test on epsilon existance
                                if(pit->subs(it->get_w()).has(get_symbol("eps"))) 
                                {
                                    ex eps_prime = lsolve(pit->subs(it->get_w()) ==n,get_symbol("eps") );
                                    lst w_in_F  = it->has_w(*pit);
                                    if(w_in_F.nops()>0) 
                                    {
                                        cout<<endl<<endl<<endl<<"ASDASDASDASDSD "<<*pit<<" "<<ca.test_single(pit->subs(get_symbol("eps") == 0))<<endl<<endl<<endl;
                                        
                                        cout<<endl<<"LEVEL "<<it->get_level()<<" Epsilon continue from eps_i = "<<eps_i<<" to "<<eps_prime<<endl<<endl;
                                        BOOST_ASSERT_MSG(abs(ex_to<numeric>(eps_i).to_double())>=abs(ex_to<numeric>(eps_prime).to_double()), "Bad continuation");
                                        
                                        //             decide what var to get res
                                        ex var_to_get_res = 0;
                                        for(lst::const_iterator vgit = w_in_F.begin(); vgit != w_in_F.end();++vgit) 
                                        {
                                            if(pit->coeff(*vgit,1) == 1) 
                                            {
                                                var_to_get_res = *vgit;
                                                break;
                                            }
                                            if(pit->coeff(*vgit,1) == -1)
                                                var_to_get_res = *vgit;
                                        }
                                        if( var_to_get_res ==0) var_to_get_res = w_in_F.op(w_in_F.nops()-1);
                                    
                                        cout<<" POLE: " << *pit<< "       var to get res   "<<var_to_get_res<<endl;

                                        MBintegral res_int = it->res(var_to_get_res ==lsolve(*pit==n,var_to_get_res),*pit,get_symbol("eps")==eps_prime);
                                        res_int.set_level(1+it->get_level());
                                        res_int*=(2*Pi*I*csgn(pit->coeff(var_to_get_res))*csgn(fEpsI-fEps0));
                                        //if(ca.test_single(pit->subs(get_symbol("eps") == 0)))res_int.set_optimizable(true);
                                        res_int.barnes1();
                                        res_int.barnes2();
                                        //R.push_back(res_int);
                                        //  if(eps_prime !=0)
                                        lastChildIt = C.append_child(it,res_int);
                                        nChildrenAdded++;
                                    }
                                    // else BOOST_ASSERT_MSG(false,"EEEEEERRRRRRRROOOORR: no W dependence in pole");
                                }
                                else BOOST_ASSERT_MSG(false,"EEEEEERRRRRRRROOOORR: no eps dependence in pole");
                                //cout<<endl<<endl<<"EEEEEERRRRRRRROOOORR: no W dependence in pole"<<endl<<endl;
                                // if n>max
                            }
                        }//else
                    }
                }//same depth
            }
            //      O = R;
        } while(nChildrenAdded > 0);
        cout<<"Continue get "<<C.size()<<" integrals"<<endl;
        // kptree::print_tree_bracketed(C);

        return C;
    }catch(std::exception &p)
    {
        throw std::logic_error(std::string("In function \"MBcontinue\":\n |___> ")+p.what());
    }
}



ex expand_and_integrate(MBintegral& int_in, lst num_subs, int expansion_order) // up to O(eps^1) 
{
  try
    {
      ex out_ex;

      //            cout<<(it->get_pole_lst().subs(it->get_w())).subs(get_symbol("eps")==0)<<endl;
      //         cout<<endl <<"series:  " << it->get_gamma_expr().series(get_symbol("eps")==0,expansion_order)<<endl<<endl;
      lst w_lst = int_in.get_w_lst();
      if(int_in.w_size()>0)// expanding and integrating
        {
          int_in.barnes1();
          int_in.barnes2();
          out_ex = series_to_poly( int_in.get_expr().series(get_symbol("eps"),expansion_order) ).subs(num_subs);
          //out_ex = series_to_poly( int_in.get_expr().series(int_in.get_eps(),expansion_order) ).subs(num_subs);
          // loop over W_i, converting integration contour
          for(lst::const_iterator wit = w_lst.begin(); wit != w_lst.end(); ++wit)
            {
              ex c_i = wit->subs(int_in.get_w());
              out_ex = (I*out_ex.subs((*wit)==c_i - I*log( (*wit)/( 1 - (*wit) ) ) ) ) / (*wit)/(1- (*wit));
            }
          cout<<"current : "<<out_ex<<endl;
          ex wo_eps_part = out_ex;
          ex vegas_ex = 0;
          for(int i = out_ex.ldegree( get_symbol("eps") ); i <expansion_order/* out_ex.degree( get_symbol("eps") )*/; i++)
            {
              cout<<"Ord( "<<i<<" ) coeff : "<< out_ex.coeff(get_symbol("eps"),i)<<endl;
              ex int_expr =  out_ex.coeff(get_symbol("eps"),i);
              RoMB::FUNCP_CUBA2 fp_real;
              RoMB::compile_ex_real(lst(evalf(int_expr)),int_in.get_w_lst(), fp_real);

              // ----------------------------------- Vegas integration-------------------------
              int  NDIM  = int_in.w_size();
              //#define NCOMP 1
#define USERDATA NULL
#define EPSREL 1e-3
#define EPSABS 1e-9
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 100000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
              //SUAVE
#define NNEW 1000
#define FLATNESS 25.

              //#define KEY1 47
              //#define KEY2 1
              //#define KEY3 1
              //#define MAXPASS 5
              //#define BORDER 0.
              //#define MAXCHISQ 10.
              //#define MINDEVIATION .25
              //#define NGIVEN 0
              //#define LDXGIVEN NDIM
              //#define NEXTRA 0

#define KEY 0
              const int NCOMP = 1;
              int comp,  neval, fail;
              double integral_real[NCOMP], error[NCOMP], prob[NCOMP];
                   
              printf("-------------------- Vegas test --------------------\n");
	      /*              
			      Vegas(NDIM, NCOMP, fp,
			      EPSREL, EPSABS, VERBOSE, 
			      MINEVAL, MAXEVAL,
			      NSTART, NINCREASE,
			      &neval, &fail, integral, error, prob);
              */     
	     
	      Vegas(NDIM, NCOMP, fp_real, USERDATA,
		    EPSREL, EPSABS, VERBOSE, SEED,
                    MINEVAL, MAXEVAL, 
		    NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE,
		    &neval, &fail, integral_real, error, prob);


          
              // printf("VEGAS RESULT:\tneval %d\tfail %d\n",
              
              for( comp = 0; comp < NCOMP; ++comp )
                printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                       integral_real[comp], error[comp], prob[comp]);
                       
              // ----------------------------------- Vegas integration-------------------------              
              vegas_ex+=pow(get_symbol("eps"),i)*(integral_real[0]);//+I*integral_imag[0]);
            }
          out_ex = vegas_ex;
        }
      else // expanding only
        {
          out_ex = series_to_poly( int_in.get_expr().series(get_symbol("eps"),expansion_order) ).subs(num_subs);
        }
      cout<<endl<<"gp "<<out_ex<<endl;
      return out_ex;
    }catch(std::exception &p)
    {
      throw std::logic_error(std::string("In function \"Expand_and_integrate\":\n |___> ")+p.what());
    }
}


ex expand_and_integrate_complex(MBintegral& int_in, lst num_subs, int expansion_order) // up to O(eps^1) 
{
  try
    {
      ex out_ex;

      //            cout<<(it->get_pole_lst().subs(it->get_w())).subs(get_symbol("eps")==0)<<endl;
      //         cout<<endl <<"series:  " << it->get_gamma_expr().series(get_symbol("eps")==0,expansion_order)<<endl<<endl;
      lst w_lst = int_in.get_w_lst();
      if(int_in.w_size()>0)// expanding and integrating
        {
          int_in.barnes1();
          int_in.barnes2();
          out_ex = series_to_poly( int_in.get_expr().series(get_symbol("eps"),expansion_order) ).subs(num_subs);
          //out_ex = series_to_poly( int_in.get_expr().series(int_in.get_eps(),expansion_order) ).subs(num_subs);
          // loop over W_i, converting integration contour
          for(lst::const_iterator wit = w_lst.begin(); wit != w_lst.end(); ++wit)
            {
              ex c_i = wit->subs(int_in.get_w());
              out_ex = (I*out_ex.subs((*wit)==c_i - I*log( (*wit)/( 1 - (*wit) ) ) ) ) / (*wit)/(1- (*wit));
            }
          cout<<"current : "<<out_ex<<endl;
          ex wo_eps_part = out_ex;
          ex vegas_ex = 0;
          for(int i = out_ex.ldegree( get_symbol("eps") ); i <expansion_order/* out_ex.degree( get_symbol("eps") )*/; i++)
            {
              cout<<"Ord( "<<i<<" ) coeff : "<< out_ex.coeff(get_symbol("eps"),i)<<endl;
              ex int_expr =  out_ex.coeff(get_symbol("eps"),i);
              RoMB::FUNCP_CUBA2 fp_real,fp_imag;
              RoMB::compile_ex_real(lst(int_expr),int_in.get_w_lst(), fp_real);
              RoMB::compile_ex_imag(lst(int_expr),int_in.get_w_lst(), fp_imag);
              
              // ----------------------------------- Vegas integration-------------------------
              int  NDIM  = int_in.w_size();
              //#define NCOMP 1
#define USERDATA NULL
#define EPSREL 1e-3
#define EPSABS 1e-9
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 100000
              
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
              //SUAVE
#define NNEW 1000
#define FLATNESS 25.
              
              //#define KEY1 47
              //#define KEY2 1
              //#define KEY3 1
              //#define MAXPASS 5
              //#define BORDER 0.
              //#define MAXCHISQ 10.
              //#define MINDEVIATION .25
              //#define NGIVEN 0
              //#define LDXGIVEN NDIM
              //#define NEXTRA 0

#define KEY 0
              const int NCOMP = 1;
              int comp, nregions, neval, fail;
              double integral[NCOMP],integral_real[NCOMP],integral_imag[NCOMP], error[NCOMP], prob[NCOMP];
                   
              printf("-------------------- Vegas test --------------------\n");
	      /*              
			      Vegas(NDIM, NCOMP, fp,
			      EPSREL, EPSABS, VERBOSE, 
			      MINEVAL, MAXEVAL,
			      NSTART, NINCREASE,
			      &neval, &fail, integral, error, prob);
              */     
	     
	      Vegas(NDIM, NCOMP, fp_real, USERDATA,
		    EPSREL, EPSABS, VERBOSE, SEED,
                    MINEVAL, MAXEVAL, 
		    NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE,
		    &neval, &fail, integral_real, error, prob);


          
              // printf("VEGAS RESULT:\tneval %d\tfail %d\n",
              
              for( comp = 0; comp < NCOMP; ++comp )
                printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                       integral_real[comp], error[comp], prob[comp]);
                       
              /*   TEMPORARY no IMAGE PART
                   Vegas(NDIM, NCOMP, fp_imag, USERDATA,
                   EPSREL, EPSABS, VERBOSE, SEED,
                   MINEVAL, MAXEVAL, 
                   NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE,
                   &neval, &fail, integral_imag, error, prob);
                   for( comp = 0; comp < NCOMP; ++comp )
                   printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                   integral_imag[comp], error[comp], prob[comp]);
              */

            
              /*
                printf("\n-------------------- Suave test real--------------------\n");
               
                Suave(NDIM, NCOMP, fp_real, USERDATA,
                EPSREL, EPSABS, VERBOSE | LAST, SEED,
                MINEVAL, MAXEVAL, NNEW, FLATNESS,
                &nregions, &neval, &fail, integral_real, error, prob);
                             
                printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
                nregions, neval, fail);
                                   
                for( comp = 0; comp < NCOMP; ++comp )
                printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                integral_real[comp], error[comp], prob[comp]);


                printf("\n-------------------- Suave test imag--------------------\n");
               
                Suave(NDIM, NCOMP, fp_imag, USERDATA,
                EPSREL, EPSABS, VERBOSE | LAST, SEED,
                MINEVAL, MAXEVAL, NNEW, FLATNESS,
                &nregions, &neval, &fail, integral_imag, error, prob);
                             
                printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
                nregions, neval, fail);*/
                                   
             
              // ----------------------------------- Vegas integration-------------------------              
              vegas_ex+=pow(get_symbol("eps"),i)*(integral_real[0]);//+I*integral_imag[0]);
            }
          out_ex = vegas_ex;
        }
      else // expanding only
        {
          out_ex = series_to_poly( int_in.get_expr().series(get_symbol("eps"),expansion_order) ).subs(num_subs);
        }
      cout<<endl<<"gp "<<out_ex<<endl;
      return out_ex;
    }catch(std::exception &p)
    {
      throw std::logic_error(std::string("In function \"Expand_int_complex\":\n |___> ")+p.what());
    }
}

