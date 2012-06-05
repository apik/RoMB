//#define GMPRATIONAL
#include <boost/assert.hpp>
#include "mbintegral.h"
#include "utils.h"

#include "constracc.h"

//#include <setoper.h>

//#include <ppl.hh>
//using namespace Parma_Polyhedra_Library;
//using namespace Parma_Polyhedra_Library::IO_Operators;
//#include "ppl_interface.h"



//#include "tree_util.hh"
/**
 *
 *  Construction from F, U=1
 *
 */

MBintegral::MBintegral(
                       UFXmap fx_in,
                       lst nu,
                       GiNaC::numeric l,
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
      //   exmap new_w_current(w_current);
      //   cout<<" Not modif:  "<<new_w_current<<endl;

      //new_w_current.erase(w_relation.lhs());
      MBintegral resINT(lst(cut_w_vec.begin(),cut_w_vec.end()),res_loran,new_eps,tree_level+1);
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

bool compExmapSecond(exmap::value_type a, exmap::value_type b)
{
    return (a.second < b.second);
}




class CompExmapFirst : public std::binary_function<GiNaC::numeric,GiNaC::numeric,bool>
{
public:
    bool operator()(const GiNaC::numeric& lh, const GiNaC::numeric& rh)
        {
//            if(is_a<GiNaC::numeric>(lh) && is_a<GiNaC::numeric>(rh))
                return (lh < rh);
                //           else
                //  return lh.compare(rh);
        }
};

template <typename T>
bool alwaysTrue (T i) { return true; }




