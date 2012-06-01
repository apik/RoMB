#include "romb.h"
#include "utils.h"
#include "shift.h"
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

extern "C"
{
#include <gsl/gsl_integration.h>
}
/**
 *
 *  loop momentums,propagator expressions,
 *  invariants substitutions,propagator powers,number of loops
 *
 \param k_lst loop momentums list
 \param p_lst propagator expressions list
 \param subs_lst invariants substitutions list
 \param nu propagator powers list
 \param l number of loops 
   
 \return 
 *
 */


RoMB_loop_by_loop:: RoMB_loop_by_loop(
				      lst k_lst,
				      lst p_lst,
				      lst subs_lst,
                                      lst nu,
                                      bool subs_U
                                      )
{
  try
    {



      /* 
	 empty integral
      */
      MBintegral MBlbl_int(lst(),lst(),1);
      /* 
	 Full set of unused propagators, will change 
      */
      exlist input_prop_set;//( p_lst.begin(),p_lst.end());
      /* 
	 map for propagator powers
      */
      exmap prop_pow_map;
      for(lst::const_iterator Pit = p_lst.begin(); Pit != p_lst.end(); ++Pit)
	{
          input_prop_set.push_back(Pit->expand());
          prop_pow_map[Pit->expand()] = nu.op(std::distance(p_lst.begin(),Pit));
	}
	
      cout<<"INPSET: "<<input_prop_set<<endl;

      /* 
	 Iterate over momentums k1,k2,k3,etc.
      */
      unsigned int displacement_x = 0;
      unsigned int displacement_w = 0;
      for(lst::const_iterator kit = k_lst.begin(); kit != k_lst.end(); ++kit)
	{
          // Integral Normalization coefficient 
        //            MBlbl_int *= pow(I,k_lst.nops());
           MBlbl_int *= 1/tgamma(1+get_symbol("eps"));
          //MBlbl_int *= pow(Pi,2-get_symbol("eps"));
          //MBlbl_int *= exp(Euler*get_symbol("eps"));	  
          cout<<"PROP_POW_MAP "<<prop_pow_map<<endl;
	  /*
	    temporary set of propagators, with all momentum,except deleted
	  */
	  exlist tmp_p_lst(input_prop_set.begin(), input_prop_set.end()); 
	  /*
	    temporary set of propagators, with KIT momentum
	  */
	  lst P_with_k_lst;
	  BOOST_FOREACH(ex prop_tmp, tmp_p_lst)
	    {
	      if(prop_tmp.has(*kit))
		{
		  P_with_k_lst.append(prop_tmp);
		  input_prop_set.remove(prop_tmp);
		}
	    }
            
	  cout<< "Set wo k_i "<<input_prop_set<<endl;
	  cout<<" PWKlst "<<P_with_k_lst<<endl;
	  bool direct_formula_applied = false;
          // if only one term in PWKLST use well known formulas
          // [Smirnov A.1]
          if(!direct_formula_applied && (P_with_k_lst.nops() == 1))
            {
              ex pr_t = P_with_k_lst.op(0);
              ex nu_t = prop_pow_map[pr_t];
              exmap repls;
              BOOST_ASSERT_MSG(pr_t.match(-pow(*kit,2) + wild(2)),"ONE PROP");
              if(pr_t.match(-pow(*kit,2) + wild(2),repls))cout<<"repls: "<<repls<<endl;
              ex mass_tadpole = (tgamma(nu_t+get_symbol("eps")-2)/tgamma(nu_t)*pow(wild(2).subs(repls),-nu_t-get_symbol("eps")+2));
              cout<<mass_tadpole<<endl;
              MBlbl_int *= mass_tadpole;
              MBlbl_int.add_pole(nu_t+get_symbol("eps")-2);
              direct_formula_applied = true;
            }
          if(!direct_formula_applied && (P_with_k_lst.nops() == 2)) 
            {
              //TWO terms in PWK_LST, [Smirnov A.4]
              exmap repls_tad;
              if((P_with_k_lst.nops()==2) &&
        	 (( (P_with_k_lst.op(0).match(-pow(*kit,2))) && (P_with_k_lst.op(1).match(-pow(*kit,2)+wild())))||
                  ( (P_with_k_lst.op(1).match(-pow(*kit,2))) && (P_with_k_lst.op(0).match(-pow(*kit,2)+wild()))))
                 && !wild().has(*kit)
                 ) 
                {
                  cout<<"Two prop tadpole "<<wild()<<endl;
                  exmap r1,r2;
                  ex mm,lmb1,lmb2;
                  if( (P_with_k_lst.op(0).match(-pow(*kit,2))) && (P_with_k_lst.op(1).match(-pow(*kit,2)+wild(),r1)))
                    {
                      lmb1 = prop_pow_map[P_with_k_lst.op(1)];
                      lmb2 = prop_pow_map[P_with_k_lst.op(0)];
                      mm=wild().subs(r1);
                    }
                  else if( (P_with_k_lst.op(1).match(-pow(*kit,2))) && (P_with_k_lst.op(0).match(-pow(*kit,2)+wild(),r2)))
                    {
                      lmb1 = prop_pow_map[P_with_k_lst.op(0)];
                      lmb2 = prop_pow_map[P_with_k_lst.op(1)];
                      mm=wild().subs(r2);
                    }
                  else throw std::logic_error(std::string("Wrong two prop topology to use eq [Smir:A.4]"));
                  ex mass_tadpole = tgamma(lmb1+lmb2+get_symbol("eps")-2)*tgamma(-lmb2-get_symbol("eps")+2)/tgamma(lmb1)/tgamma(2-get_symbol("eps"))*pow(mm,-lmb1-lmb2-get_symbol("eps")+2);
                  cout<<mass_tadpole<<endl;
                }
            }
          if(!direct_formula_applied)
            {
              //          cout<< " coe: "<<coe_prop_lst<<endl;
              /*
                lexi sort of input prop list, and it's modification
              */            
              
              
              // uf and then MB represenatation construction
              // subs only in F for last momentum
              UFXmap inUFmap;
              if(boost::next(kit) == k_lst.end())
                    inUFmap = UF(lst(*kit),P_with_k_lst,subs_lst,displacement_x);
                  else
                    inUFmap = UF(lst(*kit),P_with_k_lst,subs_lst,displacement_x); // no substitution!!!
                  displacement_x +=fusion::at_key<UFX::xlst>(inUFmap).nops(); 

                  lst nu_into;
                  for(lst::const_iterator nuit = P_with_k_lst.begin(); nuit != P_with_k_lst.end(); ++nuit )
                    nu_into.append(prop_pow_map[*nuit]);
                  cout<<" Powers list before input: "<<nu_into<<endl;
                  /*
                    MBintegral Uint(
                    fusion::make_map<UFX::F,UFX::xlst>(fusion::at_key<UFX::F>(inUFmap),
                    fusion::at_key<UFX::xlst>(inUFmap)
                    ),nu_into,1,displacement_w);
                  */

                  MBintegral Uint(inUFmap,nu_into,1,subs_U,displacement_w);
                  displacement_w+=Uint.w_size();
                  cout<<"ui9nt eps (no gamma) : "<<Uint.get_expr().subs(tgamma(wild()) == 1)<<endl;
                  cout<<"ui9nt eps : "<<Uint.get_expr()<<endl;
                  /*
                    expression to mul root integral
                    where to subs prop(k_prev)==1
                  */
                  ex expr_k_to_subs_1= Uint.get_expr();

                  ex mom_find = Uint.get_expr();
                  cout<< "where find props: "<<mom_find<<endl;
                  if(is_a<mul>(mom_find))
                    {
                      // set of a^b?, need to have momentums from *kit to *k_lst.end()
                      exset found_prop_raw,found_prop;
                      mom_find.find(pow(wild(1),wild(2)),found_prop_raw);
                      cout<<" is a mul raw "<<found_prop_raw<<endl;
                      // really props
                      BOOST_FOREACH(ex px_c,found_prop_raw)
                        {
                          bool is_a_p = false;
                          for(lst::const_iterator kpi = kit; kpi != k_lst.end(); ++kpi)
                            if(px_c.has(*kpi))
                              {
                                is_a_p = true;
                                break;
                              }
                      
                          if(is_a_p) found_prop.insert(px_c);
                        }
                      cout<<" is a mul  "<<found_prop<<endl;
                      BOOST_FOREACH(ex propex_c,found_prop)
                        {
                          // next momentum in loop momentum list
                          ex next_k ;
                          for(lst::const_iterator nkit = kit; nkit != k_lst.end(); ++nkit)
                            if(propex_c.has(*nkit))
                              {
                                next_k = *nkit;
                                break;
                              }                         
                      
                      
                          cout<<"before subs kex : "<<expr_k_to_subs_1.subs(tgamma(wild()) == 1)<<endl;
                          expr_k_to_subs_1 = expr_k_to_subs_1.subs(propex_c == 1);
                          cout<<"after subs kex : "<<expr_k_to_subs_1.subs(tgamma(wild()) == 1)<<endl;
                          /*
                            converting prop to form -p^2+m^2
                          */
                      
                          ex p_power;
                          ex p_expr;
                          ex p_not_corr = ex_to<power>(propex_c).op(0);
                          ex coeff_ksq = p_not_corr.expand().coeff(next_k,2); // coeff infront of K^2
                          if( coeff_ksq != -1 )
                            {
                              p_not_corr /=coeff_ksq;
                              cout<<"koeff_ksq "<<coeff_ksq<<endl;
                              MBlbl_int*= pow(coeff_ksq,ex_to<power>(propex_c).op(1));
                              //                  propex = pow(p_not_corr,ex_to<power>(propex_c).op(1));
                              p_power = ex_to<power>(propex_c).op(1);
                              p_expr = p_not_corr.expand();
                            }
                          else
                            {
                              p_power = ex_to<power>(propex_c).op(1);
                              p_expr = ex_to<power>(propex_c).op(0).expand();
                            }
                          
                          /*
                            Search for duplications in prop set
                          */
                          cout<<"where to find props: "<<input_prop_set<<endl;
                          cout<< input_prop_set.size()<<endl;
                          cout<<"PWK_MAP to modiff"<<prop_pow_map<<endl;
                          if(prop_pow_map.count(p_expr) > 0)
                            {
                              BOOST_ASSERT_MSG(count(input_prop_set.begin(),input_prop_set.end(),p_expr) > 0,"Propagator not found in prop set");
                              cout<<"PPM bef: "<< prop_pow_map[p_expr]<<endl;
                              prop_pow_map[p_expr] -=p_power;
                              cout<<"PPM aft: "<< prop_pow_map[p_expr]<<endl;
                            }
                          else
                            {
                              prop_pow_map[p_expr] = (-1)*p_power;
                              input_prop_set.push_back(p_expr);
                            }
                          cout<<"PWK_MAP after modiff"<<prop_pow_map<<endl;
                        }
                    }
                  //cout<<"needed props "<<prop_pow_lst<<endl;
              
              
              
                  cout<<"ya tut"<<endl;            
              
              
                  MBlbl_int*=expr_k_to_subs_1;
                  cout<<"ya tut"<<endl;            
                  //          cout<<"HAS INT :    "<<MBlbl_int.get_expr().subs(tgamma(wild(4)) == 0)<<endl;
                  MBlbl_int+=Uint;
                  cout<<"bad"<<endl;
                  //            MBlbl_int.insert_w_lst(Uint.get_w_lst());
                  // MBlbl_int.insert_pole_lst(Uint.get_pole_lst());
                  //}// else more then one prop            
                }// two prop formula
          
            }

          MBlbl_int.fix_inv();
          cout<<"expr: "<<MBlbl_int.get_expr()<<endl;
          
          cout<<"Constructed integral with:"<<endl;
          //cout<<"Poles: "<<MBlbl_int.get_poles_set()<<endl;
          //MBlbl_int.set_poles_set(MBlbl_int.poles_from_ex(MBlbl_int.get_expr()));
          cout<<"Poles true: "<<MBlbl_int.poles_from_ex(MBlbl_int.get_expr())<<endl;
          cout<<"Poles: "<<MBlbl_int.get_poles()<<endl;
          cout<<"W's : "<<MBlbl_int.get_w_lst()<<endl;
          cout<<"Expr : "<<MBlbl_int.get_expr()<<endl;



// Constraints needed for start point 
          
          constraints_ =  ConstrAcc(MBlbl_int);

//          cout << 
          constraints_.PrintPoint();// <<endl;

          MBlbl_int.SetEps(get_symbol("eps") == constraints_.GetPoint().find(get_symbol("eps"))->second);


          if(MBlbl_int.w_size()>0)          
            {
                
//              MBlbl_int.newPoint();
//              w_shared = MBlbl_int.get_w();
                
                w_shared = constraints_.GetWs();
                    
              MBlbl_int.barnes1();
            }
          print_mathematica(MBlbl_int);
          // MB only if w's existsts or eps!=0
//          if((MBlbl_int.w_size() > 0) && (MBlbl_int.get_eps().rhs() != 0))

          if((MBlbl_int.w_size() > 0) && (constraints_.GetPoint().find(get_symbol("eps"))->second != 0 ))
          {
            
            
          // setting shared point
         


          cout<< endl<<" Ready for MBcontinue?  [Y/n]: ";
          char in_ch;
          std::cin>>in_ch;
          if(in_ch=='n')  exit(0);//assert(false);



          //                assert(false);

          //        MBlst int_lst = MBcontinue(MBlbl_int);
          

          int_lst = MBcontinue(MBlbl_int);

          //exit(5);//assert(false);
          
          cout<<int_lst.size()<<endl;
          /*BOOST_FOREACH(MBintegral mbint, int_lst)
            {
            cout<< "Integral:"<<endl
            << mbint.get_expr()<<endl
            <<"Point:"<<mbint.get_w()<<endl
            <<"Wlst: "<<mbint.get_w_lst()<<endl
            <<"EPS:"<<mbint.get_eps().rhs()<<endl
            <<"LEVEL: "<<mbint.get_level()<<endl<<endl;
            }*/
          cout<<"Integrals MBC "<<int_lst.size()<<endl;
//          cout<<"Integrals MBT "<<inttr.size()<<endl;
          merge();
          cout<< endl<<" Next step?  [Y/n]: ";
  

          std::cin>>in_ch;
          if(in_ch=='n')  exit(0);//assert(false);
          /*
            ex int_expr_out = 0;
            for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
            int_expr_out+=expand_and_integrate(*it,1);
        
            cout<<" FRESULT anl : "<<"          = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
            cout<<" FRESULT num: "<<"          = "<<evalf(int_expr_out.expand().collect(get_symbol( "eps" )))<<endl;
          */
          //    lst U_lst = ex_to<lst>(fusion::at_key<UFX::U>(inUFmap));
          //  cout<<"U: "<< fusion::at_key<UFX::U>(inUFmap)<<endl;//<<" U_lst: "<<U_lst<<endl;


          //        ex U_sum =fusion::at_key<UFX::U>(inUFmap).collect(  fusion::at_key<UFX::xlst>(inUFmap),true); 
          // cout<<"U_sum: "<<U_sum<<endl;
          // MBintegral Uint(
          //                fusion::make_map<UFX::F,UFX::xlst>(fusion::at_key<UFX::F>(inUFmap),
          //                                                   fusion::at_key<UFX::xlst>(inUFmap)
          //                                                  ),nu,l);
          //     cout<<"ui9nt eps : "<<Uint.get_eps().lhs()<<endl;
          // Uint.new_point();
        
          // MBlst int_lst = MBconinue(Uint);
          //ex int_expr_out = 0;
          //for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
          // int_expr_out+=expand_and_integrate(*it,1);
        
          //cout<<" RESULT : "<<endl
          //    <<"               = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
          }
          else // no contour integration or no continuation in eps
          {
              cout << "No continuation, yet" << endl;
            int_lst.push_back(MBlbl_int);
            merge();
          }
        }catch(std::exception &p)
           {
             throw std::logic_error(std::string("In function \"RoMB\":\n |___> ")+p.what());
           }
    }


  struct mathematica_replace : public map_function
  {
    ex operator()(const ex &e)
    {
      exmap repls;
      if(e.match(log(wild()),repls))
      {
          stringstream str;
          str<< "Log[" << wild().subs(repls) <<"]";
          return get_symbol(str.str());
      }     
      else if(e.match(tgamma(wild()),repls))
        {
          stringstream str;
//          str<<"Gamma["<<wild().subs(repls)<<"]";
          str << "Gamma[";
          str << (wild().subs(repls)).map(*this);
          str<<"]";

          return get_symbol(str.str());
        }
      else if(e.match(exp(wild()),repls))
        {
          stringstream str;
          str<<"Exp["<<wild().subs(repls)<<"]";
          return get_symbol(str.str());
        }
      else if(e.match(psi(wild()),repls))
      {
          stringstream str;
          str<<"PolyGamma[0,"<<wild().subs(repls).map(*this)<<"]";
          return get_symbol(str.str());
      }
      else if(e.match(psi(wild(),wild(1)),repls))
      {
          stringstream str;
          str<<"PolyGamma[";
          str<<wild().subs(repls) << "," << wild(1).subs(repls).map(*this) <<"]";
          return get_symbol(str.str());
      }
      else if (e.info(info_flags::polynomial))
      {
          stringstream str;
          str << e;
          return get_symbol(str.str());
      }

      else return e.map(*this);
      
    }
  };
/** 
 * 
 * 
 * @param mb_in 
 */
void print_mathematica_ex(ex mfin)
  {
      //mfin = mb_in.get_expr();
//    mfin = mfin.subs(w_to_z);
    mfin = mfin.subs(Euler == get_symbol("EulerGamma"));
    mathematica_replace math_map;
    ex rpm = math_map(mfin);

 

    cout<<endl<<"(********* Begin of Mathematica output **********)"<<endl;
    cout<<"fin="<<rpm<<endl;
    cout<<"(*********  End of Mathematica output   **********)"<<endl<<endl;

  }


void RoMB_loop_by_loop::print_mathematica(MBintegral mb_in)
  {
//    exmap w_c(mb_in.get_w());

      exmap w_c(this->constraints_.GetPoint());


    exmap w_to_z;
    MBintegral::w_lst_type w_l(mb_in.get_w_lst());
    size_t cntr = 1;
    BOOST_FOREACH(ex w_inl,w_l)
      {
        string st1("z");
      
        w_to_z[w_inl] = get_symbol(st1+lexical_cast<string>(cntr));
        cntr++;
      }
    //  cout<<"W TO Z: "<<w_to_z<<endl;
    ex mfin = mb_in.get_expr();
    mfin = mfin.subs(w_to_z);
    mfin = mfin.subs(Euler == get_symbol("EulerGamma"));
    mathematica_replace math_map;
    ex rpm = math_map(mfin);

    ex eps_value = mb_in.get_eps();

    cout<<endl<<"(********* Begin of Mathematica output **********)"<<endl;
    cout<<"<< MB.m"<<endl<< "<<AMBREv1.2.m"<<endl;
    cout<<"rules={{eps->"<<eps_value.rhs()<<"},{";
    for(exmap::iterator it =w_c.begin(); it != w_c.end(); ++it )
      if(boost::next(it) == w_c.end())
        cout<<w_to_z[it->first]<<"->"<<it->second<<"}}"<<endl;
      else
        cout<<w_to_z[it->first]<<"->"<<it->second<<",";
    cout<<"fin="<<rpm<<endl;
    cout<<"tcont = MBcontinue[fin, eps -> 0, rules, Residues -> False];"<<endl;
    cout<<"MBexpand[{tcont}, Exp[2 EulerGamma eps], {eps, 0, -1}]"<<endl;
    cout<<"MBintegrate[%, {ms -> 1}]"<<endl;
    cout<<"(*********  End of Mathematica output   **********)"<<endl<<endl;

  }
/*

  void RoMB_loop_by_loop::integrate(lst number_subs_list, int exp_order)
  {
    ex int_expr_out = 0;
    for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
      int_expr_out+=expand_and_integrate(*it, number_subs_list, exp_order);
    cout<<" FRESULT for parameters: "<<number_subs_list<<endl<<endl;
    cout<<" FRESULT anl : "<<"          = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
    cout<<" FRESULT num: "<<"          = "<<evalf(int_expr_out.expand().collect(get_symbol( "eps" )))<<endl;
  }
*/

  void RoMB_loop_by_loop::merge()
  {
    for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
      {
      it->barnes1();
        MBintegral::w_lst_type wl(it->get_w_lst());
        wl.sort();
        intmap[wl] += it->get_expr();
      }
    cout<<">>>>> Merged to : "<<intmap.size()<<" integrals" <<endl;
  }


  void RoMB_loop_by_loop::integrate_map(lst number_subs_list, int exp_order)
  {
    typedef  std::pair<MBintegral::w_lst_type,ex> intpair;
    ex int_expr_out = 0;
    std::map<int,ex> squared_error;

    BOOST_FOREACH(intpair i_in_m, intmap)
      {
        ex vegex,veger;
        tie(vegex,veger) = expand_and_integrate_map(i_in_m.second,i_in_m.first,w_shared, number_subs_list, exp_order);
        int_expr_out += vegex;
        for(int i = veger.ldegree( get_symbol("eps") ); i <=veger.degree( get_symbol("eps") ); i++)
          {
            squared_error[i] += pow(veger.coeff(get_symbol("eps"),i), 2);
          }
      }
    cout<<" FRESULT for parameters: "<<number_subs_list<<endl<<endl;
    cout<<" FRESULT anl : "<<"          = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
    cout<<" FRESULT num: "<<"          = "<<evalf(int_expr_out.expand().collect(get_symbol( "eps" )))<<endl;

    for(int i = int_expr_out.ldegree( get_symbol("eps") ); i <=int_expr_out.degree( get_symbol("eps") ); i++)
      {
        cout<<" eps^"<<i<<" term: "<< int_expr_out.expand().collect(get_symbol( "eps" )).coeff(get_symbol("eps"),i)
            <<" +/- "<<sqrt(squared_error[i])<<endl;
      }


  }

template<class CubaFunc> 
struct  GSL1dintAdapter 
{
  
  static double f (double x, void * params) 
  { 
    
    CubaFunc function = reinterpret_cast< CubaFunc > (params);
    // get pointer to data from gsl_vector
    
    //    typedef int (*FUNCP_CUBA2) (const int*, const double[], const int*, double[],void*);
//    double x_in[1];
//    x_in[0] = x;

    double x_in[] = {x};
//    x_in[0] = x;

    int n_d = 1;
    int n_d_f = 1;
    double f_out[1];
    (*function)( &n_d, x_in, &n_d_f,f_out, NULL ); 
      return f_out[0];
  }
};



std::pair<ex,ex> RoMB_loop_by_loop::expand_and_integrate_map(ex int_in,MBintegral::w_lst_type w_lst,exmap w_curr, lst num_subs, int expansion_order) // up to O(eps^1) 
{
    try
    {
        ex out_ex;

        //            cout<<(it->get_pole_lst().subs(it->get_w())).subs(get_symbol("eps")==0)<<endl;
        //         cout<<endl <<"series:  " << it->get_gamma_expr().series(get_symbol("eps")==0,expansion_order)<<endl<<endl;

        if(w_lst.size() > 0)// expanding and integrating
        {
            //int_in.barnes1();
            // int_in.barnes2();
            
            
            /////////////////////////////////////
            // shifting contours
            /////////////////////////////////////
            
            
            
            out_ex = series_to_poly( int_in.series(get_symbol("eps"),expansion_order) ).expand().subs(num_subs);
            if (!out_ex.is_zero())
            {
             
                //out_ex = series_to_poly( int_in.get_expr().series(int_in.get_eps(),expansion_order) ).subs(num_subs);
                // loop over W_i, converting integration contour
                //          for(lst::const_iterator wit = w_lst.begin(); wit != w_lst.end(); ++wit)
           
//            cout<<"current : "<<out_ex<<endl;
                ex wo_eps_part = out_ex;
                ex vegas_ex = 0;
                ex vegas_err = 0;
                for(int i = out_ex.ldegree( get_symbol("eps") ); i <expansion_order/* out_ex.degree( get_symbol("eps") )*/; i++)
                {
                

                    cout<<"eps^ "<<i<<"  : "<<endl;//<< out_ex.coeff(get_symbol("eps"),i)<<endl;
                    ex int_expr =  out_ex.coeff(get_symbol("eps"),i);

                    cout <<" \n\n\n\t INT EXPR\n" ;


                    print_mathematica_ex(int_expr);

                    cout <<" \n\t INT EXPR\n" ;


                    lst psl;
                    exset f_gammapsi;
                    // GAMMA(Z)
                    exset fnd_tgamma;
                    int_expr.find(tgamma(wild()),fnd_tgamma);
                    //cout<<" TOSHIFT: "<<fnd.subs(w_curr)<<endl<<endl;
                    f_gammapsi.insert(fnd_tgamma.begin(), fnd_tgamma.end());
                    BOOST_FOREACH(ex Fe,fnd_tgamma)
                    {
                        exmap repls_fe;
                        if(Fe.match(tgamma(wild()),repls_fe))
                        {

                            psl.append(wild().subs(repls_fe).subs(get_symbol("eps")==0));
                        }
                    }
                    // PSI(Z)
                    exset fnd_psi;
                    int_expr.find(psi(wild()),fnd_psi);
                    //cout<<" TOSHIFT: "<<fnd.subs(w_curr)<<endl<<endl;
                    f_gammapsi.insert(fnd_psi.begin(), fnd_psi.end());
                    BOOST_FOREACH(ex Fe,fnd_psi)
                    {
                        exmap repls_fe;
                        if(Fe.match(psi(wild()),repls_fe))
                        {

                            psl.append(wild().subs(repls_fe).subs(get_symbol("eps")==0));
                        }
                    }
                    // PSI(K,Z)
                    exset fnd_psi_i;
                    int_expr.find(psi(wild(1),wild()),fnd_psi_i);
                    //cout<<" TOSHIFT: "<<fnd.subs(w_curr)<<endl<<endl;
                    f_gammapsi.insert(fnd_psi_i.begin(), fnd_psi_i.end());
                    BOOST_FOREACH(ex Fe,fnd_psi_i)
                    {
                        exmap repls_fe;
                        if(Fe.match(psi(wild(1),wild()),repls_fe))
                        {

                            psl.append(wild().subs(repls_fe).subs(get_symbol("eps")==0));
                        }
                    }

                    exmap sc_map = shift_contours(f_gammapsi,psl,w_curr);



                    lst w_for_pointer;
                    BOOST_FOREACH(ex wf,w_lst)
                    {
                        w_for_pointer.append(wf);
                        //ex c_i = wf.subs(w_curr);
                          ex c_i = wf.subs(sc_map);
                    
                        int_expr = (I*int_expr.subs(wf==c_i - I*log( wf/( 1 - wf ) ) ) ) / wf/(1- wf);
                    }


///////////////////////////
                    symbol x1("x");
                    ex myex = sin(x1)/x1;
                
                    RoMB::FUNCP_CUBA2 fp;
                    RoMB::compile_ex_real(lst(myex), lst(x1), fp);
                    int nn = 0, nb = 0;
                    double to_c1[] = {0};
                    double to_c2[] = {0.2};
                    double to_c3[] = {0.3};
                    double to_c4[] = {0.4};
                    double to_c5[] = {0.5};
                    double to_c6[] = {0.6};
                    double to_c7[] = {0.7};
                    double to_c8[] = {0.8};
                    double to_c9[] = {1};
                    double f_c[1];
                    fp(&nn,to_c1,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;
                
///////////////////////////

                
                    RoMB::FUNCP_CUBA2 fp_real;
                    std::string int_c_f(boost::filesystem::current_path().string());
                    int_c_f+="/int_c_f";
                    cout<<"\n\n\n"<<evalf(int_expr)<<"\n\n\n\n"<<endl;

                    print_mathematica_ex(evalf(int_expr));

                    cout<< w_for_pointer << endl;
                    RoMB::compile_ex_real(lst(evalf(int_expr)),w_for_pointer, fp_real);//,int_c_f);
                    fp_real(&nn,to_c1,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;

                    fp_real(&nn,to_c2,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;

                    fp_real(&nn,to_c3,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;

                    fp_real(&nn,to_c4,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;

                    fp_real(&nn,to_c5,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;

                    fp_real(&nn,to_c6,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;

                    fp_real(&nn,to_c7,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;

                    fp_real(&nn,to_c8,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;

                    fp_real(&nn,to_c9,&nb,f_c,0);
                    std::cout << f_c[0] << std::endl;


                    // ----------------------------------- Vegas integration-------------------------
                    const int MAXCUHREDIM = 4;
                    int  NDIM  = w_lst.size();
                    //#define NCOMP 1
#define USERDATA___ NULL
#define EPSREL___ 1e-4
#define EPSABS___ 1e-12
#define VERBOSE___ 0
#define LAST___ 4
#define SEED___ 0
#define MINEVAL___ 0
#define MAXEVAL___ 1000000

#define NSTART___ 1000
#define NINCREASE___ 1000
#define NBATCH___ 1000
#define GRIDNO___ 0
#define STATEFILE___ NULL

#define KEY___ 9
                    const int NCOMP = 1;
                    int comp,  neval, fail,nregions;
                    double integral_real[NCOMP], error[NCOMP], prob[NCOMP];
                    if (NDIM == 1)
                    {
                        printf("-------------------- QAGS  --------------------\n");
                        gsl_integration_workspace * w 
                            = gsl_integration_workspace_alloc (1000);
                        double result1d, error1d;
                    
                        gsl_function F;
                        //                    F.function = &f;
                        F.function = GSL1dintAdapter<RoMB::FUNCP_CUBA2>::f;
                        F.params = reinterpret_cast<void *>(fp_real);
                    
                        gsl_integration_qags (&F, 0, 1, 1e-12, 1e-7, 1000,
                                              w, &result1d, &error1d); 
                    
                        printf("QAGS RESULT:\t%.8f +- %.8f\t intervals = %d\n",
                               result1d, error1d, int(w->size));
                    
                        gsl_integration_workspace_free (w);
                        vegas_ex+=pow(get_symbol("eps"),i)*(result1d);//+I*integral_imag[0]);
                        vegas_err+=pow(get_symbol("eps"),i)*(error1d);//+I*integral_imag[0]);
                    }
                    else if(NDIM > MAXCUHREDIM )
                    { 

                        printf("-------------------- Vegas  --------------------\n");
                        Vegas(NDIM, NCOMP, fp_real, USERDATA___,
                              EPSREL___, EPSABS___, VERBOSE___, SEED___,
                              MINEVAL___, MAXEVAL___, 
                              NSTART___, NINCREASE___, NBATCH___, GRIDNO___, STATEFILE___,
                              &neval, &fail, integral_real, error, prob);
                    
                    
                    
                        // printf("VEGAS RESULT:\tneval %d\tfail %d\n",
                    
                        for( comp = 0; comp < NCOMP; ++comp )
                            printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                                   integral_real[comp], error[comp], prob[comp]);
                    
                        // ------------------------------ Vegas integration--------------------              
                        vegas_ex+=pow(get_symbol("eps"),i)*(integral_real[0]);//+I*integral_imag[0]);
                        vegas_err+=pow(get_symbol("eps"),i)*(error[0]);//+I*integral_imag[0]);
                    }
                    else
                    {
                        printf("-------------------- Cuhre  --------------------\n");
                        Cuhre(NDIM, NCOMP, fp_real, USERDATA___,
                              EPSREL___, EPSABS___, VERBOSE___,
                              MINEVAL___, MAXEVAL___, KEY___,
                              &nregions, &neval, &fail, integral_real, error, prob);
                   
                        for( comp = 0; comp < NCOMP; ++comp )
                            printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
                                   integral_real[comp], error[comp], prob[comp]);
                        cout<<"Hi there"<<endl;

                        // ----------------------------------- Cuhre integration-------------------------              
                        vegas_ex+=pow(get_symbol("eps"),i)*(integral_real[0]);//+I*integral_imag[0]);
                        vegas_err+=pow(get_symbol("eps"),i)*(error[0]);//+I*integral_imag[0]);
                        
                    }
                    cout<< " END int " << endl;
                
                }
                return std::make_pair(vegas_ex,vegas_err);
            }
            else return std::make_pair(0,0);
        }
        else // expanding only
        {
            return  std::make_pair(series_to_poly( int_in.series(get_symbol("eps"),expansion_order) ).subs(num_subs), 0);
        }
    }catch(std::exception &p)
    {
        throw std::logic_error(std::string("In function \"Expand_and_integrate_map\":\n |___> ")+p.what());
    }
}


NearestPoleParams RoMB_loop_by_loop::GetLeadingEps(MBintegral mbIn, numeric epsCurrent, numeric eps0)
{
    try
    {

        cout << "\n\t Get Nearest Eps for CurrentEps:\t" << epsCurrent <<"("<<evalf(epsCurrent)<<")"<< endl; 
        
        NearestPoleParams poleParams;

        poleParams.isContinued = true;
    
        if(epsCurrent == eps0) 
            throw std::logic_error( string("Continue from EpsI = Eps0") );
        else
          poleParams.EpsilonValue = eps0;
        lst polesWithEps(mbIn.poles_with_eps());
    
        for(lst::const_iterator pit  = polesWithEps.begin(); pit != polesWithEps.end(); ++pit) 
        {
        
            numeric poleValue;
        
            numeric fEps0 = ex_to<numeric>(pit->subs( this->constraints_.GetWs() ).subs( get_symbol("eps") == eps0));
            numeric fEpsI = ex_to<numeric>(pit->subs( this->constraints_.GetWs() ).subs( get_symbol("eps") == epsCurrent ));
        
            cout << "Fi\t" << fEpsI << " (" << evalf(fEpsI)<< ")" <<endl;
            cout << "F0\t" << fEps0 << " (" << evalf(fEps0)<< ")" <<endl;

	    // direction of movement of Gamma arguement 
            int dir = csgn(fEps0 - fEpsI);



// Important FF

            bool hasPole = false;

	    if ( ( dir > 0 &&  fEpsI > 0 ) || ( (dir < 0) && ( fEpsI > 0) && (fEps0 > 0) )) hasPole = false; 
	    if ( (dir < 0) && ( fEpsI > 0) && (fEps0 < 0) ) { hasPole = true; poleValue = 0; }
	    if ( fEpsI < 0  && dir > 0 && ceil(fEpsI.to_double()) < fEps0) { hasPole = true; poleValue = ceil(fEpsI.to_double()); }
	    if ( fEpsI < 0  && dir < 0 && floor(fEpsI.to_double()) > fEps0) { hasPole = true; poleValue = floor(fEpsI.to_double()); }

            if (hasPole && (poleValue <= 0)) 
            {
                   cout << setw(5) << right << poleValue <<  " --------------->  " 
                     << setw(25) << left << *pit << endl;
                numeric eps_pole_sol = ex_to<numeric>(lsolve(pit->subs(this->constraints_.GetWs()) == poleValue,get_symbol("eps") ));

                cout << "SOlution: " << eps_pole_sol << endl;

                  if(abs(epsCurrent - poleParams.EpsilonValue) > abs(epsCurrent - eps_pole_sol))
                    {
                        poleParams.EpsilonValue = eps_pole_sol;
                        poleParams.PoleValue = poleValue;
                        poleParams.Arg = (*pit);
                        poleParams.isContinued = false;
                    }
            }
        }
    
        poleParams.Print(); 
        return poleParams;
    }
    catch(std::exception &p)
    {
        throw std::logic_error(std::string("In function \"GetLeadingEps\":\n |___> ")+p.what());
    }
}





MBlst RoMB_loop_by_loop::MBcontinue(MBintegral rootint,ex eps0)
{
    try
    {
        // integrals to be continued
        MBlst O(1,rootint);
        
        // final set of continued integrals
        MBlst C;
        size_t nRej = 0;

        while(O.size()>0)
        {
            
            // integral from residues
            MBlst R;
            
            for(MBlst::iterator it = O.begin();it!=O.end();++it)
            {
                C.push_back(*it);


//         Eps starting value

                ex eps_i = ex(get_symbol("eps")).subs(it->get_eps());

                
//         Sliding eps value

                ex epsSliding = eps_i;


//         Global list of (*it) poles

                lst polesWithEps(it->poles_with_eps());

                bool mayBeContinued = true;


//         Iterate over all residues 
//         of current integral (*it)

                while(mayBeContinued)
                {




//         Find nearest pole to continue, return EPS in "EpsilonValue"
//         and gamma argument in "PoleValue"
//         if no pole exists - return                  

                  cout << ">> Step in Loop ____________________________________" << endl;
//         Nearest pole               

                  cout << ">> Step in Loop _______(sliding eps)____________" << epsSliding << endl;

          cout<< endl<<" Ready for continue?  [Y/n]: ";
          char in_ch;
          std::cin>>in_ch;

                    NearestPoleParams nearestPoleParams = GetLeadingEps(*it, ex_to<numeric>(epsSliding), ex_to<numeric>(eps0));

                    nearestPoleParams.Print();

// IC =false - pole exists and returned.


//
// Place test on reduce HERE
//
                    
                    if(!nearestPoleParams.isContinued) // test on pole existance
                    {

                      cout << "\tPOLE MOVED====================================" << endl;                    
                      cout << "\tPOLE MOVED========was  "<< constraints_.GetPoint() <<  endl;                    
                        bool isRejected = constraints_.Restrict(nearestPoleParams);
                        if (isRejected) nRej++;
                      cout << "\tPOLE MOVED=======become ==== " << nRej << " ==========" << constraints_.GetPoint() <<  endl;                    
                        
                        ex poleF =  nearestPoleParams.Arg;
                        
                        lst w_in_F  = it->has_w(poleF );
                    
                        //epsSliding = nearestPoleParams.EpsilonValue;
                        
                        epsSliding = get_symbol("eps");
                        epsSliding = epsSliding.subs(constraints_.GetPoint());
                        
//                        mayBeContinued = nearestPoleParams.isContinued;
                        mayBeContinued = true;
                    
                        cout << "EpsSl: " << epsSliding <<endl;

                        cout << "TEST IF: " << endl;
                        cout << "\t " << w_in_F.nops() << " " << nearestPoleParams.isContinued << " " << isRejected <<endl;
                    
                        if(w_in_F.nops()>0 && !nearestPoleParams.isContinued && !isRejected) 
                        {


//                              
                            //             decide what var to get res
                            ex var_to_get_res = 0;
                            for(lst::const_iterator vgit = w_in_F.begin(); vgit != w_in_F.end();++vgit) 
                            {
                                if(poleF.coeff(*vgit,1) == 1) 
                                {
                                    var_to_get_res = *vgit;
                                    break;
                                }
                                if(poleF.coeff(*vgit,1) == -1)
                                    var_to_get_res = *vgit;
                            }
                            if( var_to_get_res ==0) var_to_get_res = w_in_F.op(w_in_F.nops()-1);
                        
                            cout<<" POLE: " << poleF<< "       var to get res   "<<var_to_get_res<<endl;
                        
                            ex fEps0 = poleF.subs(constraints_.GetWs()).subs(get_symbol("eps") == eps0);
                            ex fEpsI = poleF.subs(constraints_.GetWs()).subs(get_symbol("eps") == epsSliding);
    
                            MBintegral res_int = 
                                it->res(var_to_get_res ==
                                        lsolve(poleF == nearestPoleParams.PoleValue, var_to_get_res),
                                        poleF,get_symbol("eps")==epsSliding);
                        
                            res_int.set_level(1+it->get_level());
                            res_int*=(2*Pi*I*csgn(poleF.coeff(var_to_get_res))*csgn(fEpsI - fEps0));
                        
                                             
                            R.push_back(res_int);
                        
                            mayBeContinued = true;
                        }
                    
                    }// if may be continued
                    else mayBeContinued = false;
                } // WHILE(mayBeContinued)
                
            }   
            
            // continue integrals only from the last level of residue tree
            O = R;
        }
        cout << "REJ= " << nRej <<endl;
        return C;
    } 
    catch(std::exception &p)
    {
        throw std::logic_error(std::string("In function \"MBcontinue\":\n |___> ")+p.what());
    }
    
}
    

MBtree RoMB_loop_by_loop::MBcontinue_tree(MBintegral rootint,ex eps0)
{
    using namespace mbtree;
    try 
    {
        rootint.barnes1();
        rootint.barnes2();
        
        //           accumulator for constraints
        //  cout<<"start constrainta : "<<rootint.get_poles_set()<< " "<<rootint.get_w_eps_set()<<endl;
//        ConstrAcc ca(rootint.get_poles(),rootint.get_w_eps());

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

                    cout<< " $$$" <<endl;
                    cout<< " $$$" <<endl;
                    cout<< " $$$" <<endl;
                    cout<< " $$$" <<endl;
                    cout<< " $$$" <<endl;
                    cout<< " $$$" <<endl;
                    cout<< " $$$$$$$$" <<endl;
                    cout<< " $$$$$$$$" <<endl;
                    //   cout<<std::setw(15+it->get_level())<<std::right<<"shifted on "<<it->get_level()<<endl;
                    //C.push_back(*it);//need review, multiple entries C=C U I
                    MBintegral::pole_iterator pit,pit_end;
                    ex eps_i = get_symbol("eps");
                    //          cout<<"after barness lemas "<<it->get_eps()<<endl;
                    eps_i = eps_i.subs(it->get_eps());


                    ex epsSliding = eps_i;

                    bool mayBeContinued = true;

                    // finding nearest EPS

                    typedef multimap<ex,ex> EpsPolesMap;
                    typedef EpsPolesMap::iterator mapIter;


                    lst polesWithEps(it->poles_with_eps());
                      
                    
                    while(mayBeContinued)
                    {

                        // By defualt if no residues added integral can't be continued
                        mayBeContinued = false;



                        EpsPolesMap epsPoles;

                        //             Iterate over gamma arguments with eps dependence only!!!!!!!
                        lst polesMinEps;
                        for(lst::const_iterator pit  = polesWithEps.begin(); pit != polesWithEps.end(); ++pit) 
                        {

                            ex poleValue;

                            ex fEps0 = pit->subs( this->constraints_.GetWs() ).subs( get_symbol("eps") == eps0) ;
                            ex fEpsI =  pit->subs( this->constraints_.GetWs() ).subs( get_symbol("eps") == epsSliding ) ;

                            if(fEpsI > fEps0) poleValue = int(floor(ex_to<numeric>(fEpsI).to_double()));
                            else if(fEpsI < fEps0) poleValue = int(ceil(ex_to<numeric>(fEpsI).to_double()));
                            if (poleValue <= 0) 
                            {
                                cout << setw(5) << right << poleValue <<  " --------------->  " 
                                     << setw(25) << left << *pit << endl;
                                ex eps_pole_sol = lsolve(pit->subs(this->constraints_.GetWs()) == poleValue,get_symbol("eps") );
                                polesMinEps.append(eps_pole_sol);
                                epsPoles.insert ( pair<ex,ex>(eps_pole_sol,*pit) );
                            }
                        }

                        cout << " ______________POLES NEAREST_______________  " << endl;
                        cout << "|                                          | " << endl;
                        cout << "| " << setw(40) << std::left <<S(this->constraints_.GetWs()) << endl;
                        cout << "|        " << polesMinEps << endl;
                        cout << "|__________________________________________| " << endl;
                        cout << endl;//"*** min elem: "<<*min_element(polesMinEps.begin(),polesMinEps.end())<<endl;
                    
                        mapIter m_it, s_it;
                    
                        for (m_it = epsPoles.begin();  m_it != epsPoles.end();  m_it = s_it) 
                        {
                            ex theKey = (*m_it).first;
                        
                            cout << endl;
                            cout << "  key = '" << theKey << "'" << endl;

// Storing max eps to Sliding Eps
                            epsSliding = theKey;
                            
                            ex poleF =  m_it->second.op(0);
                            
                            lst w_in_F  = it->has_w(poleF );
                            
                            if(w_in_F.nops()>0) 
                            {
//                                   cout<<endl<<endl<<endl<<"ASDASDASDASDSD "<<*pit<<" "<<ca.test_single(pit->subs(get_symbol("eps") == 0))<<endl<<endl<<endl;
                                
//                                cout<<endl<<"LEVEL "<<it->get_level()<<" Epsilon continue from eps_i = "<<eps_i<<" to "<<eps_prime<<endl<<endl;
//                                BOOST_ASSERT_MSG(abs(ex_to<numeric>(eps_i).to_double())>=abs(ex_to<numeric>(eps_prime).to_double()), "Bad continuation");
                                
                                //             decide what var to get res
                                ex var_to_get_res = 0;
                                for(lst::const_iterator vgit = w_in_F.begin(); vgit != w_in_F.end();++vgit) 
                                {
                                    if(poleF.coeff(*vgit,1) == 1) 
                                    {
                                        var_to_get_res = *vgit;
                                        break;
                                    }
                                    if(poleF.coeff(*vgit,1) == -1)
                                        var_to_get_res = *vgit;
                                }
                                if( var_to_get_res ==0) var_to_get_res = w_in_F.op(w_in_F.nops()-1);
                                    
                                cout<<" POLE: " << poleF<< "       var to get res   "<<var_to_get_res<<endl;

                                MBintegral res_int = it->res(var_to_get_res ==lsolve(poleF==poleF.subs(constraints_.GetWs()).subs(get_symbol("eps") == epsSliding), var_to_get_res),poleF,get_symbol("eps")==epsSliding);
                                res_int.set_level(1+it->get_level());
                                res_int*=(2*Pi*I*csgn(poleF.coeff(var_to_get_res))*csgn(epsSliding-poleF.subs(constraints_.GetWs()).subs(get_symbol("eps") == 0)));
                                //if(ca.test_single(poleF.subs(get_symbol("eps") == 0)))res_int.set_optimizable(true);
                                res_int.barnes1();
                                res_int.barnes2();
                                //R.push_back(res_int);
                                //  if(eps_prime !=0)
                                lastChildIt = C.append_child(it,res_int);
                                nChildrenAdded++;
                            }


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
                        
                        if (polesMinEps.nops() == 0) mayBeContinued = false;
                        else cout << "Poles MIN eps for test: " << polesMinEps << endl << endl;

                    } // WHILE(mayBeContinued)

/*
  LOOP over all gamma arguments to find nearest poles
*/
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


