#include "romb.h"
#include "utils.h"
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

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

          if(MBlbl_int.w_size()>0)          
            {
              MBlbl_int.new_point();
              w_shared = MBlbl_int.get_w();
            }
          print_mathematica(MBlbl_int);
          // MB only if w's existsts or eps!=0
          if((MBlbl_int.w_size() > 0) && (MBlbl_int.get_eps().rhs() != 0))
          {
            
            
          // setting shared point
         


          cout<< endl<<" Ready for MBcontinue?  [Y/n]: ";
          char in_ch;
          std::cin>>in_ch;
          if(in_ch=='n')  exit(0);//assert(false);



          //                assert(false);

          //        MBlst int_lst = MBcontinue(MBlbl_int);
          MBtree inttr = MBcontinue_tree(MBlbl_int);
          int opt_sum = 0;
          MBtree::breadth_first_queued_iterator last_it;
          for(MBtree::breadth_first_queued_iterator bfit = inttr.begin_breadth_first(); bfit != inttr.end_breadth_first(); ++bfit)
            {
              if(bfit->get_optimizable())
                {
                  cout<<"SIZE: "<<inttr.size(bfit)<<endl;
                  last_it  = bfit;
                  //inttr.erase(bfit);
                  cout<< inttr.size()<<endl;
                  //break;   
		
                  //opt_sum +=inttr.size(bfit);
                }
            }
          /*     
                 MBtree::sibling_iterator sit_gen(last_it);
                 // loop over neighbours
                 std::vector<MBtree::iterator> v_it;
                 for(MBtree::iterator sit = sit_gen.range_first(); sit != sit_gen.range_last(); ++ sit)
                 {
                 if(sit->get_optimizable())
                 {
                 v_it.push_back(sit);
                 }
                 }

                 //combination code
                 const int r = 3;
                 const int n = 10;
                 std::vector<int> v_int(n);
                 for (int i = 0; i < n; ++i) { v_int[i] = i; }
                 int N = 0;
                 do {
                 ++N;
                 if (N < 10 || N > 117) {
                 std::cout << "[ " << v_int[0];
                 for (int j = 1; j < r; ++j) { std::cout << ", " << v_int[j]; }
                 std::cout << " ]" << std::endl;
                 } else if (N == 10) {
                 std::cout << "  . . ." << std::endl;
                 }
                 } while (next_combination(v_int.begin(), v_int.begin() + r, v_int.end()));
                 // end of combination code
                 */
          cout<< "Min int: "<<inttr.size()-opt_sum<<endl;
          int_lst.assign(inttr.begin(),inttr.end());
          //exit(5);//assert(false);
          //int_lst = MBcontinue(MBlbl_int);
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
          cout<<"Integrals MBT "<<inttr.size()<<endl;
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
        
          // MBlst int_lst = MBcontinue(Uint);
          //ex int_expr_out = 0;
          //for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
          // int_expr_out+=expand_and_integrate(*it,1);
        
          //cout<<" RESULT : "<<endl
          //    <<"               = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
          }
          else // no contour integration or no continuation in eps
          {
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
      if(e.match(tgamma(wild()),repls))
        {
          stringstream str;
          str<<"Gamma["<<wild().subs(repls)<<"]";
          return get_symbol(str.str());
        }
      else if(e.match(exp(wild()),repls))
        {
          stringstream str;
          str<<"Exp["<<wild().subs(repls)<<"]";
          return get_symbol(str.str());
        }
      else return e.map(*this);
    }
  };

  void print_mathematica(MBintegral mb_in)
  {
    exmap w_c(mb_in.get_w());
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
  void RoMB_loop_by_loop::integrate(lst number_subs_list, int exp_order)
  {
    ex int_expr_out = 0;
    for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
      int_expr_out+=expand_and_integrate(*it, number_subs_list, exp_order);
    cout<<" FRESULT for parameters: "<<number_subs_list<<endl<<endl;
    cout<<" FRESULT anl : "<<"          = "<<int_expr_out.expand().collect(get_symbol( "eps" ))<<endl;
    cout<<" FRESULT num: "<<"          = "<<evalf(int_expr_out.expand().collect(get_symbol( "eps" )))<<endl;
  }


  void RoMB_loop_by_loop::merge()
  {
    for(MBlst::iterator it = int_lst.begin();it!= int_lst.end();++it)
      {
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

  std::pair<ex,ex> expand_and_integrate_map(ex int_in,MBintegral::w_lst_type w_lst,exmap w_curr, lst num_subs, int expansion_order) // up to O(eps^1) 
  {
    try
      {
        ex out_ex;

        //            cout<<(it->get_pole_lst().subs(it->get_w())).subs(get_symbol("eps")==0)<<endl;
        //         cout<<endl <<"series:  " << it->get_gamma_expr().series(get_symbol("eps")==0,expansion_order)<<endl<<endl;

        if(w_lst.size() > 0)// expanding and integrating
          {
            //          int_in.barnes1();
            // int_in.barnes2();
            out_ex = series_to_poly( int_in.series(get_symbol("eps"),expansion_order) ).expand().subs(num_subs);
            //out_ex = series_to_poly( int_in.get_expr().series(int_in.get_eps(),expansion_order) ).subs(num_subs);
            // loop over W_i, converting integration contour
            //          for(lst::const_iterator wit = w_lst.begin(); wit != w_lst.end(); ++wit)
            lst w_for_pointer;
            BOOST_FOREACH(ex wf,w_lst)
              {
                w_for_pointer.append(wf);
                ex c_i = wf.subs(w_curr);
                out_ex = (I*out_ex.subs(wf==c_i - I*log( wf/( 1 - wf ) ) ) ) / wf/(1- wf);
              }
            //cout<<"current : "<<out_ex<<endl;
            ex wo_eps_part = out_ex;
            ex vegas_ex = 0;
            ex vegas_err = 0;
            for(int i = out_ex.ldegree( get_symbol("eps") ); i <expansion_order/* out_ex.degree( get_symbol("eps") )*/; i++)
              {
                cout<<"eps^ "<<i<<"  : "<<endl;//<< out_ex.coeff(get_symbol("eps"),i)<<endl;
                ex int_expr =  out_ex.coeff(get_symbol("eps"),i);
                RoMB::FUNCP_CUBA2 fp_real;
                std::string int_c_f(boost::filesystem::current_path().string());
                int_c_f+="/int_c_f";
                RoMB::compile_ex_real(lst(evalf(int_expr)),w_for_pointer, fp_real);//,int_c_f);

                // ----------------------------------- Vegas integration-------------------------
                int  NDIM  = w_lst.size();
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
                vegas_err+=pow(get_symbol("eps"),i)*(error[0]);//+I*integral_imag[0]);
              }
            return std::make_pair(vegas_ex,vegas_err);
          }
        else // expanding only
          {
            return  std::make_pair(series_to_poly( int_in.series(get_symbol("eps"),expansion_order) ).subs(num_subs), 0);
          }
      }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"Expand_and_integrate\":\n |___> ")+p.what());
      }
  }
