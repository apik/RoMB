#include "romb.h"
#include "utils.h"
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
				      lst nu)
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
      exlist input_prop_set( p_lst.begin(),p_lst.end());
      cout<<"INPSET: "<<input_prop_set<<endl;

      /* 
	 map for propagator powers
      */
      exmap prop_pow_map;
      for(lst::const_iterator Pit = p_lst.begin(); Pit != p_lst.end(); ++Pit)
	{
	  prop_pow_map[*Pit] = nu.op(std::distance(p_lst.begin(),Pit));
	}
      /* 
	 Iterate over momentums k1,k2,k3,etc.
      */
      unsigned int displacement_x = 0;
      unsigned int displacement_w = 0;
      for(lst::const_iterator kit = k_lst.begin(); kit != k_lst.end(); ++kit)
	{
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
            
          ex mul_int(1);
          //          exmap mul_subs_map;
          lst coe_prop_lst;
          for(lst::const_iterator kkit = P_with_k_lst.begin(); kkit != P_with_k_lst.end(); ++kkit)
            {
              if(kkit->expand().coeff(*kit,2) <0 )
                {
                  coe_prop_lst.append((-1)*(*kkit));
                  mul_int *= pow(-1,prop_pow_map[*kkit]);
                }
              else
                coe_prop_lst.append(*kkit);
            }
	  cout<< "Set wo k_i "<<input_prop_set<<endl;
	  cout<<" PWKlst "<<P_with_k_lst<<endl;
          cout<< " coe: "<<coe_prop_lst<<endl;
	  /*
	    lexi sort of input prop list, and it's modification
	  */            
 

	  // uf and then MB represenatation construction
	  // subs only in F for last momentum
	  UFXmap inUFmap;
	  if(boost::next(kit) == k_lst.end())
            inUFmap = UF(lst(*kit),coe_prop_lst,subs_lst,displacement_x);
	  else
            inUFmap = UF(lst(*kit),coe_prop_lst,subs_lst,displacement_x); // no substitution!!!
	  displacement_x +=fusion::at_key<UFX::xlst>(inUFmap).nops(); 

	  lst nu_into;
	  for(lst::const_iterator nuit = P_with_k_lst.begin(); nuit != P_with_k_lst.end(); ++nuit )
	    nu_into.append(prop_pow_map[*nuit]);
	  cout<<" Powers list before input: "<<nu_into<<endl;
	  MBintegral Uint(
			  fusion::make_map<UFX::F,UFX::xlst>(fusion::at_key<UFX::F>(inUFmap),
							     fusion::at_key<UFX::xlst>(inUFmap)
							     ),nu_into,1,displacement_w);
	  displacement_w+=Uint.w_size();
	  cout<<"ui9nt eps : "<<Uint.get_expr()<<endl;
          Uint *= mul_int;
	  /*
	    expression to mul root integral
	    where to subs prop(k_prev)==1
	  */
	  ex expr_k_to_subs_1= Uint.get_expr();

	  ex mom_find = Uint.get_expr();
	  cout<< "where find props: "<<mom_find<<endl;
	  if(is_a<mul>(mom_find))
	    {
	      exset found_prop;
	      mom_find.find(pow(wild(1),wild(2)),found_prop);
	      cout<<" is a mul  "<<found_prop<<endl;
	      //                lst prop_pow_lst;
	      if(boost::next(kit)!=k_lst.end())
		{
		  ex mom_to_find = *(boost::next(kit));
		  BOOST_FOREACH(ex propex,found_prop)
		    {
		      cout<<"kit "<<mom_to_find<<" propex "<<propex.has(mom_to_find)<<endl;                        
                        
		      if(propex.has(mom_to_find))
			{

			  cout<<"before subs kex : "<<expr_k_to_subs_1<<endl;
			  expr_k_to_subs_1 = expr_k_to_subs_1.subs(propex == 1);
			  cout<<"after subs kex : "<<expr_k_to_subs_1<<endl;
			  /*
			    Search for duplications in prop set
			  */
			  cout<<"where to find props: "<<input_prop_set<<endl;
			  bool dupl_found = false;
			  cout<< input_prop_set.size()<<endl;
			  BOOST_FOREACH(ex pr, input_prop_set)
			    {
			      exmap repls;
			      if(is_a<power>(propex))
				{
				  if((ex_to<power>(propex).op(0)).match(wild()*pr,repls))
				    {
				      cout<<"MATCH WILD "<<repls<<endl;
				      ex mul_ex = pow(wild(),ex_to<power>(propex).op(1)).subs(repls);
				      cout<<"MUL_EX: "<<mul_ex<<endl;
				      expr_k_to_subs_1*=mul_ex;

				      prop_pow_map[pr] += (-1)*ex_to<power>(propex).op(1);
				      dupl_found = true;
				    }
				  if(ex_to<power>(propex).op(0).match(pr))
				    {
				      cout<<"add"<<endl;
				      prop_pow_map[pr] += (-1)*ex_to<power>(propex).op(1);
				      dupl_found = true;
				    }
				}
			      else  throw std::logic_error(string("incompitible power"));
			    }
			  if(!dupl_found)
			    {
			      if(is_a<power>(propex))
				input_prop_set.push_back(ex_to<power>(propex).op(0));
			      else throw std::logic_error(string("incompitible power"));
			      //prop_pow_lst.append(ex_to<power>(propex).op(1));
			      prop_pow_map[ex_to<power>(propex).op(0)] = (-1)*ex_to<power>(propex).op(1);
			    }
			}
		    }
		  //cout<<"needed props "<<prop_pow_lst<<endl;
		}
	    }
            


	  MBlbl_int*=expr_k_to_subs_1;

	  MBlbl_int+=Uint;
	  cout<<"bad"<<endl;
	  //            MBlbl_int.insert_w_lst(Uint.get_w_lst());
	  // MBlbl_int.insert_pole_lst(Uint.get_pole_lst());
            
            
	}

      MBlbl_int.fix_inv();
      cout<<"Constructed integral with:"<<endl;
      //cout<<"Poles: "<<MBlbl_int.get_poles_set()<<endl;
      //MBlbl_int.set_poles_set(MBlbl_int.poles_from_ex(MBlbl_int.get_expr()));
    cout<<"Poles true: "<<MBlbl_int.poles_from_ex(MBlbl_int.get_expr())<<endl;
    cout<<"Poles: "<<MBlbl_int.get_poles()<<endl;
    cout<<"W's : "<<MBlbl_int.get_w_lst()<<endl;
      cout<<"Expr : "<<MBlbl_int.get_expr()<<endl;

      cout<< endl<<" Ready for MBcontinue?  [Y/n]: ";
      char in_ch;
      std::cin>>in_ch;
      if(in_ch=='n')  exit(0);//assert(false);
      MBlbl_int.new_point();        


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
    }catch(std::exception &p)
    {
      throw std::logic_error(std::string("In function \"RoMB\":\n |___> ")+p.what());
    }
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
