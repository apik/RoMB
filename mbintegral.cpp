#include "mbintegral.h"
#include "utils.h"
#include "constants.h"
/**
 *
 *  Construction from F, U=1
 *
 */

MBintegral::MBintegral(FXmap fx_in,lst nu,numeric l, unsigned int displacement):tree_level(0) // lst nu is a list of powers of propagators and l is a number of loops
  {
    try
      {
        eps_current = (get_symbol("eps")==0);
        ex N = accumulate(nu.begin(),nu.end(),ex(0));

        //    cout<<"N(nu)= "<<N<<endl;
        lst  x_lst(fusion::at_key<UFX::xlst>(fx_in));
        // summ of feynman parameters satisfying summ(X_i)=1
        //        ex x_summ = accumulate(x_lst.begin(),x_lst.end(),ex(0));

        //relational delta_subs(x_lst.op(x_lst.nops()-1),lsolve(x_summ==1,x_lst.op(x_lst.nops()-1)));// relation delta function
        // cout<<"xsumm "<<delta_subs.lhs()<<" == "<<delta_subs.rhs()<<endl;
        
        //ex F = fusion::at_key<UFX::F>(fx_in).collect(x_lst,true);///< distributed polynom factorizing X(i)X(j)*(...)
        ex F = (fusion::at_key<UFX::F>(fx_in)).expand();///< distributed polynom factorizing X(i)X(j)*(...)
        //F=F.subs(delta_subs);//applying delta function relation on F-polynom
        //x_lst.remove_last();

        // assuming 1/(U^a*F^b)
        
        ex F_pow = (N-l*(2-get_symbol("eps")));
	exset f_set;
	if(F.find(wild(1)*wild(1)+2*wild(1)*wild(2)+wild(2)*wild(2), f_set))cout<< "HAVE square"<<endl;
        cout<<setw(30)<<std::internal<<"+++INTEGRAL PARAMETERS+++"<<endl;
        cout<<setw(30)<<std::left<<"** Number of loops L=   "<<std::internal<<l<<endl;
        cout<<setw(30)<<std::left<<"** Summ of powers N=   "<<N<<endl;
        //cout<<setw(35)<<std::left<<"** Dimension D=   "<<D.subs(D_subs)<<endl;
        cout<<setw(35)<<std::left<<"** F =   "<<F<<endl;
        cout<<setw(35)<<std::left<<"** F power=   "<<F_pow<<endl;
        cout<<setw(30)<<std::internal<<"+++++++++++++++++++++++++"<<endl;
        
        	lst coe_l,xsq_l;
        	tie(coe_l,xsq_l) = collect_square(F,x_lst);
        	cout<<">>> Found "<<coe_l.nops()<<" full squares in F polynomial"<<endl;
        	cout<< coe_l<<" * "<<xsq_l<<endl;
        	cout<<" F qad " <<F<<endl;
        	
        	F = F.collect(x_lst,true);
		//	assert(false);

        //lst gamma_lst; //Gamma with poles left of contour
        //    lst w_lst;
        //  ex gamma_right; // and right of contour
        exmap x_power_map; // map of X(j) powers
        for(lst::const_iterator it1=x_lst.begin();it1!=x_lst.end();++it1)
          {
            x_power_map[*it1] = nu.op(std::distance(x_lst.begin(),it1))-1;
          }

        cout<<"x_map_start "<<x_power_map<<endl;

		
        //    ex coeff = 1;                             // numerical coeeficient independent of X(j)
        ex coeff =tgamma(F_pow)* pow(exp(get_symbol("eps")*Euler),l);///pow(I*pow(Pi,2 - get_symbol("eps")),l);
        // important if power = 0????
        for(lst::const_iterator nui = nu.begin();nui!=nu.end();++nui)
          coeff/=tgamma(*nui);
        coeff/=tgamma(F_pow);
        //    ex out_ex = 1;///pow(2*Pi*I,U.nops()-1)/tgamma(al_pow)/pow(2*Pi*I,F.nops()-1)/tgamma(al_pow); //need review

        //working with F-term \Gamma(\nu-L*D/2) contractedx
        ex w_sum = 0;  //F-term generates only integrations in W

	//--------------------------------
	//      MB for full squares
	//--------------------------------
	size_t z_idx = 0;
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
		string str = "w"+boost::lexical_cast<string>(displacement);
		displacement++;
                //symbol w(str);
		ex w_i = get_symbol(str);
                w_lst.append(w_i);
                coeff*=tgamma(-w_i)*pow(coe_l.op(i),w_i)/(2*Pi*I);
                
                cout<<"w_i_power "<<w_i<<endl;
                gamma_poles.append(-w_i); //!!!! review
                w_sum+=w_i;

		// subMB construction
		ex sq_lst(xsq_l.op(i));
		cout<<"SQLST : "<<sq_lst<<endl;
		coeff /= (pow(2*Pi*I,sq_lst.nops()-1)*tgamma(-2*w_i));
		gamma_poles.append(-2*w_i);
		ex z_sum = 0;
		for(const_iterator x_it = sq_lst.begin(); x_it != sq_lst.end(); ++x_it)
		  {
		    ex a_power;
		    if(sq_lst.end() == boost::next(x_it)) // X_k expr
		      {
			coeff *= tgamma(-2*w_i + z_sum);
			a_power = 2*w_i - z_sum;
			gamma_poles.append(-2*w_i + z_sum);
		      }
		    else // ordinary exprs
		      {
			string z_str = "z_"+boost::lexical_cast<string>(z_idx);
			z_idx++;
			ex z_i = get_symbol(z_str);
			cout << z_i<<endl;
			w_lst.append(z_i);
			coeff*=tgamma(-z_i);
			a_power = z_i;
			gamma_poles.append(-z_i); //!!!! review
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
		      }
		  }
	      }
	  }
        // x_lst -> to list<symbol> for comparision
        std::list<symbol> x_sym_list;
        for(lst::const_iterator sit = x_lst.begin(); sit != x_lst.end(); ++sit)
          if(is_a<symbol>(*sit)) x_sym_list.push_back(ex_to<symbol>(*sit));
       
        ex_is_lesseq_degrevlex F_comp(x_sym_list);


        exlist Fl(F.begin(),F.end());
        //   cout<<"FEX "<<Fex<<endl;

        // sorting lexicographicaly
        Fl.sort(F_comp);
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
                gamma_poles.append(F_pow+w_sum);
                cout<<"end achieved"<<endl;
              }
            else
              {
                string str = "w"+boost::lexical_cast<string>(displacement + w_index);
                //symbol w(str);
                w_lst.append(get_symbol(str));
                coeff*=tgamma(-get_symbol(str));
                x_power = get_symbol(str);
                cout<<"x_power "<<x_power<<endl;
                gamma_poles.append(-x_power); //!!!! review
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
        cout<<"Gammas after MB: "<<endl<<gamma_poles<<endl;
        cout<<"X powers list:"<<"  "<<x_power_map<<endl;
        //!!!!!!!!!!!!!!!!!!
        //        assert(false);
        // applying X integration
        ex gamma_den = 0; // gamma in denominator
        for(exmap::const_iterator mi = x_power_map.begin();mi!=x_power_map.end();++mi)
          {
            cout<<(*mi).first<<" "<<(*mi).second<<endl;
            gamma_poles.append((*mi).second+1);
            coeff*=tgamma((*mi).second+1);
            gamma_den+=((*mi).second+1);
          }
        cout<<"GAMMA_DEN: "<<gamma_den<<endl;
        bool gamma_den_has_w = false;
        for(lst::const_iterator wit = w_lst.begin(); wit != w_lst.end(); ++wit)
          if(gamma_den.has(*wit))gamma_den_has_w = true;
        if(gamma_den_has_w) gamma_poles.append(gamma_den);
        coeff/=tgamma(gamma_den);
        coeff/=(pow(2*Pi*I,w_lst.nops()));
        cout<<"New gamma list:"<<endl<<gamma_poles<<endl;
        full_int_expr = coeff;
        cout<<w_lst<<endl;
      }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"MBintegral)UFX)\":\n |___> ")+p.what());
      }

  }


MBintegral MBintegral::res(relational w_relation,ex pole,relational new_eps)
  {
    try
      {
        // remove W-residue from w-list
        exvector cut_w_vec(w_lst.nops()-1);
        std::remove_copy(w_lst.begin(),w_lst.end(),cut_w_vec.begin(),w_relation.lhs());
        cout<<"removed "<<w_relation.lhs()<<"  in list "<<lst(cut_w_vec.begin(),cut_w_vec.end())<<endl;

        lst new_gamma_pole_list;
        for(lst::const_iterator it = gamma_poles.begin();it!=gamma_poles.end();++it)
          if(*it != pole && has_w(it->subs(w_relation),w_lst).nops() > 0) new_gamma_pole_list.append(it->subs(w_relation));
        cout<<"removed pole list  "<<pole.subs(w_relation)<<"  in list "<<new_gamma_pole_list<<endl;
        // !!! IMPORTANT!!! no 2*pi*i multiplication and no sign multiplication 

        if(full_int_expr.denom().has(tgamma(pole))) assert(false);
        cout<<"  TAKING RES ON:  "<<full_int_expr<<endl;

        cout<<"RES LORAN: "<< full_int_expr.series(w_relation,0).coeff(w_relation.lhs(),-1)<<endl;
        /*
          REsidue by Loran serties !!!!!
         */
        ex res_loran = full_int_expr.series(w_relation,0).coeff(w_relation.lhs(),-1);

        ex new_no_gamma_part = (full_int_expr.subs(tgamma(pole)==pow(-1,-pole.subs(w_relation))/factorial(-pole.subs(w_relation)))).subs(w_relation);
        // new_no_gamma_part  = pow(-1,pole.subs(w_relation))/factorial(pole.subs(w_relation))*full_int_expr.subs(w_relation);
        cout<< new_no_gamma_part<<endl;
        exmap new_w_current(w_current);
        cout<<" Not modif:  "<<new_w_current<<endl;

        new_w_current.erase(w_relation.lhs());
        cout<<" modif:  "<<new_w_current<<endl;
        cout<<"CHECK:  "<<new_gamma_pole_list<<endl;
        //        MBintegral resINT(lst(cut_w_vec.begin(),cut_w_vec.end()),new_gamma_pole_list,new_no_gamma_part,new_w_current,new_eps);

        MBintegral resINT(lst(cut_w_vec.begin(),cut_w_vec.end()),new_gamma_pole_list,res_loran,new_w_current,new_eps);
        return resINT;
      }catch(std::exception &p)
      {
        throw std::logic_error(std::string("In function \"MBintegral.res\":\n |___> ")+p.what());
      }
  }

lst MBintegral::gamma_args_with_eps()
  {
    lst lst_with_eps;
    for(lst::const_iterator it = gamma_poles.begin(); it != gamma_poles.end(); ++it)
      if(it->has(get_symbol("eps"))) lst_with_eps.append(*it);
    lst_with_eps.unique();
    return lst_with_eps;
  }


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
    gamma_poles = set2lst(poles_set);
        cout<<"Update poles "<<poles_set<<endl;
  }

  exmap MBintegral::new_point()
  {
    lst var_list(w_lst);
    var_list.append(get_symbol("eps"));
    eps_w_current=start_point_diff_w(get_pole_lst(),var_list);
    cout<<"Int: "<<  interior_point(get_pole_lst(),eps_w_current)<<endl;
    BOOST_ASSERT_MSG(interior_point(get_pole_lst(),eps_w_current),"Not a convex polyhedron interior point");
    //eps_w_current = start_point_convex(get_pole_lst(),var_list);
    //eps_w_current=findinstance(get_pole_lst(),var_list);
    for(lst::const_iterator it = w_lst.begin();it!=w_lst.end();++it)
      w_current[*it] = it->subs(eps_w_current);
    eps_current = (get_symbol("eps")==eps_w_current[get_symbol("eps")]);
    return eps_w_current;
  }





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

          cout<<endl<<"Epsilon continue from eps_i = "<<eps_i<<endl<<endl;

          // Iterate over gamma arguments with eps dependence only!!!!!!!
          lst with_eps_lst(it->gamma_args_with_eps());
          for(lst::const_iterator pit  = with_eps_lst.begin(); pit != with_eps_lst.end(); ++pit)
            {
              //   cout<<"F(eps_i) "<<pit->subs(it->get_w()).subs(it->get_eps())<<"F(eps=0) "<<pit->subs(it->get_w()).subs(get_symbol("eps")==eps0)<<"   min  "<<std::min(pit->subs(it->get_w()).subs(it->get_eps()),pit->subs(it->get_w()).subs(get_symbol("eps")==eps0))<<endl;
             
              ex F_eps0 = pit->subs( it->get_w() ).subs( get_symbol("eps") == eps0) ;
              ex F_epsi =  pit->subs( it->get_w() ).subs( it->get_eps() ) ;

              if(F_eps0==F_epsi) 
                cout<<"Terminating eps=0 achieved  "<<std::min(F_eps0,F_epsi)<<endl;

              for(int n = 0;n>std::min(F_eps0,F_epsi);n--)
                {
                  // cout<<pit->subs(it->get_w()) <<endl;
                  // test on epsilon existance
                  if(pit->subs(it->get_w()).has(get_symbol("eps")))
                    {
                      ex eps_prime = lsolve(pit->subs(it->get_w()) ==n,get_symbol("eps") );
                      //       cout<<"solve"<<endl;
                      // cout<<"F= "<<*pit<<endl;
                      // cout<<"eps_i: "<<lsolve(pit->subs(it->get_w()) ==n,get_symbol("eps") )<<endl;
                      // cout<<" Poles of Gamma on eps_i: "<<it->get_pole_lst().subs(it->get_w()).subs(get_symbol("eps")==eps_prime)<<endl;
                      lst w_in_F  = has_w(*pit,it->get_w_lst());
                      if(w_in_F.nops()>0)
                        {
                          // cout<<lsolve(*pit==n,w_in_F.op(0))<<endl;
                          //   cout<<"sign(z) = "<<csgn(pit->coeff(w_in_F.op(0)))<<"     sign(F_i-F_0) = "<<csgn(F_epsi-F_eps0)<<endl;
                      
                          //MBintegral newi(it->get_w_lst(),it->get_pole_lst(),it->get_expr(),it->get_w(),it->get_eps());
                          //      cout<<"NEW INT: "<< newi.get_expr()<<endl;
                          // cout<<"debug fepsi"
                          //   <<"1: "<<lsolve(*pit==n,w_in_F.op(0))<<endl
                          //   <<"2: "
                          //   <<endl;
                          cout<<"BEFORE RESIDUE!: "<<it->get_w_lst()<<endl
                              <<it->get_expr()<<endl;
                          MBintegral res_int = it->res(w_in_F.op(w_in_F.nops()-1)==lsolve(*pit==n,w_in_F.op(w_in_F.nops()-1)),*pit,get_symbol("eps")==eps_prime);
                          res_int.set_level(1+it->get_level());
                          cout<<"after RESIDUE!: "<<res_int.get_w_lst()<<endl
                              <<res_int.get_expr()<<endl;
                          res_int.update_poles_from_ex();

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
                }

            }
        }
      O = R;
    }
  cout<<"Continue get "<<C.size()<<" integrals"<<endl;
  cout<< endl<<" Next step?  [Y/n]: ";
  
  char in_ch;
  std::cin>>in_ch;
  if(in_ch=='n')  exit(0);//assert(false);

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

ex expand_and_integrate(MBintegral& int_in, lst num_subs, int expansion_order) // up to O(eps^1) 
{
  try
    {
      ex out_ex;

      //            cout<<(it->get_pole_lst().subs(it->get_w())).subs(get_symbol("eps")==0)<<endl;
      //         cout<<endl <<"series:  " << it->get_gamma_expr().series(get_symbol("eps")==0,expansion_order)<<endl<<endl;
      lst w_lst = int_in.get_w_lst();
      if(int_in.get_w_lst().nops()>0)// expanding and integrating
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
          for(int i = out_ex.ldegree( get_symbol("eps") ); i <= out_ex.degree( get_symbol("eps") ); i++)
            {
              cout<<"Ord( "<<i<<" ) coeff : "<< out_ex.coeff(get_symbol("eps"),i)<<endl;
              ex int_expr =  out_ex.coeff(get_symbol("eps"),i);
              RoMB::FUNCP_CUBA2 fp_real,fp_imag;
              RoMB::compile_ex_real(lst(int_expr),int_in.get_w_lst(), fp_real);
              RoMB::compile_ex_imag(lst(int_expr),int_in.get_w_lst(), fp_imag);

              // ----------------------------------- Vegas integration-------------------------
               int  NDIM  = int_in.get_w_lst().nops();
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
      throw std::logic_error(std::string("In function \"RoMB\":\n |___> ")+p.what());
    }
}
