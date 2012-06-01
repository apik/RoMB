
#include "constracc.h"
#include "utils.h"

#include <boost/numeric/interval.hpp> 
//#include <boost/numeric/interval/compare/lexicographic.hpp>
#include <boost/numeric/interval/compare/explicit.hpp>



/*constr_acc::constr_acc(lst constr_lst,lst w_lst_in) :
  constraints_(constr_lst),ws_(w_lst_in)
{
}
*/

ConstrAcc::ConstrAcc()
{
}

 // Constructor  from MBIntegral

ConstrAcc::ConstrAcc(const MBintegral& mbIn)
{
 
/*   BOOST_FOREACH(ex ce,mbIn.get_poles())
    {
        constraints_.append(ce);
    }
    BOOST_FOREACH(ex we, mbIn.get_w_lst())
    {
        ws_.append(we);
    }

*/

    *this = ConstrAcc(mbIn.get_poles(), mbIn.get_w_lst());
}



// Constructor  from explicit poles and arg lists

ConstrAcc::ConstrAcc(const MBintegral::p_lst_type& constr_lst,
                     const MBintegral::w_lst_type& w_lst_in)
{
try
{

  BOOST_FOREACH(ex we,w_lst_in)
    {
        ws_.append(we);
  
    }



/* 
 *
 *   DANGER!!! eps by hand!!!
 *   No eps in w_lst_in
 *
 *
 */ 
 


        wsAndEps_ = w_lst_in;
    wsAndEps_.push_back(get_symbol("eps"));


// push only w's and eps

    size_t varN = 0;
    BOOST_FOREACH(ex w, wsAndEps_)
    {
        Variable vNew(varN);
        varVector.push_back(vNew);
        varN++;
    }
    
    
// Sphere radius
    
    R_ = new Variable(varN);

    Constraint_System cs,csWithR;
    
    BOOST_FOREACH(ex p, constr_lst)
    {
        
        Linear_Expression l = ExToLe(p);
        cout <<"================================================" << endl;
        Linear_Expression lWithR = ExToLeMinusA(p);

        cs.insert(l <= 0);
        csWithR.insert(lWithR <= 0);
        

    }
    cout << " System created " << cs << endl;
    ph_ = NNC_Polyhedron(cs);
    phR_ = NNC_Polyhedron(csWithR);




/*
  ==============================================================================
  =                                                                            =
  =                                                                            =
  =                                                                            =
  ==============================================================================
*/

//    epsAndWsCurrent_ = chebyshevSphere(wsAndEps, constr_lst);

    epsAndWsCurrent_ = chebyshevSphere();

    

    cout << "test" << epsAndWsCurrent_ << endl;
    PrintPoint();

    
    BOOST_FOREACH(ex ce,constr_lst)
    {
        constraints_.append(ce);
    }
    BOOST_FOREACH(ex we,w_lst_in)
    {
//        ws_.append(we);

        if(we != get_symbol("eps"))
            WsCurrent_[we] = epsAndWsCurrent_[we]; 
    }


    }catch(std::exception &p)
    {
        throw std::logic_error(std::string("In function \"ConstrAcc\":\n |___> ")+p.what());
    }
}


const exmap& ConstrAcc::GetPoint() const
{
    return epsAndWsCurrent_;
}

const exmap& ConstrAcc::GetWs() const
{
    return WsCurrent_;
}

void ConstrAcc::PrintPoint() 
{
    cout << "Print in constr\n\t" <<  epsAndWsCurrent_ << endl;
}

void ConstrAcc::PrintWs() 
{
    cout << "Print in constr\n\t" <<  WsCurrent_ << endl;
}


/*
bool ConstrAcc::add_single(const ex& constr)
{
  lst add_lst(constraints_);
  add_lst.append(constr);
  if(zero_volume(add_lst,ws_))return false;
  else
    {
      constraints_.append(constr);
      return true;
    }
}
bool ConstrAcc::test_single(const ex& constr)
{
  lst add_lst(constraints_);
  add_lst.append(constr);
  return !zero_volume(add_lst,ws_);
}
bool ConstrAcc::test_lst(lst& cl)
{
  lst add_lst(constraints_);
  for(lst::const_iterator lit = cl.begin(); lit != cl.end(); ++lit)
    add_lst.append(*lit);
  return !zero_volume(add_lst,ws_);
}
*/



exmap ConstrAcc::chebyshevSphere(MBintegral::w_lst_type wIn, MBintegral::p_lst_type pIn) 
{

    try
    {

        size_t varN = 0;
        BOOST_FOREACH(ex w, wIn)
        {
            Variable vNew(varN);
            varVector.push_back(vNew);
            varN++;
        }

    
// Sphere radius

        R_ = new Variable(varN);

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
//                mpz_class 


                    GMP_Integer gmpInt(ex_to<GiNaC::numeric>(linCoeff).to_long());
                    l = sub_mul_assign(l, gmpInt, varVector[dimN]);
                    dimN++;
                }
                else throw std::logic_error(std::string("Not an integer coefficient in Linear_Expression" ));
            }
// Add B_i        
            if(bI.info(info_flags::integer))
            {
                GMP_Integer gmpInt(ex_to<GiNaC::numeric>(bI).to_long());
                l -= gmpInt;
            }
            else throw std::logic_error(std::string( "Not an integer B_i in Linear_Expression" ));

// Add radius
        
            GMP_Integer aNormSqrt;
            if(aNormSquared.info(info_flags::integer))
            {

                GMP_Integer gmpInt(ex_to<GiNaC::numeric>(aNormSquared).to_long());

// Assign integer value of sqare root: 
//
//      int(sqrt) < sqrt
//
                sqrt_assign(aNormSqrt, gmpInt);
                aNormSqrt += 1;
                cout << "SQRT: " << aNormSqrt << endl;
            }
            else throw std::logic_error(std::string( "Not an integer B_i in Linear_Expression" ));

            l = add_mul_assign(l, aNormSqrt, *R_);
        
            cout << l << endl;
            cs.insert(l <= 0);
        }
        cout << " System created " << cs << endl;
        ph_ = NNC_Polyhedron(cs);

        bool maxVal;
        Coefficient supN,supD;
        Generator g  = closure_point();
        if( ph_.maximize(*R_,supN,supD,maxVal,g) )
        {
            std::cout<< maxVal<<" then " << supN<<"/"
                     <<supD<< " Generator: " << g << std::endl;
       
        }
        else throw std::logic_error(std::string( "No solution for Chebyshev Sphere" ));

        ex rEx = GiNaC::numeric(supN.get_si(),supD.get_si());

        Generator_System gs;
        gs.insert(g);
        NNC_Polyhedron ph2(gs);
        if(ph_.strictly_contains(ph2)) cout<< " Contains" <<endl;
    
    
// 168 first prime numbers
        const unsigned int primeNumbers[168] = {
// no first 28
//        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
            109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 
            181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 
            257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 
            337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 
            419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 
            491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 
            587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 
            659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 
            751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 
            839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 
            937, 941, 947, 953, 967, 971, 977, 983, 991, 997 
        };
      
        const unsigned int* primeNumbersEnd = primeNumbers + 167 - 28;


        cout<< " dimension " << g.space_dimension() << endl;
    
// no r variable
        long int largestNumerator = 0;
        for(dimension_type i = 0; i < g.space_dimension() - 2; i++)
            if(largestNumerator < abs(g.coefficient(varVector[i]).get_si()))
                largestNumerator = abs(g.coefficient(varVector[i]).get_si());
    
        cout << largestNumerator << endl;        
        const unsigned int* primeToDivide = upper_bound (primeNumbers, primeNumbersEnd, largestNumerator);

        if (primeToDivide == primeNumbersEnd)
            throw std::logic_error(string("Unavailable prime number to divide"));

        exmap intPtest;
        size_t cncn = 0;
        BOOST_FOREACH(ex w, wIn)
        {
            intPtest[w] = GiNaC::numeric(g.coefficient(varVector[cncn]).get_si(),
                                  g.divisor().get_si()) + rEx/(*primeToDivide);
            primeToDivide++;
            cncn++;
        }
        
//       if (binary_search (primeNumbers, primeNumbersEnd, g.coefficient(varVector[i]).get_si()))

//        cout << "found!\n"; else cout << "not found.\n";
    
        cout<< intPtest <<endl;
    
//     BOOST_ASSERT_MSG(interior_point(pIn,intPtest),"Not a convex polyhedron interior point");


        Parma_Polyhedra_Library::restore_pre_PPL_rounding();
        return intPtest;
    }catch(std::exception &p)
    {
        throw std::logic_error(std::string("In function \"chebyshevSphere\":\n |___> ")+p.what());
    }
}


struct classcomp {
    bool operator() (const boost::numeric::interval<mpq_class>& lhs, const boost::numeric::interval<mpq_class>& rhs) const
        {return  boost::numeric::interval_lib::cerlt(lhs,rhs);}
};


exmap ConstrAcc::DiffContours(const Generator& generator, const mpq_class& r)
{
    using namespace boost::numeric;

        
//    cout << "In diff: " << generator  << " R " << r <<endl;

    typedef interval<mpq_class> IntervalQ;
    typedef std::multimap<interval<mpq_class>, ex,classcomp> IntervalMap;

    IntervalMap contourIntervals;
 
    size_t cncn = 0;
    BOOST_FOREACH(ex w, wsAndEps_)
//        for(lst::const_iterator wii = ws_.begin(); wii != ws_.end(); ++wii)
    {
//            intPtest[w] = numeric(g.coefficient(varVector[cncn]).get_si(),
//                                  g.divisor().get_si()) + rEx/(*primeToDivide);
        
        //intPtest[w] = GiNaC::numeric(g.coefficient(varVector[cncn]).get_str().c_str())/
        //    GiNaC::numeric(g.divisor().get_str().c_str()) + rEx/(*primeToDivide);

// rational coordinate
        mpq_class w_mpq(generator.coefficient(varVector[cncn]),generator.divisor());
        

        IntervalMap::iterator it = contourIntervals.insert(IntervalMap::value_type(IntervalQ(w_mpq - r, w_mpq + r), w));

//        cout << "W: " << width(it->first) << " [" << it->first.lower() << ", " << it->first.upper() << "]" <<endl;

        //       cout << " ORDERED:" << endl;
        
        
//        primeToDivide++;
        cncn++;
    }

/*for (IntervalMap::iterator im = contourIntervals.begin(); im != contourIntervals.end(); ++im)
        {
//            cout << "W: " << width(im->first) << " [" << im->first.lower() << ", " << im->first.upper() << "]" <<endl;

        }
*/

// New map for non intersecting intervals
    
//    IntervalMap contourIntervalsDiff;
    
    typedef IntervalMap::iterator mapIter;
    

    exmap outContours;

// Upper bound of last interval

    mpq_class lastMax;

    mapIter m_it, s_it;
                    
    for (m_it = contourIntervals.begin();  m_it != contourIntervals.end();  m_it = s_it) 
    {
        IntervalMap::key_type theKey = m_it->first;
                        
        //       cout << endl;
        //  cout << "  key = '" << theKey << "'" << endl;

        // working with key here:

        size_t nWithKey = contourIntervals.count(theKey);                      
                      
        pair<mapIter, mapIter> keyRange = contourIntervals.equal_range(theKey);
                      
        // Iterate over all map elements with key == theKey
        
        size_t nInRange = 1;
                      
        for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it) 
        {
            // absolute lower bound
            if(s_it == contourIntervals.begin())
                lastMax = s_it->first.lower();

            
            mpq_class lb = std::max(lastMax, s_it->first.lower());

            lastMax = lb +  width(s_it->first) * mpq_class(nInRange,nWithKey);

//            cout << "lb: " << lb << " ub: " << lastMax<< " w: " << s_it->second << endl;
            
            
            mpq_class median_q((lb + lastMax)/2);
            
            //           cout << "Med " << median_q <<endl;
            
            outContours[s_it->second] = GiNaC::numeric(median_q.get_num().get_str().c_str())/
                GiNaC::numeric(median_q.get_den().get_str().c_str());
      
            nInRange++;
            
        }
    }

//    cout << outContours <<endl;
    return outContours;   
}



exmap ConstrAcc::chebyshevSphere()
{

    try
    {

        bool maxVal;
        Coefficient supN,supD;
        Generator g  = closure_point();
        if( phR_.maximize(*R_,supN,supD,maxVal,g) )
        {
            std::cout<< maxVal<<" then " << supN<<"/"
                     <<supD<< " Generator: " << g << std::endl;
       
        }
        else throw std::logic_error(std::string( "No solution for Chebyshev Sphere" ));

        ex rEx = GiNaC::numeric(supN.get_si(),supD.get_si());

        mpq_class r_mpq(supN,supD);

        exmap ommap =  DiffContours(g, r_mpq);

        Parma_Polyhedra_Library::restore_pre_PPL_rounding();
        return  ommap;

        Generator_System gs;
        gs.insert(g);
        NNC_Polyhedron ph2(gs);
        if(phR_.strictly_contains(ph2)) cout<< " Contains" <<endl;
    
    
// 168 first prime numbers
        const unsigned int primeNumbers[168] = {
// no first 28
//        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
            109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 
            181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 
            257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 
            337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 
            419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 
            491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 
            587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 
            659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 
            751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 
            839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 
            937, 941, 947, 953, 967, 971, 977, 983, 991, 997 
        };
      













        const unsigned int* primeNumbersEnd = primeNumbers + 167 - 28;


        cout<< " dimension " << g.space_dimension() << endl;
    
// no r variable
        long int largestNumerator = 0;
        for(dimension_type i = 0; i < g.space_dimension() - 2; i++)
            if(largestNumerator < abs(g.coefficient(varVector[i]).get_si()))
                largestNumerator = abs(g.coefficient(varVector[i]).get_si());
    
        cout << largestNumerator << endl;        
        const unsigned int* primeToDivide = upper_bound (primeNumbers, primeNumbersEnd, largestNumerator);

        if (primeToDivide == primeNumbersEnd)
            throw std::logic_error(string("Unavailable prime number to divide"));

        exmap intPtest;
        size_t cncn = 0;
        BOOST_FOREACH(ex w, wsAndEps_)
//        for(lst::const_iterator wii = ws_.begin(); wii != ws_.end(); ++wii)
        {
//            intPtest[w] = numeric(g.coefficient(varVector[cncn]).get_si(),
//                                  g.divisor().get_si()) + rEx/(*primeToDivide);

            intPtest[w] = GiNaC::numeric(g.coefficient(varVector[cncn]).get_str().c_str())/
                GiNaC::numeric(g.divisor().get_str().c_str()) + rEx/(*primeToDivide);


            primeToDivide++;
            cncn++;
        }
        
//       if (binary_search (primeNumbers, primeNumbersEnd, g.coefficient(varVector[i]).get_si()))

//        cout << "found!\n"; else cout << "not found.\n";
    
        cout<< intPtest <<endl;
    
//     BOOST_ASSERT_MSG(interior_point(pIn,intPtest),"Not a convex polyhedron interior point");


     
        return intPtest;
    }catch(std::exception &p)
    {
        throw std::logic_error(std::string("In function \"chebyshevSphere\":\n |___> ")+p.what());
    }
}

GMP_Integer exToGmp(ex e)
{
    stringstream ss;
    if(e.info(info_flags::integer))
        ss << e;
    else throw std::invalid_argument(string("Not an integer number conversion to mpz_class"));
    return GMP_Integer(ss.str());
}


Linear_Expression ConstrAcc::ExToLe(const ex& pIn)
{
try
{
    ex p(pIn);

// Converting coefficients to integer

// Temporary value of LCM

    ex lcmBuf = 1;
    ex bBuf = p;
    // for(lst::const_iterator wi = ws_.begin(); wi != ws_.end(); ++wi)
   
    BOOST_FOREACH(ex w,wsAndEps_)
    {
        cout << "W " << w <<endl;
        if(p.coeff(w,1).info(info_flags::rational))
            lcmBuf = lcm(lcmBuf, denom(p.coeff(w,1)));
        else
            throw std::invalid_argument( string("Not rational coefficient in constraint") );
        bBuf = bBuf.coeff(w,0);
    }
    if(bBuf.info(info_flags::rational))
        lcmBuf = lcm(lcmBuf, denom(bBuf));
    else
        throw std::invalid_argument( string("Not rational coefficient in constriant") );
    
    cout <<"p " <<  p << endl << "LCM " << lcmBuf << endl;


// To conserve sign of inequality
// in p stored inequality with integer coefficients
    
    p *= abs(lcmBuf);

    Linear_Expression l;

    size_t dimN = 0;
    ex bI = p;
//    for(lst::const_iterator wi = ws_.begin(); wi != ws_.end(); ++wi)
    BOOST_FOREACH(ex w,wsAndEps_)
    {
        ex linCoeff = p.coeff(w,1);
        bI = bI.coeff(w,0);
        if(linCoeff.info(info_flags::integer))
        {
//                mpz_class 
            
//            GMP_Integer gmpInt(ex_to<GiNaC::numeric>(linCoeff).to_long());
            GMP_Integer gmpInt(exToGmp(linCoeff));
            l = sub_mul_assign(l, gmpInt, varVector[dimN]);
            dimN++;
        }
        else throw std::logic_error(std::string( "Not an integer coefficient in Linear_Expression" ));
    }
// Add B_i        
    if(bI.info(info_flags::integer))
    {
        //       GMP_Integer gmpInt(ex_to<GiNaC::numeric>(bI).to_long());
        GMP_Integer gmpInt(exToGmp(bI));
        l -= gmpInt;
    }
    return l;
    }catch(std::exception &p)
    {
        throw std::logic_error(std::string("In function \"ExToLe\":\n |___> ")+p.what());
    }
}


Linear_Expression ConstrAcc::ExToLeMinusA(const ex& pIn)
{
try
{
    ex p(pIn);

// Converting coefficients to integer

// Temporary value of LCM

    ex lcmBuf = 1;
    ex bBuf = p;
//    for(lst::const_iterator wi = ws_.begin(); wi != ws_.end(); ++wi)
    BOOST_FOREACH(ex w, wsAndEps_)
    {
        if(p.coeff(w,1).info(info_flags::rational))
            lcmBuf = lcm(lcmBuf, denom(p.coeff(w,1)));
        else
            throw std::invalid_argument( string("Not rational coefficient in constriant") );
        bBuf = bBuf.coeff(w,0);
    }
    if(bBuf.info(info_flags::rational))
        lcmBuf = lcm(lcmBuf, denom(bBuf));
    else
        throw std::invalid_argument( string("Not rational coefficient in constriant") );
    
    cout <<"p " <<  p << endl << "LCM " << lcmBuf << endl;


// To conserve sign of inequality
// in p stored inequality with integer coefficients
    
    p *= abs(lcmBuf);

    Linear_Expression l;

    size_t dimN = 0;
    ex bI = p;

    ex aNormSquared = 0;
    //for(lst::const_iterator wi = ws_.begin(); wi != ws_.end(); ++wi)
    BOOST_FOREACH(ex w, wsAndEps_)
    {
        ex linCoeff = p.coeff(w,1);
        bI = bI.coeff(w,0);
        aNormSquared += pow(linCoeff,2);
        if(linCoeff.info(info_flags::integer))
        {
//                mpz_class 
            
//            GMP_Integer gmpInt(ex_to<GiNaC::numeric>(linCoeff).to_long());
            GMP_Integer gmpInt(exToGmp(linCoeff));
            l = sub_mul_assign(l, gmpInt, varVector[dimN]);
            dimN++;
        }
        else throw std::logic_error(std::string( "Not an integer coefficient in Linear_Expression" ));
    }
// Add B_i        
    if(bI.info(info_flags::integer))
    {
        //       GMP_Integer gmpInt(ex_to<GiNaC::numeric>(bI).to_long());
        GMP_Integer gmpInt(exToGmp(bI));
        l -= gmpInt;
    }


// Add radius
        
            GMP_Integer aNormSqrt;
            if(aNormSquared.info(info_flags::integer))
            {

                GMP_Integer gmpInt(exToGmp(aNormSquared));

// Assign integer value of sqare root: 
//
//      int(sqrt) < sqrt
//
                sqrt_assign(aNormSqrt, gmpInt);
                aNormSqrt += 1;
                cout << "SQRT: " << aNormSqrt << endl;
            }
            else throw std::logic_error(std::string( "Not an integer B_i in Linear_Expression" ));

            l = add_mul_assign(l, aNormSqrt, *R_);

            return l;
    }catch(std::exception &p)
    {
        throw std::logic_error(std::string("In function \"ExToLeMinusA\":\n |___> ")+p.what());
    }
}


bool ConstrAcc::Restrict(const NearestPoleParams& nearestPoleParams)
{
    ex cs = nearestPoleParams.Arg;
    cs = cs.subs(get_symbol("eps") == nearestPoleParams.EpsilonValue);

    Linear_Expression le = ExToLe(cs);
    cout << "Lin ex " <<cs << " PPL " << le << endl;   
    
    Constraint constr = (le < 0);

// Test strict intersection of constraint and existing polyhedron
    if (ph_.relation_with(constr) == Poly_Con_Relation::strictly_intersects())
    {
        Linear_Expression leR = ExToLeMinusA(cs);
        cout << "Lin ex " << leR << endl;   
        
        Constraint constrR = (leR < 0);
        
        phR_.add_constraint(constrR);

        ph_.add_constraint(constr);

        epsAndWsCurrent_ = chebyshevSphere();
        cout << "   CONSTRAINT MAY DISCARDED!!!"<<endl;    
        return true;
    }
    else return false;
}


bool interior_point(MBintegral::p_lst_type ineq_lst, exmap subs_map)
{
  for(lst::const_iterator it = ineq_lst.begin(); it != ineq_lst.end(); ++it)
    {
        cout<<"int point: " << *it << " --> " <<it->subs(subs_map)<<endl;
       if((it->subs(subs_map) <0)) return false;
    }
  return true;
}
 
