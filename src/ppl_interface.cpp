#include "ppl_interface.h"

//#undef PPL_NO_AUTOMATIC_INITIALIZATION 
#include <ppl.hh>
using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;


//include "ppl_c.h"

exmap chebyshevSphere(MBintegral::w_lst_type wIn, MBintegral::p_lst_type pIn) 
{

try
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
//                mpz_class 


                GMP_Integer gmpInt(ex_to<numeric>(linCoeff).to_long());
                l = sub_mul_assign(l, gmpInt, varVector[dimN]);
                dimN++;
            }
            else throw std::logic_error(std::string( "Not an integer coefficient in Linear_Expression" ));
        }
// Add B_i        
        if(bI.info(info_flags::integer))
        {
            GMP_Integer gmpInt(ex_to<numeric>(bI).to_long());
            l -= gmpInt;
        }
        else throw std::logic_error(std::string( "Not an integer B_i in Linear_Expression" ));

// Add radius
        
        GMP_Integer aNormSqrt;
        if(aNormSquared.info(info_flags::integer))
        {

            GMP_Integer gmpInt(ex_to<numeric>(aNormSquared).to_long());

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
    cout << " System created " << cs << endl;
    NNC_Polyhedron ph(cs);
    bool maxVal;
    Coefficient supN,supD;
    Generator g  = closure_point();
    if( ph.maximize(r,supN,supD,maxVal,g) )
    {
        std::cout<< maxVal<<" then " << supN<<"/"<<supD<< "Generator: " << g << std::endl;
       
    }
    else throw std::logic_error(std::string( "No solution for Chebyshev Sphere" ));

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
        intPtest[w] = numeric(g.coefficient(varVector[cncn]).get_si(),g.divisor().get_si()) + rEx/(*primeToDivide);
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


