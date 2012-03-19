#ifndef __CONSTRACC_H__
#define __CONSTRACC_H__
#include <ginac/ginac.h>
#include "mbintegral.h"

#include <ppl.hh>
using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;

using namespace GiNaC;




class ConstrAcc
{

    lst constraints_;
    
    lst ws_;

    MBintegral::w_lst_type wsAndEps_;

    exmap epsAndWsCurrent_;

    exmap WsCurrent_;


// PPL part

    // Polyhedron with constrints
    // Ax>b

    NNC_Polyhedron ph_;


// Polyhedron for LP on Chebyshev sphere
// Ax+An*r>b
    
    NNC_Polyhedron phR_;

// All variables except radius(R)
    std::vector<Variable> varVector;

    Variable* R_;

// Converts ex to linear expression

    Linear_Expression ExToLe(const ex&);

// Converts ex to LE -sqrt(a^2)

    Linear_Expression ExToLeMinusA(const ex&);
    
    

    exmap chebyshevSphere(MBintegral::w_lst_type, MBintegral::p_lst_type);
    exmap chebyshevSphere();
    
public:
    

    ConstrAcc();
    
    ConstrAcc(const MBintegral::p_lst_type&, const MBintegral::w_lst_type&);

    ConstrAcc(const MBintegral&);

    const exmap& GetPoint() const;

    const exmap& GetWs() const;

    void PrintPoint();
        
    void PrintWs();

    bool Restrict(const NearestPoleParams&);

    

    //  constr_acc(MBintegral::w_lst_type,MBintegral::w_lst_type);
    
    bool add_single(const ex&);
    
    bool test_single(const ex&);
    
    size_t add_lst(lst&);
    
    bool test_lst(lst& cl);
    
};



// function test exmap is suitable solution of system of inequalities

bool interior_point(MBintegral::p_lst_type , exmap );

#endif // __CONSTRACC_H__
