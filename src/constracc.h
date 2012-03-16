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

    exmap epsAndWsCurrent_;

    exmap WsCurrent_;


// PPL part

    NNC_Polyhedron ph_;

    std::vector<Variable> varVector;

// Converts ex to linear expression

    Linear_Expression ExToLe(ex);

    exmap chebyshevSphere(MBintegral::w_lst_type, MBintegral::p_lst_type);
    
public:
    

    ConstrAcc();
    
    ConstrAcc(const MBintegral::p_lst_type&, const MBintegral::w_lst_type&);

    ConstrAcc(const MBintegral&);

    const exmap& GetPoint() const;

    const exmap& GetWs() const;

    void PrintPoint();
        
    void PrintWs();

    bool Restrict(NearestPoleParams);

    

    //  constr_acc(MBintegral::w_lst_type,MBintegral::w_lst_type);
    
    bool add_single(const ex&);
    
    bool test_single(const ex&);
    
    size_t add_lst(lst&);
    
    bool test_lst(lst& cl);
    
};



// function test exmap is suitable solution of system of inequalities

bool interior_point(MBintegral::p_lst_type , exmap );

#endif // __CONSTRACC_H__
