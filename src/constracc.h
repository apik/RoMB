#ifndef __CONSTRACC_H__
#define __CONSTRACC_H__
#include <ginac/ginac.h>
#include "mbintegral.h"
using namespace GiNaC;
class ConstrAcc
{

    lst constraints;
    
    lst w_lst;
    
public:
    

    
    ConstrAcc(MBintegral::p_lst_type,MBintegral::w_lst_type);

    ConstrAcc(const MBintegral&);
    
    //  constr_acc(MBintegral::w_lst_type,MBintegral::w_lst_type);
    
    bool add_single(const ex&);
    
    bool test_single(const ex&);
    
    size_t add_lst(lst&);
    
    bool test_lst(lst& cl);
    
};


exmap chebyshevSphere(MBintegral::w_lst_type, MBintegral::p_lst_type);


// function test exmap is suitable solution of system of inequalities

bool interior_point(MBintegral::p_lst_type , exmap );

#endif // __CONSTRACC_H__
