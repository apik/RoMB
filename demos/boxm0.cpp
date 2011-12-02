#include "romb.h"
int main()
{
try
  {

     symbol k("k"),q("q"),p("p"),p1("p1"),p2("p2"),p3("p3"),ms("ms"),l("l"),s("s"),m1s("m1s"),m2s("m2s"),m3s("m3s");
    symbol l1("l1"),l2("l2"),l3("l3"),l4("l4"),t("t"),p4("p4"),p5("p5"),tp("tp"),v1("v1"),v2("v2"),l5("l5");
    symbol k1("k1"),k2("k2"),k3("k3"),k4("k4"),k5("k5"),ms1("ms1"),ms2("ms2"),ms3("ms3"),ms4("ms4");

                                      // M=0 with factor tgamma(1-eps)^2/tgamma(1-2eps)
                                      /*      PJfry
                                      *	1/eps^-2 :(-0.111111,0)
                                      *	1/eps^-1 :(0.0856421,0)
                                      *	1/eps^0 :(0.0513422,0)  -3.28987
                                      
                                      */
                                      
    
          RoMB_loop_by_loop box1loopm0(lst(k),lst(-pow(k,2),-pow(k+p1,2),-pow(k+p1+p2,2),-pow(k+p1+p2+p4,2)),
                                  lst(pow(p1,2)==0,pow(p2,2)==0,pow(p4,2)==0,
                                      p1*p2==s/2,//
                                      
                                      p1*p4==-s/2-t/2,//
                                      
                                      p2*p4==t/2 //
                                      ),
                                       lst(1,1,1,1),true);
    box1loopm0.integrate_map(lst(s==-5,t==-1),3);
   

  }
  catch(std::exception &p)
    {
      std::cerr<<"******************************************************************"<<endl;
      std::cerr<<"   >>>ERROR:  "<<p.what()<<endl;
      std::cerr<<"******************************************************************"<<endl;
      return 1;
    }
  return 0;
}
