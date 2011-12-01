#include "romb.h"
int main()
{
  try
    {
  
      symbol k("k"),q("q"),p("p"),p1("p1"),p2("p2"),p3("p3"),ms("ms"),l("l"),s("s"),m1s("m1s"),m2s("m2s"),m3s("m3s");
      symbol l1("l1"),l2("l2"),l3("l3"),l4("l4"),t("t"),p4("p4"),p5("p5"),tp("tp"),v1("v1"),v2("v2"),l5("l5");
      symbol k1("k1"),k2("k2"),k3("k3"),k4("k4"),k5("k5"),ms1("ms1"),ms2("ms2"),ms3("ms3"),ms4("ms4");
      symbol s12("s12"),s23("s23"),s34("s34"),s45("s45"),s51("s51"),s13("s13"),s15("s15");

      lst inv_l =  lst(
                       p1*p1 == 0, 
                       p2*p2 == 0, 
                       p3*p3 == 0, 
                       p4*p4 == 0, 
                       //p5*p5 == 0, 

                       p1*p2 == s12/2,
                       p1*p3 ==  (-s12-s23+s45)/2,
                       p1*p4 == (-s15+s23-s45)/2,
                       p1*p5 == s15/2,
                       p2* p3 ==s23/2,
                       p2* p4 ==(s15-s23-s34)/2,
                       p2* p5 ==(-s12-s15+s34)/2,
                       p3*p4 == s34/2,
                       p3*p5 == (s12-s34-s45)/2,
                       p4*p5 == s45/2);
      
      RoMB_loop_by_loop pent(lst(k1),lst(-pow(p1 + k1,2),-pow(p1 + p2 + k1,2),
                                         -pow(p1 + p2 + p3 + k1,2),-pow(p1 + p2 + p3 + p4 + k1,2),
                                         -pow(k1,2)),
                             inv_l,
                             lst(1,1,1,1,1),false);
      pent.integrate_map(lst(s12==-2,s23==-3,s34==-4,s45==-5,s15==-6,ms==0));
       
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
