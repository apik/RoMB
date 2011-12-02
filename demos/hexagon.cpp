#include "romb.h"
int main()
{
  try
    {
      symbol k("k"),q("q"),p("p"),p1("p1"),p2("p2"),p3("p3"),ms("ms"),l("l"),s("s"),m1s("m1s"),m2s("m2s"),m3s("m3s");
      symbol l1("l1"),l2("l2"),l3("l3"),l4("l4"),t("t"),p4("p4"),p5("p5"),p6("p6"),tp("tp"),v1("v1"),v2("v2"),l5("l5");
      symbol k1("k1"),k2("k2"),k3("k3"),k4("k4"),k5("k5"),ms1("ms1"),ms2("ms2"),ms3("ms3"),ms4("ms4");
      symbol s12("s12"),s23("s23"),s34("s34"),s45("s45"),s51("s51"),s13("s13"),s15("s15"),s56("s56"),s16("s16"),s123("s123"),s234("s234"),s345("s345");
      lst inv_l;
      inv_l.append(p1*p1 == 0);
      inv_l.append( p2*p2 == 0);inv_l.append( p3*p3  ==  0);inv_l.append( p4*p4  ==  0);inv_l.append( p5*p5  ==  0);inv_l.append( p6*p6  ==  0);
      inv_l.append(p1* p2  ==  s12/2);inv_l.append( p2* p3  ==  s23/2);inv_l.append( p3* p4  ==  s34/2);inv_l.append( p4* p5  ==  s45/2);
      inv_l.append(p5* p6  ==  s56/2);inv_l.append( p1* p6  ==  s16/2);inv_l.append( p1* p3  ==  (-s12 + s123 - s23)/2);
      inv_l.append(p2* p4  ==  (-s23 + s234 - s34)/2);
      inv_l.append( p3* p5  ==  (-s34 + s345 - s45)/2);
      inv_l.append(p1* p4  ==  (-s123 + s23 - s234 + s56)/2);
      inv_l.append(p1* p5  ==  (-s16 + s234 - s56)/2);
      inv_l.append( p2* p5  ==  (s16 - s234 + s34 - s345)/2);
      inv_l.append( p2* p6  ==  (-s12 - s16 + s345)/2);
      inv_l.append( p3* p6  ==  (s12 - s123 - s345 + s45)/2);
      inv_l.append( p4* p6  ==  (s123 - s45 - s56)/2);
      
      
      RoMB_loop_by_loop hexag(lst(k1),
                              lst(-pow(p1 + k1,2),-pow(p1 + p2 + k1,2),
                                  -pow(p1 + p2 + p3 + k1,2),
                                  -pow(p1 + p2 + p3 + p4 + k1,2),
                                  -pow(p1+p2+p3+p4+p5+k1,2),-pow(k1,2)),
                              inv_l,
                              lst(1,1,1,1,1,1),true);
      hexag.integrate_map(lst(s12 == -1, s23 == -2, s34 == -3, s45 == -4, s56 == -5, s16 == -6, s123 == -7, s234 == -8, s345 == -9));
      
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
