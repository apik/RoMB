#include "romb.h"
int main()
{
try
  {

     symbol k("k"),q("q"),p("p"),p1("p1"),p2("p2"),p3("p3"),ms("ms"),l("l"),s("s"),m1s("m1s"),m2s("m2s"),m3s("m3s");
    symbol l1("l1"),l2("l2"),l3("l3"),l4("l4"),t("t"),p4("p4"),p5("p5"),tp("tp"),v1("v1"),v2("v2"),l5("l5");
    symbol k1("k1"),k2("k2"),k3("k3"),k4("k4"),k5("k5"),ms1("ms1"),ms2("ms2"),ms3("ms3"),ms4("ms4");
  // oneloop box
  //      UFXmap l45 = UF(lst(k),lst(pow(k,2),pow(k+p1,2),pow(k+p1+p2,2),pow(k+p1+p2+p3,2)),lst(pow(p1,2)==0,pow(p2,2)==0));
  // MBintegral root_int(l45,lst(1,1,1,1),1);

   //two loop box bubble
  // UFXmap l45 = UF(lst(k,l),lst(pow(k,2),pow(k+p1,2),pow(k+p1+p2,2),pow(l+p1+p2,2),pow(l+p1+p2+p3,2),pow(l,2),pow(k-l,2)),lst(pow(p1,2)==0,pow(p2,2)==0,pow(p3,2)==0));
  //MBintegral root_int(l45,lst(1,1,1,1,1,1,1),2);

  // B0
  //    UFXmap l45 = UF(lst(k),lst(ms-pow(k,2),ms-pow(-k,2)),lst(ms==1));
  // MBintegral root_int(l45,lst(1,1),1);

  // 2 loop sunrise
  //UFXmap l45 = UF(lst(k,q),lst(ms-pow(k,2),ms-pow(-q-k,2),ms-pow(q,2)),lst(ms==1));
  //MBintegral root_int(l45,lst(1,1,1),2);


  //RoMB_planar box2loop(lst(k,l),lst(pow(k,2),pow(k+p1,2),pow(k+p1+p2,2),pow(l+p1+p2,2),pow(l+p1+p2+p3,2),pow(l,2),pow(k-l,2)),lst(pow(p1,2)==0,pow(p2,2)==0,pow(p3,2)==0),lst(1,1,1,1,1,1,1),2);

  //  RoMB_planar  box1loop(lst(k),lst(pow(k,2),pow(k+p1,2)-ms,pow(k+p1+p2,2),pow(k+p1+p2+p3,2)),lst(pow(p1,2)==0,pow(p2,2)==0,pow(p3,2)==0,p1==0,p2==0,p3==0,ms==1),lst(1,1,1,1),1);

  //  RoMB_planar B0_1loop(lst(k),lst(pow(k,2)-ms,pow(p+k,2)-ms),lst(ms==0,pow(p,2)==1),lst(1,1),1);

  //  RoMB_planar C0_1loop(lst(k),lst(pow(k,2)-ms,pow(p1+k,2)-ms,pow(p1+p2+k,2)),lst(ms==1,pow(p1,2)==0,pow(p2,2)==0,p1*p2==50),lst(1,1,1),1);
//cout<<" new point "<<endl<<root_int.new_point()<<endl;
// cout<<" saved point "<<endl<<root_int.get_point()<<endl;
//  MBcontinue(root_int);
  //cout<<MB_lst(l45,lst(1,1,1,1),1).expand()<<endl;


  // RoMB_loop_by_loop box2loop(lst(k,l),lst(pow(k,2),pow(k+p1,2),pow(k+p1+p2,2),pow(l+p1+p2,2),pow(l+p1+p2+p3,2),pow(l,2),pow(k-l,2)),lst(pow(p1,2)==0,pow(p2,2)==0,pow(p3,2)==0),lst(1,1,1,1,1,1,1));
  //      RoMB_loop_by_loop t2(lst(k,l), lst(pow(k,2),pow(p+k,2),pow(p+k+l,2),pow(l,2),pow(k+l,2)),lst(pow(p,2)==1),lst(1,1,1,1,1));


  // works!!!
  //        RoMB_loop_by_loop sunset(lst(k,l), lst(pow(k,2)-1,pow(p-k-l,2)-4,pow(l,2)-5),lst(pow(p,2)==s),lst(1,1,1));
      
  //   RoMB_loop_by_loop sunset(lst(k,l), lst(pow(k,2)-m1s,pow(-k-l,2)-m2s,pow(l,2)-m3s),lst(pow(p,2)==s),lst(1,1,1));
  //                 sunset.integrate(lst(m1s==1,m2s==1,m3s==1,s==0),0);

    //     bubble sunset 2=loop
    //               RoMB_loop_by_loop sunset_bub(lst(k,l), lst(-pow(k,2)+ms,-pow(-k-l,2)+ms,-pow(l,2)+ms),lst(pow(p,2)==0),lst(1,1,1));
    //  sunset_bub.integrate_map(lst(ms==1,m2s==1,m3s==1,s==0),1);

    //     bubble sunset 3=loop
//#define TOPO 1
//#if TOPO==1    
/******************************************************************    
 *    FRESULT for parameters: {ms==1,m2s==1,m3s==1,s==0}
 *   
 *    FRESULT anl :           = 21.308685443306456902+23/3*eps^(-2)+2*eps^(-3)+35/2*eps^(-1)
 *    FRESULT num:           = 21.308685443306456902+(7.6666666666666666665)*eps^(-2)+(2.0)*eps^(-3)+(17.5)*eps^(-1)
 *    eps^-3 term: 2 +/- 0
 *    eps^-2 term: 23/3 +/- 0
 *    eps^-1 term: 35/2 +/- 0
 *    eps^0 term: 21.308685443306456902 +/- 0.01814768000077260732
 ***************************************************************/
//           RoMB_loop_by_loop sunset_bub(lst(p,k,l), lst(-pow(p,2)+ms,-pow(k,2)+ms,-pow(l,2)+ms,-pow(-p-k-l,2)+ms),lst(pow(l3,2)==s),lst(1,1,1,1));
//           sunset_bub.integrate_map(lst(ms==1,m2s==1,m3s==1,s==0),1);
//#elif TOPO==2
// RoMB_loop_by_loop sunset_bub_d(lst(l1,l2,l3), lst(-pow(l1,2)+ms,-pow(l2,2)+ms,-pow(l3,2)+ms,-pow(l1+l2,2)+ms,-pow(l1+l2+l3,2)+ms),lst(pow(p,2)==s),lst(1,1,1,1,1));
// sunset_bub_d.integrate_map(lst(ms==1,m2s==1,m3s==1,s==0),1);
//#elif TOPO==3
// RoMB_loop_by_loop sunset_bub_e(lst(l1,l2,l3), lst(-pow(l1,2)+ms,-pow(l2,2)+ms,-pow(l3,2)+ms,-pow(l1-l2,2)+ms,-pow(l2-l3,2)+ms,-pow(l3-l1,2)+ms),lst(pow(p,2)==s),lst(1,1,1,1,1,1));  
// sunset_bub_e.integrate_map(lst(ms==1,m2s==1,m3s==1,s==0),1);
//#endif
    //bubble 4-loop
    //            RoMB_loop_by_loop sunset_bub(lst(k,l1,l2,l3), lst(-pow(k,2)+ms,-pow(l2,2)+ms,-pow(l1,2)+ms,-pow(l3,2)+ms,-pow(k+l1+l2+l3,2)+ms),lst(pow(p,2)==s),lst(1,1,1,1,1));
    // sunset_bub.integrate(lst(ms==1,m2s==1,m3s==1,s==0),0);
    
//bubble 5-loop
/*
FRESULT for parameters: {ms==1,m2s==1,m3s==1,s==0}

 FRESULT anl :           = 274.5475357301444122+1247/24*eps^(-3)+6/5*Pi^4+(125.67152533053854918)*eps^(-2)+3*eps^(-5)+33/2*eps^(-4)+(259.98755698571087874)*eps^(-1)-110/3*zeta(3)
  FRESULT num:           = 347.3630251884288798+(51.958333333333333332)*eps^(-3)+(125.67152533053854918)*eps^(-2)+(3.0)*eps^(-5)+(16.5)*eps^(-4)+(259.98755698571087874)*eps^(-1)
   eps^-5 term: 3 +/- 0
    eps^-4 term: 33/2 +/- 0
     eps^-3 term: 1247/24 +/- 0
      eps^-2 term: 125.67152533053854918 +/- 5.2760713655226570643E-5
       eps^-1 term: 259.98755698571087874 +/- 9.888628922902401464E-6
        eps^0 term: 274.5475357301444122+6/5*Pi^4-110/3*zeta(3) +/- 0.043609817405085687474
        
*/
//     RoMB_loop_by_loop sunset_bub5(lst(l3,k,l1,l2,l4), lst(-pow(l3,2)+ms,-pow(k,2)+ms,-pow(l1,2)+ms,-pow(l2,2)+ms,-pow(l4,2)+ms,-pow(k+l1+l2+l3+l4,2)+ms),lst(pow(p,2)==s),lst(1,1,1,1,1,1));
//     sunset_bub5.integrate_map(lst(ms==1,m2s==1,m3s==1,s==0),1);



//    RoMB_loop_by_loop sunset_bubC2(lst(l1,l2,l3,l4,l5), lst(-pow(l3,2)+ms,-pow(l5,2)+ms,-pow(l1,2)+ms,-pow(l2,2)+ms,-pow(l4,2)+ms,-pow(l5+l1+l2,2)+ms,-pow(l5+l3+l4,2)+ms),lst(pow(p,2)==s),lst(1,1,1,1,1,1,1));
//     sunset_bubC2.integrate_map(lst(ms==1,m2s==1,m3s==1,s==0),0);


//    RoMB_loop_by_loop sunset_bubC2(lst(l1,l2,l3,l4,l5), lst(-pow(l3,2)+ms,-pow(l5,2)+ms,-pow(l1,2)+ms,-pow(l2,2)+ms,-pow(l4,2)+ms,-pow(l3+l4+l5,2)+ms,-pow(l1+l2+l5+l3+l4,2)+ms),lst(pow(p,2)==s),lst(1,1,1,1,1,1,1));
//     sunset_bubC2.integrate_map(lst(ms==1,m2s==1,m3s==1,s==0),0);


//    RoMB_loop_by_loop sunset_bubC1(lst(l1,l2,l3,l4,l5), lst(-pow(l3,2)+ms,-pow(l5,2)+ms,-pow(l1,2)+ms,-pow(l2,2)+ms,-pow(l4,2)+ms,-pow(l2+l3+l4+l5,2)+ms,-pow(l1+l2+l5+l3+l4,2)+ms),lst(pow(p,2)==s),lst(1,1,1,1,1,1,1));
//     sunset_bubC1.integrate_map(lst(ms==1,m2s==1,m3s==1,s==0),0);



/*
MEGA 5-LOOP BUBBLE with 12 propagators
*/

//RoMB_loop_by_loop l5p12(lst(k5,k2,k1,k4,k3),lst(-pow(k1,2)+ms,-pow(k2,2)+ms,-pow(k3,2)+ms,-pow(k4,2)+ms,-pow(k5,2)+ms,
//-pow(k1-k3,2)+ms,-pow(k1-k4,2)+ms,-pow(k3-k2,2)+ms,-pow(k2-k4,2)+ms,-pow(k5+k3-k1,2)+ms,-pow(k5+k3-k2,2)+ms,-pow(k5+k3-k4,2)+ms),lst(pow(p,2)==0),lst(1,1,1,1,1,1,1,1,1,1,1,1));
//l5p12.integrate(lst(ms==0),0);      

//RoMB_loop_by_loop l5c1(lst(k2,k5,k3,k4,k1),lst(-pow(k3,2)+ms,-pow(k2,2)+ms,-pow(k1,2)+ms,-pow(k4,2)+ms,-pow(k5,2)+ms,
//-pow(k1+k3+k4,2)+ms,-pow(k2+k5-k3-k4,2)+ms),lst(pow(p,2)==0),lst(1,1,1,1,1,1,1));
//l5c1.integrate_map(lst(ms==1),0);      



    // RoMB_loop_by_loop t2loop(lst(l,k), lst(-pow(k,2)+ms,-pow(p+k,2)+ms,-pow(p+k+l,2)+ms,-pow(k+l,2)+ms,-pow(l,2)+ms),lst(pow(p,2)==s,ms==0),lst(1,1,1,1,1));
    //t2loop.integrate_map(lst(s==-1,ms == 0),3);
    
    /*     RoMB_loop_by_loop bubble_five_loop(lst(k,l1,l2,l3,l4), 
	   lst(pow(k,2)-ms,pow(l1,2)-ms,pow(l2,2)-ms,pow(l3,2)-ms,pow(l4,2)-ms,pow(k+l1,2)-ms,pow(k+l1+l2,2)-ms,pow(k+l1+l2+l3,2)-ms,pow(k+l1+l2+l3+l4,2)-ms,pow(k+l1+l2+l3,2)-ms,pow(k+l1+l2,2)-ms,pow(k+l1,2)-ms),
	   lst(ms==1),
	   lst(1,1,1,1,1,1,1,1,1,1,1,1));
    */

    // works!!!
    //             RoMB_loop_by_loop B0_1loop_lbl(lst(k),lst(pow(k,2)-2-ms,pow(p+k,2)-ms),lst(ms==0,pow(p,2)==1),lst(2,1));
    
    //     RoMB_loop_by_loop B0_1loop_lbl(lst(k),lst(pow(k,2)-m1s,pow(p+k,2)-m2s),lst(pow(p,2)==s),lst(1,1));
    //    B0_1loop_lbl.integrate(lst(s==-1,m1s==1,m2s==1));
    
    //MB works???
    //                                RoMB_loop_by_loop C0_1loop_lbl(lst(k),lst(pow(k,2),pow(k+p1,2)-m1s,pow(k-p2,2)-m2s),lst(ms==1,pow(p1,2)==m1s,pow(p2,2)==m2s,p1*p2==(s-m1s-m2s)/2),lst(1,1,1));
    //          C0_1loop_lbl.integrate(lst(m1s==1,m2s==1,s==-100));
    

  //MB works???
  
    /*
          RoMB_loop_by_loop box1loopm0(lst(k),lst(-pow(k,2),-pow(k+p1,2),-pow(k+p1+p2,2),-pow(k+p1+p2+p4,2)),
                                  lst(pow(p1,2)==0,pow(p2,2)==0,pow(p4,2)==0,
                                      p1*p2==-s/2,//
                                      
                                      p1*p4==s/2+t/2,//
                                      
                                      p2*p4==-t/2 //
                                      ),
                                       lst(1,1,1,1),false);
    box1loopm0.integrate_map(lst(s==3,t==1));
   
    box1loopm0.integrate(lst(s==5,t==2));
    */
//MASIVE BOX LBL
    /*   
    RoMB_loop_by_loop box1loopm(lst(k),lst(-pow(k,2)+ms,-pow(k+p1,2)+ms,-pow(k+p1+p2,2)+ms,-pow(k+p1+p2+p4,2)+ms),
                                    lst(pow(p1,2)==0,pow(p2,2)==0,pow(p4,2)==0,
                                       p1*p2==-s/2,//
                                        
                                      p1*p4==(s/2+t/2),//
                                        
                                        p2*p4==-t/2 //
                                        ),
                                lst(1,1,1,1),false);
        box1loopm.integrate_map(lst(ms1==1,ms2==1,ms3==1,ms4==1,ms==1,s==30,t==5));
    */

    //triple box
/*
    RoMB_loop_by_loop tribox1loopm(lst(k1,k2,k3),lst(-pow(k1,2)+ms,-pow(k1+p1,2),-pow(k1+p1+p2,2)+ms,
                                                     -pow(k1-k2,2),-pow(k2,2)+ms,-pow(k2+p1+p2,2)+ms,
                                                     -pow(k2-k3,2),-pow(k3,2)+ms,-pow(k3+p1+p2,2)+ms,
                                                     -pow(k3-p3,2)),

                                   lst(pow(p1,2)==ms,pow(p2,2)==ms,pow(p3,2)==ms,pow(p4,2)==ms,
                                       p1*p2==s/2-ms,//
                                        
                                      p1*p3==t/2-ms,//
                                        
                                       p2*p3==ms-(s+t)/2 //
                                        ),
                                   lst(1,1,1,1,1,1,1,1,1,1),true);
    tribox1loopm.integrate_map(lst(ms1==1,ms2==1,ms3==1,ms4==1,ms==1,s==-1/2,t==-3));
*/



    //double box
/*
FRESULT for parameters: {ms1==1,ms2==1,ms3==1,ms4==1,ms==1,s==-2,t==-3}

 FRESULT anl :           = -1.1308856615264005236+(0.09635412596108643146)*eps^(-2)-(0.106390522306243254236)*eps^(-1)
  FRESULT num:           = -1.1308856615264005236+(0.09635412596108643146)*eps^(-2)-(0.106390522306243254236)*eps^(-1)
   eps^-2 term: 0.09635412596108643146 +/- 9.516754647190050904E-6
    eps^-1 term: -0.106390522306243254236 +/- 6.3521193751299645005E-5
     eps^0 term: -1.1308856615264005236 +/- 0.0015177257237448666386
     
     */    
    RoMB_loop_by_loop dobox1loopm(lst(k2,k1),lst(-pow(k1,2)+ms,-pow(k1+p1,2),-pow(k1+p1+p2,2)+ms,
                                                     -pow(k1-k2,2),-pow(k2,2)+ms,-pow(k2+p1+p2,2)+ms,
                                                     -pow(k2-p3,2)),

                                   lst(pow(p1,2)==ms,pow(p2,2)==ms,pow(p3,2)==ms,pow(p4,2)==ms,
                                       p1*p2==s/2-ms,//
                                        
                                      p1*p3==t/2-ms,//
                                        
                                       p2*p3==-(s+t)/2+ms //
                                        ),
                                   lst(1,1,1,1,1,1,1),true);
   dobox1loopm.integrate_map(lst(ms1==1,ms2==1,ms3==1,ms4==1,ms==1,s==-2,t==-3));
    
        
    /*
      4-loop  tadpole
    */

       /*
         RoMB_loop_by_loop tad4(lst(l1, l2, l3, l4),lst(-pow(l1,2)+ ms,-pow(l2,2)+ ms,-pow(l3 ,2)+ ms,-pow(l4,2),-pow(l1+l2+l3+l4,2)),lst(),lst(1,1,1,1,1));
         tad4.integrate_map(lst(ms == 1),1);

         */
    
     /*
       Pentagon
     */
  /* 
            RoMB_loop_by_loop pent(lst(k1),lst(-pow(p1 + k1,2)+ ms,-pow(p1 + p5 + k1,2),
                                        -pow(p1 + p5 + p4 + k1,2)+ ms,-pow(p1 + p5 + p4 + p3 + k1,2)+ ms,
                                        -pow(k1,2)),
                            lst(



p1*p1 == ms, p2*p2 == ms, p3*p3 == 0, p4*p4 == ms, p5*p5 == ms, 
p1*p2 == 1/2* (tp - 2* ms), p1*p3 == 1/2* (t - tp - v1), 
p1*p4 == ms - 1/2* (s + t - v1), p1*p5 == 1/2* (s - 2* ms), 
p2* p3 == 1/2* v1, p2* p4 == 1/2* (s - 2* ms - v1 - v2), 
p2* p5 == ms - 1/2* (s + tp - v2), p3* p4 == 1/2* v2, 
p3* p5 == 1/2* (tp - t - v2), p4* p5 == 1/2* (t - 2* ms)),
lst(1,1,1,1,1));
pent.integrate_map(lst(s==-2,t==-3,v2==-4,tp==-5,v1==-6,ms==1));

*/
    /*
        RoMB_loop_by_loop pent(lst(k1),lst(-pow(p1 + k1,2)+ ms,-pow(p1 + p5 + k1,2),
                                           -pow(p1 + p5 + p4 + k1,2)+ ms,-pow(p1 + p5 + p4 + p3 + k1,2)+ ms,
                                           -pow(k1,2)),
                               lst(
                                   p1*p1 == ms, 
                                   p2*p2 == ms, 
                                   p3*p3 == 0, 
                                   p4*p4 == ms, 
                                   p5*p5 == ms, 
                                   p1*p2 == 1/2* (tp - 2* ms), 
                                   wild(1)*p1*p3 == wild(1)*1/2* (t - tp - v1), 
                                   wild(2)*p1*p4 == wild(2)*(ms - 1/2* (s + t - v1)), 
                                  wild(3)* p1*p5 == wild(3)*1/2* (s - 2* ms), 
                                   wild(4)*p2* p3 == wild(4)*1/2* v1, 
                                   wild(5)*p2* p4 == wild(5)*1/2* (s - 2* ms - v1 - v2), 
                                   wild(6)*p2* p5 ==wild(6)*( ms - 1/2* (s + tp - v2)), 
                                   wild(7)*p3* p4 == wild(7)*1/2* v2, 
                                   wild(8)*p3* p5 == wild(8)*1/2* (tp - t - v2), 
                                   wild()*p4* p5 == wild()*1/2* (t - 2* ms)),
                               lst(1,1,1,1,1));
        pent.integrate_map(lst(s==-2,t==-3,v2==-4,tp==-5,v1==-6,ms==1));
       
        */
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
