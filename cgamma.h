#ifndef __CGAMMA_H__
#define __CGAMMA_H__
#include <complex>
namespace RoMB
{
  //
  //  Returns gamma function or log(gamma) for complex argument 'z'.
  //
  //  OPT value       function
  //  ---------       --------
  //      0           complex gamma
  //      1           complex log(gamma)
  //
  //  Returns (1e308,0) if the real part of the argument is a negative integer
  //  or 0 or exceeds 171.
  //
  std::complex<double> cgamma(std::complex<double> z,int OPT = 0)
    {
      std::complex<double> g,z0,z1;
      double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
      double na,t,x1,y1,sr,si;
      int i,j,k;

      static double a[] = {
        8.333333333333333e-02,
        -2.777777777777778e-03,
        7.936507936507937e-04,
        -5.952380952380952e-04,
        8.417508417508418e-04,
        -1.917526917526918e-03,
        6.410256410256410e-03,
        -2.955065359477124e-02,
        1.796443723688307e-01,
        -1.39243221690590};

      x = real(z);
      y = imag(z);
      if (x > 171) return std::complex<double>(1e308,0);
      if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
        return std::complex<double>(1e308,0);
      else if (x < 0.0) {
        x1 = x;
        y1 = y;
        x = -x;
        y = -y;
      }
      x0 = x;
      if (x <= 7.0) {
        na = (int)(7.0-x);
        x0 = x+na;
      }
      q1 = sqrt(x0*x0+y*y);
      th = atan(y/x0);
      gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
      gi = th*(x0-0.5)+y*log(q1)-y;
      for (k=0;k<10;k++){
        t = pow(q1,-1.0-2.0*k);
        gr += (a[k]*t*cos((2.0*k+1.0)*th));
        gi -= (a[k]*t*sin((2.0*k+1.0)*th));
      }
      if (x <= 7.0) {
        gr1 = 0.0;
        gi1 = 0.0;
        for (j=0;j<na;j++) {
          gr1 += (0.5*log((x+j)*(x+j)+y*y));
          gi1 += atan(y/(x+j));
        }
        gr -= gr1;
        gi -= gi1;
      }
      if (x1 <= 0.0) {
        q1 = sqrt(x*x+y*y);
        th1 = atan(y/x);
        sr = -sin(M_PI*x)*cosh(M_PI*y);
        si = -cos(M_PI*x)*sinh(M_PI*y);
        q2 = sqrt(sr*sr+si*si);
        th2 = atan(si/sr);
        if (sr < 0.0) th2 += M_PI;
        gr = log(M_PI/(q1*q2))-gr;
        gi = -th1-th2-gi;
        x = x1;
        y = y1;
      }
      if (OPT == 0) {
        g0 = exp(gr);
        gr = g0*cos(gi);
        gi = g0*sin(gi);
      }
      g = std::complex<double>(gr,gi);
      return g;
    }

  std::complex<double> cpsi(std::complex<double> z)
    {
      double x,y,x0,x1,y1,psr,psi;
      double q0,q2,rr,ri,th,tn,tm,ct2;
      std::complex<double>ps;

      int n,k;
      static double a[] = {
        -0.8333333333333e-01,
        0.83333333333333333e-02,
        -0.39682539682539683e-02,
        0.41666666666666667e-02,
        -0.75757575757575758e-02,
        0.21092796092796093e-01,
        -0.83333333333333333e-01,
        0.4432598039215686};

      x = real(z);
      y = imag(z);
      x1 = 0.0;
      n = 0;
      if ((y == 0.0) && (x == (int)x) && (x <= 0.0)) {
        ps = std::complex<double> (1e308,0);
      }
      else {
        if (x < 0.0) {
          x1 = x;
          y1 = y;
          x = -x;
          y = -y;
        }
        x0 = x;
        if (x < 8.0) {
          n = 8 - (int)x;
          x0 = x + n;
        }
        if ((x0 == 0.0) && (y != 0.0))
          th = 0.5*M_PI;
        if (x0 != 0.0)
          th = atan(y/x0);
        q2 = x0*x0+y*y;
        q0 = sqrt(q2);
        psr = log(q0)-0.5*x0/q2;
        psi = th+0.5*y/q2;
        for (k=1;k<=8;k++) {
          psr += (a[k-1]*pow(q2,-k)*cos(2.0*k*th));
          psi -= (a[k-1]*pow(q2,-k)*sin(2.0*k*th));
        }
        if (x < 8.0) {
          rr = 0.0;
          ri = 0.0;
          for (k=1;k<=n;k++) {
            rr += ((x0-k)/(pow(x0-k,2.0)+y*y));
            ri += (y/(pow(x0-k,2.0)+y*y));
          }
          psr -= rr;
          psi += ri;
        }
        if (x1 < 0.0) {
          tn = tan(M_PI*x);
          tm = tanh(M_PI*y);
          ct2 = tn*tn+tm*tm;
          psr = psr+x/(x*x+y*y)+M_PI*(tn-tn*tm*tm)/ct2;
          psi = psi-y/(x*x+y*y)-M_PI*tm*(1.0+tn*tn)/ct2;
          x = x1;
          y = y1;
        }
        ps = std::complex<double>(psr,psi);
      }
      return ps;
    }
  // Helper functions for real and imagenary 
  // parts of double and complex expressions
  double Real(std::complex<double> c_in)
  {
    return c_in.real();
  }
  double Real(double d_in)
  {
    return d_in;
  }
  double Imag(std::complex<double> c_in)
  {
    return c_in.imag();
  }
  double Imag(double d_in)
  {
    return 0;
  }


}; // namespace RoMB

std::complex<double> psi(std::complex<double> z)
{
  double x,y,x0,x1,y1,psr,psi;
  double q0,q2,rr,ri,th,tn,tm,ct2;
  std::complex<double>ps;

  int n,k;
  static double a[] = {
    -0.8333333333333e-01,
    0.83333333333333333e-02,
    -0.39682539682539683e-02,
    0.41666666666666667e-02,
    -0.75757575757575758e-02,
    0.21092796092796093e-01,
    -0.83333333333333333e-01,
    0.4432598039215686};

  x = real(z);
  y = imag(z);
  x1 = 0.0;
  n = 0;
  if ((y == 0.0) && (x == (int)x) && (x <= 0.0)) {
    ps = std::complex<double> (1e308,0);
  }
  else {
    if (x < 0.0) {
      x1 = x;
      y1 = y;
      x = -x;
      y = -y;
    }
    x0 = x;
    if (x < 8.0) {
      n = 8 - (int)x;
      x0 = x + n;
    }
    if ((x0 == 0.0) && (y != 0.0))
      th = 0.5*M_PI;
    if (x0 != 0.0)
      th = atan(y/x0);
    q2 = x0*x0+y*y;
    q0 = sqrt(q2);
    psr = log(q0)-0.5*x0/q2;
    psi = th+0.5*y/q2;
    for (k=1;k<=8;k++) {
      psr += (a[k-1]*pow(q2,-k)*cos(2.0*k*th));
      psi -= (a[k-1]*pow(q2,-k)*sin(2.0*k*th));
    }
    if (x < 8.0) {
      rr = 0.0;
      ri = 0.0;
      for (k=1;k<=n;k++) {
        rr += ((x0-k)/(pow(x0-k,2.0)+y*y));
        ri += (y/(pow(x0-k,2.0)+y*y));
      }
      psr -= rr;
      psi += ri;
    }
    if (x1 < 0.0) {
      tn = tan(M_PI*x);
      tm = tanh(M_PI*y);
      ct2 = tn*tn+tm*tm;
      psr = psr+x/(x*x+y*y)+M_PI*(tn-tn*tm*tm)/ct2;
      psi = psi-y/(x*x+y*y)-M_PI*tm*(1.0+tn*tn)/ct2;
      x = x1;
      y = y1;
    }
    ps = std::complex<double>(psr,psi);
  }
  return ps;
}

std::complex<double> psi(double div,std::complex<double> z)
{
  // trigamma function
  if(1 == int(div))
    {
      double x,y,x_asmpt,y_asmpt,x0,y0,x1,y1,psr,psi;
      double q0,q2,rr,ri,th,tn,tm,ct2;
      std::complex<double>ps;

      int n,k;
      static double a[] = {
        0.1666666666666667,  // B_2
        -0.3333333333333333e-01,// B_4
        0.2380952380952381e-01, //B_6
        -0.3333333333333333e-01,// B_8
        0.7575757575757576e-01,//B_10
        -0.2531135531135531,//B_12
        1.166666666666667,//B_14
        -7.092156862745098,//B_16
        54.97117794486216,//B_18
        -529.1242424242424,//B_20
        6192.123188405797};//B_22

      x = real(z);
      y = imag(z);
      x_asmpt = 0.0; // point for 
      y_asmpt = 0.0; // assymptotic formula application
      x1 = 0.0;
      n = 0;
      if ((y == 0.0) && (x == (int)x) && (x <= 0.0)) 
        {
          ps = std::complex<double> (1e308,0);
        }
      else 
        {
          if(x < 0.0)
            {
              x0 = x;
              y0 = y;
              x = -x0;
              y = -y0;
            }
          x_asmpt = x;
          y_asmpt = y;
          if ((0.0 <= x) &&(x < 10.0)) 
            {
              n = 10 - (int)x;
              x_asmpt = x + n;
            }

          std::complex<double> z_asmpt(x_asmpt,y_asmpt);
          std::complex<double> r = pow(z_asmpt,-1)+0.5*pow(z_asmpt,-2);
          for(int m = 1;m < 12; m++)
            {
              r+= a[m-1]*pow(z_asmpt,-(2*m+1));
            }
          for(int k = 1;k <= n;k++)
            {
              r+=pow(z_asmpt - std::complex<double>(double(k),0.0),-2);
            }
          if(x0 < 0.0)
            {
              ps = -r + pow(std::complex<double>(x0,y0),-2) + M_PI*M_PI*pow(sin(M_PI*std::complex<double>(x0,y0)),-2);
            }
          else
            ps=r;
        }
      return ps;
    }
}


double const Pi = M_PI; //4*atan(1); // Pi=3.14
double const Euler = 0.577215664901532860606;  //  \gamma_E = 0.577215664901532860606 
#define tgamma RoMB::cgamma                                                                                                                                        
//#define psi RoMB::cpsi
#endif
