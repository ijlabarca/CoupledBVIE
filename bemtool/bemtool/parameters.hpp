#ifndef BEMTOOL_PARAMETERS_HPP
#define BEMTOOL_PARAMETERS_HPP

#include <fstream>
#include "tools.hpp"

using namespace bemtool;


Real epsilon0 = 1.;
Real epsilon1 = std::sqrt(2.);

Real mu0 = 1.;
Real mu1 = std::sqrt(2.);


Real kappa0  = std::sqrt(epsilon0 * mu0);
Real kappa2 = kappa0*kappa0;


Real kappa1  = std::sqrt(epsilon1 * mu1);


// ----------------------------------------------------------//
// ------------ MAXWELL TRANSMISSION PROBLEM COEFFS ---------//
// ----------------------------------------------------------//

// Real epsilon(R3 x) {return epsilon1;}
// Real mu(R3 x) {return mu1;}

// Real epsilon(R3 x) {return epsilon1 + 4.*(x[0]*(1-x[0]))*(x[1]*(1-x[1]))*(x[2]*(1-x[2]));}
// Real mu(R3 x)      {return mu1      + 4.*(x[0]*(1-x[0]))*(x[1]*(1-x[1]))*(x[2]*(1-x[2]));}

// Real epsilon(R3 x) {return epsilon1 + 8./3.*((x[0]*x[0])+(x[1]*x[1])+(x[2]*x[2]));}
// Real mu(R3 x)      {return mu1      + 0.*(x[0]*(1-x[0]))*(x[1]*(1-x[1]))*(x[2]*(1-x[2]));}

Real epsilon(R3 x) {return epsilon1 + 1.*(1+x[0])*(1-x[0])*(1+x[1])*(1-x[1])*(1+x[2])*(1-x[2]);}
Real mu(R3 x)      {return mu1      + 1.*(1+x[0])*(1-x[0])*(1+x[1])*(1-x[1])*(1+x[2])*(1-x[2]);}


Real invmu(R3 x) {return 1./mu(x);}

/*
// */
// R3 gradeps(R3 x) {R3 val;val[0] = 0.;
//                          val[1] = 0.;
//                          val[2] = 0.;
//                   return val;}
                  
                  
// R3 gradmu(R3 x) {R3 val;val[0] = 0.;
//                          val[1] = 0.;
//                          val[2] = 0.;
//                   return val;}
                  
R3 gradeps(R3 x) {R3 val;val[0] =  1.*(-2.*x[0])*(1+x[1])*(1-x[1])*(1+x[2])*(1-x[2]);
                          val[1] = 1.*(1+x[0])*(1-x[0])*(-2.*x[1])*(1+x[2])*(1-x[2]);
                          val[2] = 1.*(1+x[0])*(1-x[0])*(1+x[1])*(1-x[1])*(-2.*x[2]);
                   return val;}


R3 gradmu(R3 x) {R3 val;val[0] =   1.*(-2.*x[0])*(1+x[1])*(1-x[1])*(1+x[2])*(1-x[2]);
                          val[1] = 1.*(1+x[0])*(1-x[0])*(-2.*x[1])*(1+x[2])*(1-x[2]);
                          val[2] = 1.*(1+x[0])*(1-x[0])*(1+x[1])*(1-x[1])*(-2.*x[2]);
                   return val;}





// R3 gradeps(R3 x) {R3 val;val[0] = 4.*(1-2.*x[0])*(x[1]*(1-x[1]))*(x[2]*(1-x[2]));
//                          val[1] = 4.*(x[0]*(1-x[0]))*(1-2.*x[1])*(x[2]*(1-x[2]));
//                          val[2] = 4.*(x[0]*(1-x[0]))*(x[1]*(1-x[1]))*(1-2.*x[2]);
//                   return val;}


// R3 gradmu(R3 x) {R3 val; val[0] = 4.*(1-2.*x[0])*(x[1]*(1-x[1]))*(x[2]*(1-x[2]));
//                          val[1] = 4.*(x[0]*(1-x[0]))*(1-2.*x[1])*(x[2]*(1-x[2]));
//                          val[2] = 4.*(x[0]*(1-x[0]))*(x[1]*(1-x[1]))*(1-2.*x[2]);
//                   return val;}





// R3 gradeps(R3 x) {R3 val;val[0] = 16./3.*x[0];
//                          val[1] = 16./3.*x[1];
//                          val[2] = 16./3.*x[2];
//                   return val;}


// R3 gradmu(R3 x) {R3 val; val[0] = 0.*(1-2.*x[0])*(x[1]*(1-x[1]))*(x[2]*(1-x[2]));
//                          val[1] = 0.*(x[0]*(1-x[0]))*(1-2.*x[1])*(x[2]*(1-x[2]));
//                          val[2] = 0.*(x[0]*(1-x[0]))*(x[1]*(1-x[1]))*(1-2.*x[2]);
//                   return val;}

Real epstilde(R3 x) {return epsilon(x)/epsilon1;}
Real invepstilde(R3 x) {return epsilon1/epsilon(x);}


Real mutilde(R3 x) {return mu(x)/mu1;}
Real invmutilde(R3 x) {return mu1/mu(x);}

Real p_e(R3 x) {return (1. - epsilon(x)/epsilon1);}
Real p_m(R3 x) {return (1. - mu(x)/mu1);}


R3 tau_e(R3 x) {R3 val = gradeps(x); R3 val2;
                 val2[0] = val[0] / epsilon(x);
                 val2[1] = val[1] / epsilon(x);
                 val2[2] = val[2] / epsilon(x);
          return val2;}
               
               
R3 tau_m(R3 x) {R3 val = gradmu(x); R3 val2;
                 val2[0] = val[0] / mu(x);
                 val2[1] = val[1] / mu(x);
                 val2[2] = val[2] / mu(x);
          return val2;}

R3 gradp_e(R3 x) {R3 val = gradeps(x); R3 val2;
                       val2[0] = -val[0]/epsilon1;
                       val2[1] = -val[1]/epsilon1;
                       val2[2] = -val[2]/epsilon1;
               return val2;}
               
R3 gradp_m(R3 x) {R3 val = gradmu(x); R3 val2;
       val2[0] = -val[0]/mu1;
       val2[1] = -val[1]/mu1;
       val2[2] = -val[2]/mu1;
return val2;}
               

bool testEps(const Mesh3D& mesh){

// This function returns: 
// if constant eps is constant : 1
// else                        : 0

      int nb_elt = NbElt(mesh);

      int count = 0;
      
      for(int n = 0; n < nb_elt; n++){
        
        Elt3D elt = mesh[n];
        
        R3 X = Ctr(elt);
        
        if(std::abs(epsilon1 - epsilon(X)) < 1e-12){
          count++;
        }
        else{
          break;
        }  
      }
      
      return (count == nb_elt);
}
        
bool testEps(const Mesh2D& mesh){
// This function returns: 
// if constant eps is constant : 1
// else                        : 0


      int nb_elt = NbElt(mesh);
      int count = 0;
      
      for(int n = 0; n < nb_elt; n++){
        
        Elt2D elt = mesh[n];
        
        R3 X = Ctr(elt);
        
        if(std::abs(epsilon1 - epsilon(X)) < 1e-12){
          count++;
        }
        else{
          break;
        }  
      }
      
      return (count == nb_elt);
}       



           

bool testMu(const Mesh3D& mesh){
// This function returns: 
// if eps is constant : 1
// else               : 0


      int nb_elt = NbElt(mesh);
      int count = 0;
      
      for(int n = 0; n < nb_elt; n++){
        
        Elt3D elt = mesh[n];
        
        R3 X = Ctr(elt);
        
        if(std::abs(mu1 - mu(X)) < 1e-12){
          count++;
        }
        else{
          break;
        }  
      }
      
      return (count == nb_elt);
}
        
bool testMu(const Mesh2D& mesh){
// This function returns: 
// if mu is constant : 1
// else              : 0


      int nb_elt = NbElt(mesh);
      int count = 0;
      
      for(int n = 0; n < nb_elt; n++){
        
        Elt2D elt = mesh[n];
        
        R3 X = Ctr(elt);
        
        if(std::abs(mu1 - mu(X)) < 1e-12){
          count++;
        }
        else{
          break;
        }  
      
      }
      return (count == nb_elt);
}       
               
// ----------------------------------------------------------//
// ------------INCIDENT PLANEWAVE----------------------------//
// ----------------------------------------------------------//
Cplx f(R3 x) {return std::exp(iu * kappa0 * (-x[2]));}

Real fx(R3 x) {return 1.0;}
Real fy(R3 x) {return 1.0;}
Real fz(R3 x) {return 1.0;}


Cplx Ex(R3 x) {return  0. * iu;}
Cplx Ey(R3 x) {return  f(x);}
Cplx Ez(R3 x) {return  0. * iu;}

Cplx Hx(R3 x) {return  f(x);}
Cplx Hy(R3 x) {return  0. * iu;}
Cplx Hz(R3 x) {return  0. * iu;}

R3 func(R3 x) {R3 val; val[0] = fx(x);
                       val[1] = fy(x);
                       val[2] = fz(x);
               return val;}

// ----------------------------------------------------------//


int CheckInteraction(const int& jx, const int& jy, const Mesh3D& mesh){

    const int dimx = 3;
    const int dimy = 3;
    int rule = 0;
    

    Elt3D ex = mesh[jx];   Elt3D ey = mesh[jy];

    for(int p=0; p<dimx+1; p++){
      for(int q=rule; q<dimy+1; q++){
        if( &ex[p]==&ey[q] ){
          Swap(ex,rule,p); 
          Swap(ey,rule,q); 
          rule++; break;
        }
      }
    }

    return rule;

}



int CheckInteraction(const int& jx, const int& jy, const Mesh2D& mesh){

    const int dimx = 2;
    const int dimy = 2;
    int rule = 0;
    

    Elt2D ex = mesh[jx];   Elt2D ey = mesh[jy];

    for(int p=0; p<dimx+1; p++){
      for(int q=rule; q<dimy+1; q++){
        if( &ex[p]==&ey[q] ){
          Swap(ex,rule,p); 
          Swap(ey,rule,q); 
          rule++; break;
        }
      }
    }

    return rule;

}



int CheckInteraction(const int& jx, const int& jy, const Mesh3D& mesh3D, const Mesh2D& mesh2D){

    const int dimx = 3;
    const int dimy = 2;
    int rule = 0;
    

    Elt3D ex = mesh3D[jx];   Elt2D ey = mesh2D[jy];

    for(int p=0; p<dimx+1; p++){
      for(int q=rule; q<dimy+1; q++){
        if( &ex[p]==&ey[q] ){
          Swap(ex,rule,p); 
          Swap(ey,rule,q); 
          rule++; break;
        }
      }
    }

    return rule;

}


int CheckInteraction(const int& jx, const int& jy, const Mesh2D& mesh2D, const Mesh3D& mesh3D){

    const int dimx = 2;
    const int dimy = 3;
    int rule = 0;
    

    Elt2D ex = mesh2D[jx];   Elt3D ey = mesh3D[jy];

    for(int p=0; p<dimx+1; p++){
      for(int q=rule; q<dimy+1; q++){
        if( &ex[p]==&ey[q] ){
          Swap(ex,rule,p); 
          Swap(ey,rule,q); 
          rule++; break;
        }
      }
    }

    return rule;

}

#endif
