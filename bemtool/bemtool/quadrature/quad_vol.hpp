#ifndef BEMTOOL_QUADRATURE_QUADVOL_HPP
#define BEMTOOL_QUADRATURE_QUADVOL_HPP

#include "quad.hpp"
#include <fstream>

namespace bemtool{



  template <int,int> class QuadVol;


  template <> class QuadVol<3,3>{

  public:
    typedef R3   qp_t;

  private:
    // arrays of std::vectors
    std::vector<R3>   x_[5]; // 5 cases: disjoint, common vertex, edge, face, same
    std::vector<R3>   y_[5];
    std::vector<Real> w_[5];

  public:
    const std::vector<R3>&   x(const int& rule) const {return x_[rule];}
    const std::vector<R3>&   y(const int& rule) const {return y_[rule];}
    const std::vector<Real>& w(const int& rule) const {return w_[rule];}

    // Constructor
    QuadVol<3,3>(const int& order){

      R3 x,y; Real w,dw,dww;
      std::vector<Real> t,dt, tt, dtt;
      Quad1D(order,t,dt);
      int nt = t.size();
      int ntt = tt.size();

      


    //================================//
    //   Cas 0: triangles disjoints   //
    //================================//
      std::vector<R3>   qp;
      std::vector<Real> qw;
      Quad3D(4,qp,qw);
      int nq = qp.size();

    for(int j=0; j<nq; j++){
      for(int k=0; k<nq; k++){

 	x[0] = qp[j][0];
	x[1] = qp[j][1];
	x[2] = qp[j][2];
	y[0] = qp[k][0];
	y[1] = qp[k][1];
  y[2] = qp[k][2];	
	w    = qw[j]*qw[k];

  // Case 0 -> index 0
	x_[0].push_back(x);
	y_[0].push_back(y);
	w_[0].push_back(w);

      }
    }
      
    // for(int j0=0; j0<ntt; j0++){
    //   for(int j1=0; j1<ntt; j1++){
	  //     for(int j2=0; j2<ntt; j2++){
	  //       for(int j3=0; j3<ntt; j3++){
	  //         for(int j4=0; j4<ntt; j4++){
	  //           for(int j5=0; j5<ntt; j5++){
	                        
                          
    //           //================================//
    //           //         Case 0: Disjoint       //
    //           //================================//
    //             Real ee1  = tt[j0];
    //             Real ee2  = tt[j1];
    //             Real ee3  = tt[j2];
    //             Real ee4  = tt[j3];
    //             Real ee5  = tt[j4];
    //             Real ee6  = tt[j5];
                
                
    //       	    dww = dtt[j0]*dtt[j1]*dtt[j2]*dtt[j3]*dtt[j4]*dtt[j5];


    //             x[0] = ee1*ee2*ee3;
    //             x[1] = ee1*ee2*(1-ee3);
    //             x[2] = ee1*(1-ee2);
                
    //             y[0] = ee4*ee5*ee6;
    //             y[1] = ee4*ee5*(1-ee6);
    //             y[2] = ee4*(1-ee5);
                
    //             w    = ee1*ee1*ee2*ee4*ee4*ee5*dww;

    //             x_[0].push_back(x);
    //             y_[0].push_back(y);
    //             w_[0].push_back(w);
                
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }
      
      

    for(int j0=0; j0<nt; j0++){
      for(int j1=0; j1<nt; j1++){
	      for(int j2=0; j2<nt; j2++){
	        for(int j3=0; j3<nt; j3++){
	          for(int j4=0; j4<nt; j4++){
	            for(int j5=0; j5<nt; j5++){
	                        
                          
              //================================//
              //         Case 0: Disjoint       //
              //================================//
                Real e1  = t[j0];
                Real e2  = t[j1];
                Real e3  = t[j2];
                Real e4  = t[j3];
                Real e5  = t[j4];
                Real e6  = t[j5];
                
                
          	    // dw = dt[j0]*dt[j1]*dt[j2]*dt[j3]*dt[j4]*dt[j5];

/*
                x[0] = e1*e2*e3;
                x[1] = e1*e2*(1-e3);
                x[2] = e1*(1-e2);
                
                y[0] = e4*e5*e6;
                y[1] = e4*e5*(1-e6);
                y[2] = e4*(1-e5);
                
                w    = e1*e1*e2*e4*e4*e5*dw;

                x_[0].push_back(x);
                y_[0].push_back(y);
                w_[0].push_back(w);
  */              
                
              //================================//
              //      Case 1: Common Vertex     //
              //================================//
              
              /*
                const Real& e1  = t[j0];
                const Real& e2  = t[j1];
                const Real& e3  = t[j2];
                const Real& e4  = t[j3];
                const Real& e5  = t[j4];
                const Real& z   = t[j5];
                */
                
                
                e1  = t[j0];
                e2  = t[j1];
                e3  = t[j2];
                e4  = t[j3];
                e5  = t[j4];
                Real z   = t[j5];
                
                x[0] = z*e1*e2;
                x[1] = z*e1*(1-e2);
                x[2] = z*(1-e1);
                
                y[0] = z*e3*e4*e5;
                y[1] = z*e3*e4*(1-e5);
                y[2] = z*e3*(1-e4);
                
                w    = z*z*z*z*z*e1*e3*e3*e4*dw;
                          
                // Case 1 -> index 1
	              x_[1].push_back(x);
	              y_[1].push_back(y);
	              w_[1].push_back(w);
	              
	              //---------------------------------------//
                
                x[0] = z*e3*e4*e5;
                x[1] = z*e3*e4*(1-e5);
                x[2] = z*e3*(1-e4);
                
                y[0] = z*e1*e2;
                y[1] = z*e1*(1-e2);
                y[2] = z*(1-e1);
                
	              x_[1].push_back(x);
	              y_[1].push_back(y);
	              w_[1].push_back(w);
	              
	              
              //================================//
              //      Case 2: Common Edge       //
              //================================//
              
              /*
                const Real& e1  = t[j0];
                const Real& e2  = t[j1];
                const Real& e3  = t[j2];
                const Real& e4  = t[j3];
                const Real& z1  = t[j4];
                const Real& z2  = t[j5];
                */
                
                e1  = t[j0];
                e2  = t[j1];
                e3  = t[j2];
                e4  = t[j3];
                Real z1  = t[j4];
                Real z2  = t[j5];
                
                
                x[0] = 1-z1;
                x[1] = z1*z2*e1;
                x[2] = z1*z2*(1-e1);
                
                y[0] = 1-z1 + z1*z2*e2*e3*e4;
                y[1] = z1*z2*e2*e3*(1-e4);
                y[2] = z1*z2*e2*(1-e3);
                
                w    = z1*z1*z1*z1*z1*z2*z2*z2*z2*e2*e2*e3*dw;
                          
                // Case 2 -> index 2
	              x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1 + z1*z2*e2*e3*e4;
                x[1] = z1*z2*e2*e3*(1-e4);
                x[2] = z1*z2*e2*(1-e3);
                
                y[0] = 1-z1;
                y[1] = z1*z2*e1;
                y[2] = z1*z2*(1-e1);
                
                w    = z1*z1*z1*z1*z1*z2*z2*z2*z2*e2*e2*e3*dw;
                          
                x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1;
                x[1] = z1*z2*e3*e4;
                x[2] = z1*z2*e3*(1-e4);
                
                y[0] = 1-z1 + z1*z2*e1*e2;
                y[1] = z1*z2*e1*(1-e2);
                y[2] = z1*z2*(1-e1);
                
                w    = z1*z1*z1*z1*z1*z2*z2*z2*z2*e1*e3*dw;
                          
                x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
              
              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1 + z1*z2*e1*e2;
                x[1] = z1*z2*e1*(1-e2);
                x[2] = z1*z2*(1-e1);
                
                y[0] = 1-z1;
                y[1] = z1*z2*e3*e4;
                y[2] = z1*z2*e3*(1-e4);
                
                w    = z1*z1*z1*z1*z1*z2*z2*z2*z2*e1*e3*dw;
                          
                x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
	              
	              
              //================================//
              //      Case 3: Common Face       //
              //================================//
              
              /*
                const Real& e1  = t[j0];
                const Real& e2  = t[j1];
                const Real& e3  = t[j2];
                const Real& z1  = t[j3];
                const Real& z2  = t[j4];
                const Real& z3  = t[j5];
                */
                
                
                e1  = t[j0];
                e2  = t[j1];
                e3  = t[j2];
                z1  = t[j3];
                z2  = t[j4];
                Real z3  = t[j5];
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1;
                x[2] = z1*z2*z3*(1-e1);
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2*e3;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*e1*e2*(1-e3);
                
                w    = z1*z1*z1*z1*z1*z2*z2*z2*z2*z3*z3*z3*e1*e1*e2*dw;
                          
                // Case 3 -> index 3
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                x[2] = z1*z2*z3*e1*(1-e2);
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2*e3;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*(1-e1*e2*e3);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                x[2] = z1*z2*z3*(1-e1*e2);
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2*e3;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*e1*(1-e2*e3);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*e2*e3;
                x[2] = z1*z2*z3*e1*e2*(1-e3);
                
                y[0] = 1-z1 + z1*z2*z3*e1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*(1-e1);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*e2*e3;
                x[2] = z1*z2*z3*e1*(1-e2*e3);
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*(1-e1*e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
                
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*e2*e3;
                x[2] = z1*z2*z3*(1-e1*e2*e3);
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*e1*(1-e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2);
                x[2] = z1*z2*z3;
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2*e3;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2*(1-e3);
                y[2] = z1*z2*z3*e1*(1-e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //========================================//
	              
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2*e3;
                x[1] = z1*(1-z2);
                x[2] = z1*z2*z3*e1*e2*(1-e3);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1;
                y[2] = z1*z2*z3*(1-e1);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2*e3;
                x[1] = z1*(1-z2);
                x[2] = z1*z2*z3*(1-e1*e2*e3);
                
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                y[2] = z1*z2*z3*e1*(1-e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2*e3;
                x[1] = z1*(1-z2);
                x[2] = z1*z2*z3*e1*(1-e2*e3);
                
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                y[2] = z1*z2*z3*(1-e1*e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                
                x[0] = 1-z1 + z1*z2*z3*e1;
                x[1] = z1*(1-z2);
                x[2] = z1*z2*z3*(1-e1);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2*e3;
                y[2] = z1*z2*z3*e1*e2*(1-e3);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2;
                x[1] = z1*(1-z2);
                x[2] = z1*z2*z3*(1-e1*e2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2*e3;
                y[2] = z1*z2*z3*e1*(1-e2*e3);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
                
	              //---------------------------------------//
                
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2;
                x[1] = z1*(1-z2);
                x[2] = z1*z2*z3*e1*(1-e2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2*e3;
                y[2] = z1*z2*z3*(1-e1*e2*e3);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              //---------------------------------------//
                
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2*e3;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*e2*(1-e3);
                x[2] = z1*z2*z3*e1*(1-e2);
                
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3;
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
	              //---------------------------------------//
	              
	              
                x[0] = 1-z1 + z1*z2*z3*e1*e2*e3;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*e2*(1-e3);
                x[2] = z1*z2*z3*(1-e1*e2);
                
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*e1;
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
	              
	              
                x[0] = 1-z1;
                x[1] = z1*(1-z2);
                x[2] = z1*z2*z3*e3;
                
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*(1-e2);
                y[2] = z1*z2*z3*(1-e1);
                
                w    = z1*z1*z1*z1*z1*z2*z2*z2*z2*z3*z3*z3*e1*dw;
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              //---------------------------------------//
	              
	              
                x[0] = 1-z1 + z1*z2*z3*e1*e2;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*(1-e2);
                x[2] = z1*z2*z3*(1-e1);
                
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*e1*e3;
                
                w    = z1*z1*z1*z1*z1*z2*z2*z2*z2*z3*z3*z3*e1*e1*dw;
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
              //================================//
              //      Case 4: Identical         //
              //================================//
              
              /*
                const Real& e1  = t[j0];
                const Real& e2  = t[j1];
                const Real& z1  = t[j2];
                const Real& z2  = t[j3];
                const Real& z3  = t[j4];
                const Real& z4  = t[j5];
                */
                
                e1  = t[j0];
                e2  = t[j1];
                z1  = t[j2];
                z2  = t[j3];
                z3  = t[j4];
                Real z4  = t[j5];
              
                
                Real w1  = z1;
                Real w2  = z1*z2;
                Real w3  = z1*z2*z3;
                Real w4  = z1*z2*z3*z4;
                Real w5  = z1*z2*z3*z4*e1;
                Real w6  = z1*z2*z3*z4*e1*e2;
                
                x[0] = 1-w1;
                x[1] = w1 - w2 + w5;
                x[2] = w3 - w5;
                
                y[0] = 1-w1 +w6;
                y[1] = w1-w2;
                y[2] = w3 - w4;
                
                w    = z1*z1*z1*z1*z1*z2*z2*z2*z2*z3*z3*z3*z4*z4*e1*dw;
                          
                // Case 4 -> index 4
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1;
                x[1] = w1 - w2 + w4;
                x[2] = w3 - w4;
                
                y[0] = 1-w1 +w6;
                y[1] = w1-w2;
                y[2] = w3 - w4 + w5 - w6;
                
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1;
                x[1] = w1 - w2 + w5;
                x[2] = w3 - w4;
                
                y[0] = 1-w1 +w6;
                y[1] = w1-w2;
                y[2] = w3 - w6;
                
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1;
                x[1] = w1 - w2 + w6;
                x[2] = w3 - w6;
                
                y[0] = 1-w1 +w5;
                y[1] = w1-w2;
                y[2] = w3 - w4;
                
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1;
                x[1] = w1 - w2 + w6;
                x[2] = w3 - w4 + w5 - w6;
                
                y[0] = 1-w1 +w4;
                y[1] = w1-w2;
                y[2] = w3 - w4;
                
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1;
                x[1] = w1 - w2 + w6;
                x[2] = w3 - w4;
                
                y[0] = 1-w1 +w5;
                y[1] = w1-w2;
                y[2] = w3 - w5;
                
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1;
                x[1] = w1 - w2;
                x[2] = w3;
                
                y[0] = 1-w1 +w6;
                y[1] = w1-w2 +w5-w6;
                y[2] = w3 - w4;
                
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1;
                x[1] = w1 - w2;
                x[2] = w3 - w4;
                
                y[0] = 1-w1 +w6;
                y[1] = w1-w2 +w5-w6;
                y[2] = w3 - w5;
                
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
                //==================================//
                
                
                x[0] = 1-w1 +w6;
                x[1] = w1-w2;
                x[2] = w3 - w4;
                
                y[0] = 1-w1;
                y[1] = w1 - w2 + w5;
                y[2] = w3 - w5;
                
                
                          
                // Case 4 -> index 4
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1 +w6;
                x[1] = w1-w2;
                x[2] = w3 - w4 + w5 - w6;
                
                y[0] = 1-w1;
                y[1] = w1 - w2 + w4;
                y[2] = w3 - w4;
                
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
	              
                x[0] = 1-w1 +w6;
                x[1] = w1-w2;
                x[2] = w3 - w6;
                
                y[0] = 1-w1;
                y[1] = w1 - w2 + w5;
                y[2] = w3 - w4;
                
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
	              
                x[0] = 1-w1 +w5;
                x[1] = w1-w2;
                x[2] = w3 - w4;
                
                y[0] = 1-w1;
                y[1] = w1 - w2 + w6;
                y[2] = w3 - w6;
                
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1 +w4;
                x[1] = w1-w2;
                x[2] = w3 - w4;
                
                y[0] = 1-w1;
                y[1] = w1 - w2 + w6;
                y[2] = w3 - w4 + w5 - w6;
                
                
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-w1 +w5;
                x[1] = w1-w2;
                x[2] = w3 - w5;
                
                y[0] = 1-w1;
                y[1] = w1 - w2 + w6;
                y[2] = w3 - w4;
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                x[0] = 1-w1 +w6;
                x[1] = w1-w2 +w5-w6;
                x[2] = w3 - w4;
                
                y[0] = 1-w1;
                y[1] = w1 - w2;
                y[2] = w3;
                
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
                
                
	              //---------------------------------------//
	              
                x[0] = 1-w1 +w6;
                x[1] = w1-w2 +w5-w6;
                x[2] = w3 - w5;
                
                y[0] = 1-w1;
                y[1] = w1 - w2;
                y[2] = w3 - w4;
                
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
	              /*
                const Real& w1  = z1;
                const Real& w2  = z1*z2;
                const Real& w3  = z1*z2*z3;
                const Real& w4  = z1*z2*z3*z4;
                const Real& w5  = z1*z2*z3*z4*e1;
                const Real& w6  = z1*z2*z3*z4*e2;
                */
                
                
                w1  = z1;
                w2  = z1*z2;
                w3  = z1*z2*z3;
                w4  = z1*z2*z3*z4;
                w5  = z1*z2*z3*z4*e1;
                w6  = z1*z2*z3*z4*e2;
                
                x[0] = 1-w1;
                x[1] = w1 - w2;
                x[2] = w3 - w4 + w5;
                
                y[0] = 1-w1 +w6;
                y[1] = w1-w2 +w4-w6;
                y[2] = w3 - w4;
                
                w    = z1*z1*z1*z1*z1*z2*z2*z2*z2*z3*z3*z3*z4*z4*dw;
                          
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-w1 +w6;
                x[1] = w1-w2 +w4-w6;
                x[2] = w3 - w4;
                
                y[0] = 1-w1;
                y[1] = w1 - w2;
                y[2] = w3 - w4 + w5;
                
	              x_[4].push_back(x);
	              y_[4].push_back(y);
	              w_[4].push_back(w);
	              
	            }     
            }
          }
        }
      }
    }


    // R3 xtemp, ytemp;
    // R3x3 B; B(0,0)=1.;B(1,1)=1.;B(0,1)=-1.;B(0,2)=-1.;B(2,2)=1.;
    // for(int q=2; q<5; q++){
    //   for(int j=0; j<w_[q].size(); j++){
    //     xtemp = B*x_[q][j]; x_[q][j] = xtemp;
    //     ytemp = B*y_[q][j]; y_[q][j] = ytemp;
    //   }
    // }
    




    }

  };
  
  
  
  
  
  // -----------------------------------------------------------------------------//
  
  template <> class QuadVol<2,3>{

  public:
    typedef R3   qp_t;

  private:
    // arrays of std::vectors
    std::vector<R2>   x_[4]; // 4 cases: disjoint, common vertex, edge, face
    std::vector<R3>   y_[4];
    std::vector<Real> w_[4];

  public:
    const std::vector<R2>&   x(const int& rule) const {return x_[rule];}
    const std::vector<R3>&   y(const int& rule) const {return y_[rule];}
    const std::vector<Real>& w(const int& rule) const {return w_[rule];}

    // Constructor
    QuadVol<2,3>(const int& order){

      R2 x; R3 y; Real w,dw,dww;
      std::vector<Real> t,dt, tt, dtt;
      Quad1D(order,t,dt);
      int nt = t.size();

      
      std::vector<R3>   qpx;
      std::vector<R2>   qpy;
      std::vector<Real> qw;




    //================================//
    //   Cas 0: Disjoint              //
    //================================//
      std::vector<R2>   qp2;
      std::vector<R3>   qp3;
      std::vector<Real> qw2, qw3;

      Quad2D(10,qp2,qw2);
      Quad3D(4,qp3,qw3);
      int nq2 = qp2.size();

      int nq3 = qp3.size();

    for(int j=0; j<nq2; j++){
      for(int k=0; k<nq3; k++){

      x[0] = qp2[j][0];
      x[1] = qp2[j][1];
      y[0] = qp3[k][0];
      y[1] = qp3[k][1];
      y[2] = qp3[k][2];	
      w    = qw2[j]*qw3[k];

      // Case 0 -> index 0
      x_[0].push_back(x);
      y_[0].push_back(y);
      w_[0].push_back(w);

      }
    }
      
    // for(int j0=0; j0<ntt; j0++){
    //   for(int j1=0; j1<ntt; j1++){
	  //     for(int j2=0; j2<ntt; j2++){
	  //       for(int j3=0; j3<ntt; j3++){
	  //         for(int j4=0; j4<ntt; j4++){
	                        
                      
    //           //================================//
    //           //         Case 0: Disjoint       //
    //           //================================//
    //             Real ee1  = tt[j0];
    //             Real ee2  = tt[j1];
    //             Real ee3  = tt[j2];
    //             Real ee4  = tt[j3];
    //             Real ee5  = tt[j4];
                
                
    //       	    dww = dtt[j0]*dtt[j1]*dtt[j2]*dtt[j3]*dtt[j4];


    //             x[0] = ee4*ee5;
    //             x[1] = ee4*(1.-ee5);
                
    //             y[0] = ee1*ee2*ee3;
    //             y[1] = ee1*ee2*(1.-ee3);
    //             y[2] = ee1*(1.-ee2);
                
    //             w    = ee1*ee1*ee2*ee4*dww;

    //             x_[0].push_back(x);
    //             y_[0].push_back(y);
    //             w_[0].push_back(w);
                
                
    //           }
    //         }
    //       }
    //     }
    //   }
      
      

    for(int j0=0; j0<nt; j0++){
      for(int j1=0; j1<nt; j1++){
	      for(int j2=0; j2<nt; j2++){
	        for(int j3=0; j3<nt; j3++){
	          for(int j4=0; j4<nt; j4++){
	                        
                          
              //================================//
              //         Case 0: Disjoint       //
              //================================//
                Real e1  = t[j0];
                Real e2  = t[j1];
                Real e3  = t[j2];
                Real e4  = t[j3];
                Real z   = t[j4];
                
                
          	    dw = dt[j0]*dt[j1]*dt[j2]*dt[j3]*dt[j4];

/*
                x[0] = e1*e2*e3;
                x[1] = e1*e2*(1-e3);
                x[2] = e1*(1-e2);
                
                y[0] = e4*e5*e6;
                y[1] = e4*e5*(1-e6);
                y[2] = e4*(1-e5);
                
                w    = e1*e1*e2*e4*e4*e5*dw;

                x_[0].push_back(x);
                y_[0].push_back(y);
                w_[0].push_back(w);
  */              
                
              //================================//
              //      Case 1: Common Vertex     //
              //================================//
              
              /*
                const Real& e1  = t[j0];
                const Real& e2  = t[j1];
                const Real& e3  = t[j2];
                const Real& e4  = t[j3];
                const Real& e5  = t[j4];
                const Real& z   = t[j5];
                */
                
                
                e1  = t[j0];
                e2  = t[j1];
                e3  = t[j2];
                e4  = t[j3];
                z   = t[j4];
                
                x[0] = z*e3*e4;
                x[1] = z*e3*(1-e4);
                
                y[0] = z*e1*e2;
                y[1] = z*e1*(1-e2);
                y[2] = z*(1-e1);
                
                w    = z*z*z*z*e1*e3*dw;
                          
                // Case 1 -> index 1
	              x_[1].push_back(x);
	              y_[1].push_back(y);
	              w_[1].push_back(w);
	              
	              //---------------------------------------//
                
                x[0] = z*e1;
                x[1] = z*(1-e1);
                
                y[0] = z*e2*e3*e4;
                y[1] = z*e2*e3*(1-e4);
                y[2] = z*e2*(1-e3);
                
                w    = z*z*z*z*e2*e2*e3*dw;
                
	              x_[1].push_back(x);
	              y_[1].push_back(y);
	              w_[1].push_back(w);
	              
	              
              //================================//
              //      Case 2: Common Edge       //
              //================================//
              
              /*
                const Real& e1  = t[j0];
                const Real& e2  = t[j1];
                const Real& e3  = t[j2];
                const Real& e4  = t[j3];
                const Real& z1  = t[j4];
                const Real& z2  = t[j5];
                */
                
                e1       = t[j0];
                e2       = t[j1];
                e3       = t[j2];
                Real z1  = t[j3];
                Real z2  = t[j4];
                
                
                x[0] = 1-z1 + z1*z2*e2*e3; 
                x[1] = z1*z2*e2*(1-e3);
                
                y[0] = 1-z1;
                y[1] = z1*z2*e1;
                y[2] = z1*z2*(1-e1);
                
                w    = z1*z1*z1*z1*z2*z2*z2*e2*dw;
                          
                // Case 2 -> index 2
	              x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1 + z1*z2*e1;
                x[1] = z1*z2*(1-e1);
                
                y[0] = 1-z1;
                y[1] = z1*z2*e2*e3;
                y[2] = z1*z2*e2*(1-e3);
                
                x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1;
                x[1] = z1*z2*e1;
                
                y[0] = 1-z1 + z1*z2*e2*e3;
                y[1] = z1*z2*e2*(1-e3);
                y[2] = z1*z2*(1-e2);
                
                          
                x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
              
              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1;
                x[1] = z1*z2;
                
                y[0] = 1-z1 + z1*z2*e1*e2*e3;
                y[1] = z1*z2*e1*e2*(1-e3);
                y[2] = z1*z2*e1*(1-e2);
                
                w    = z1*z1*z1*z1*z2*z2*z2*e1*e1*e2*dw;
                          
                x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
	              
	              
              //================================//
              //      Case 3: Common Face       //
              //================================//
              
                
                
                e1       = t[j0];
                e2       = t[j1];
                z1       = t[j2];
                z2       = t[j3];
                Real z3  = t[j4];
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2;
                x[1] = z1*(1-z2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1;
                y[2] = z1*z2*z3*(1-e1);
                
                w    = z1*z1*z1*z1*z2*z2*z2*z3*z3*e1*dw;
                          
                // Case 3 -> index 3
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1 + z1*z2*z3;
                x[1] = z1*(1-z2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                y[2] = z1*z2*z3*e1*(1-e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1 + z1*z2*z3*e1;
                x[1] = z1*(1-z2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                y[2] = z1*z2*z3*(1-e1*e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*(1-e2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3;
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2);
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*(1-e2);
                y[2] = z1*z2*z3*(1-e1);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
                
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                
                y[0] = 1-z1 + z1*z2*z3*e1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*(1-e1);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3;
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*e1*(1-e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1;
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*(1-e1*e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                
                x[0] = 1-z1 + z1*z2*z3*e1;
                x[1] = z1*(1-z2)+ z1*z2*z3*(1-e1);
                
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*e2;
                
                w    = z1*z1*z1*z1*z2*z2*z2*z3*z3*dw;
                          
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              
	            
            }
          }
        }
      }
    }
    


    // R3x3 B3; B3(0,0)=1.;B3(1,1)=1.;B3(0,1)=-1.;B3(0,2)=-1.;B3(2,2)=1.;
    // R2x2 B2; B2(0,0)=1.;B2(1,1)=1.;B2(0,1)=-1.;

    // R2 xtemp;
    // R3 ytemp;

    // for(int q=2; q<4; q++){
    //   for(int j=0; j<w_[q].size(); j++){
    //     ytemp = B3*y_[q][j];
    //     y_[q][j] = ytemp;
    //   }
    // }

    // for(int q=2; q<4; q++){
    //   for(int j=0; j<w_[q].size(); j++){

    //     xtemp = B2*x_[q][j];
    //     x_[q][j] = xtemp; 
    //   }
    // }

    }

  };
  
  
  
  
    
  // -----------------------------------------------------------------------------//
  
  template <> class QuadVol<3,2>{

  public:
    typedef R3   qp_t;

  private:
    // arrays of std::vectors
    std::vector<R2>   x_[4]; // 4 cases: disjoint, common vertex, edge, face
    std::vector<R3>   y_[4];
    std::vector<Real> w_[4];

  public:
    const std::vector<R3>&   x(const int& rule) const {return y_[rule];} // HERE!
    const std::vector<R2>&   y(const int& rule) const {return x_[rule];} // HERE!
    const std::vector<Real>& w(const int& rule) const {return w_[rule];}

    // Constructor
    QuadVol<3,2>(const int& order){

      R2 x; R3 y; Real w,dw,dww;
      std::vector<Real> t,dt, tt, dtt;
      Quad1D(order,t,dt);
      int nt = t.size();

      
      std::vector<R3>   qpx;
      std::vector<R2>   qpy;
      std::vector<Real> qw;
      



    //================================//
    //   Cas 0: Disjoint              //
    //================================//
      std::vector<R2>   qp2;
      std::vector<R3>   qp3;
      std::vector<Real> qw2, qw3;

      Quad2D(10,qp2,qw2);
      Quad3D(4,qp3,qw3);
      int nq2 = qp2.size();

      int nq3 = qp3.size();

    for(int j=0; j<nq2; j++){
      for(int k=0; k<nq3; k++){

      x[0] = qp2[j][0];
      x[1] = qp2[j][1];
      y[0] = qp3[k][0];
      y[1] = qp3[k][1];
      y[2] = qp3[k][2];
      w    = qw2[j]*qw3[k];

      // Case 0 -> index 0
      x_[0].push_back(x);
      y_[0].push_back(y);
      w_[0].push_back(w);

      }
    }
    // for(int j0=0; j0<ntt; j0++){
    //   for(int j1=0; j1<ntt; j1++){
	  //     for(int j2=0; j2<ntt; j2++){
	  //       for(int j3=0; j3<ntt; j3++){
	  //         for(int j4=0; j4<ntt; j4++){
	                        
                      
    //           //================================//
    //           //         Case 0: Disjoint       //
    //           //================================//
    //             Real ee1  = tt[j0];
    //             Real ee2  = tt[j1];
    //             Real ee3  = tt[j2];
    //             Real ee4  = tt[j3];
    //             Real ee5  = tt[j4];
                
                
    //       	    dww = dtt[j0]*dtt[j1]*dtt[j2]*dtt[j3]*dtt[j4];


    //             x[0] = ee4*ee5;
    //             x[1] = ee4*(1.-ee5);
                
    //             y[0] = ee1*ee2*ee3;
    //             y[1] = ee1*ee2*(1.-ee3);
    //             y[2] = ee1*(1.-ee2);
                
    //             w    = ee1*ee1*ee2*ee4*dww;

    //             x_[0].push_back(x);
    //             y_[0].push_back(y);
    //             w_[0].push_back(w);
                
                
    //           }
    //         }
    //       }
    //     }
    //   }
      
      

    for(int j0=0; j0<nt; j0++){
      for(int j1=0; j1<nt; j1++){
	      for(int j2=0; j2<nt; j2++){
	        for(int j3=0; j3<nt; j3++){
	          for(int j4=0; j4<nt; j4++){
	                        
                          
              //================================//
              //         Case 0: Disjoint       //
              //================================//
                Real e1  = t[j0];
                Real e2  = t[j1];
                Real e3  = t[j2];
                Real e4  = t[j3];
                Real z   = t[j4];
                
                
          	    dw = dt[j0]*dt[j1]*dt[j2]*dt[j3]*dt[j4];

/*
                x[0] = e1*e2*e3;
                x[1] = e1*e2*(1-e3);
                x[2] = e1*(1-e2);
                
                y[0] = e4*e5*e6;
                y[1] = e4*e5*(1-e6);
                y[2] = e4*(1-e5);
                
                w    = e1*e1*e2*e4*e4*e5*dw;

                x_[0].push_back(x);
                y_[0].push_back(y);
                w_[0].push_back(w);
  */              
                
              //================================//
              //      Case 1: Common Vertex     //
              //================================//
              
              /*
                const Real& e1  = t[j0];
                const Real& e2  = t[j1];
                const Real& e3  = t[j2];
                const Real& e4  = t[j3];
                const Real& e5  = t[j4];
                const Real& z   = t[j5];
                */
                
                
                e1  = t[j0];
                e2  = t[j1];
                e3  = t[j2];
                e4  = t[j3];
                z   = t[j4];
                
                x[0] = z*e3*e4;
                x[1] = z*e3*(1-e4);
                
                y[0] = z*e1*e2;
                y[1] = z*e1*(1-e2);
                y[2] = z*(1-e1);
                
                w    = z*z*z*z*e1*e3*dw;
                          
                // Case 1 -> index 1
	              x_[1].push_back(x);
	              y_[1].push_back(y);
	              w_[1].push_back(w);
	              
	              //---------------------------------------//
                
                x[0] = z*e1;
                x[1] = z*(1-e1);
                
                y[0] = z*e2*e3*e4;
                y[1] = z*e2*e3*(1-e4);
                y[2] = z*e2*(1-e3);
                
                w    = z*z*z*z*e2*e2*e3*dw;
                
	              x_[1].push_back(x);
	              y_[1].push_back(y);
	              w_[1].push_back(w);
	              
	              
              //================================//
              //      Case 2: Common Edge       //
              //================================//
              
              /*
                const Real& e1  = t[j0];
                const Real& e2  = t[j1];
                const Real& e3  = t[j2];
                const Real& e4  = t[j3];
                const Real& z1  = t[j4];
                const Real& z2  = t[j5];
                */
                
                e1       = t[j0];
                e2       = t[j1];
                e3       = t[j2];
                Real z1  = t[j3];
                Real z2  = t[j4];
                
                
                x[0] = 1-z1 + z1*z2*e2*e3; 
                x[1] = z1*z2*e2*(1-e3);
                
                y[0] = 1-z1;
                y[1] = z1*z2*e1;
                y[2] = z1*z2*(1-e1);
                
                w    = z1*z1*z1*z1*z2*z2*z2*e2*dw;
                          
                // Case 2 -> index 2
	              x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1 + z1*z2*e1;
                x[1] = z1*z2*(1-e1);
                
                y[0] = 1-z1;
                y[1] = z1*z2*e2*e3;
                y[2] = z1*z2*e2*(1-e3);
                
                x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1;
                x[1] = z1*z2*e1;
                
                y[0] = 1-z1 + z1*z2*e2*e3;
                y[1] = z1*z2*e2*(1-e3);
                y[2] = z1*z2*(1-e2);
                
                          
                x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
              
              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1;
                x[1] = z1*z2;
                
                y[0] = 1-z1 + z1*z2*e1*e2*e3;
                y[1] = z1*z2*e1*e2*(1-e3);
                y[2] = z1*z2*e1*(1-e2);
                
                w    = z1*z1*z1*z1*z2*z2*z2*e1*e1*e2*dw;
                          
                x_[2].push_back(x);
	              y_[2].push_back(y);
	              w_[2].push_back(w);
	              
	              
              //================================//
              //      Case 3: Common Face       //
              //================================//
              
                
                
                e1       = t[j0];
                e2       = t[j1];
                z1       = t[j2];
                z2       = t[j3];
                Real z3  = t[j4];
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2;
                x[1] = z1*(1-z2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1;
                y[2] = z1*z2*z3*(1-e1);
                
                w    = z1*z1*z1*z1*z2*z2*z2*z3*z3*e1*dw;
                          
                // Case 3 -> index 3
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1 + z1*z2*z3;
                x[1] = z1*(1-z2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                y[2] = z1*z2*z3*e1*(1-e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1 + z1*z2*z3*e1;
                x[1] = z1*(1-z2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                y[2] = z1*z2*z3*(1-e1*e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1 + z1*z2*z3*e1*e2;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*(1-e2);
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3;
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2);
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2;
                y[1] = z1*(1-z2) + z1*z2*z3*e1*(1-e2);
                y[2] = z1*z2*z3*(1-e1);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
                
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1*e2;
                
                y[0] = 1-z1 + z1*z2*z3*e1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*(1-e1);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              //---------------------------------------//
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3;
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*e1*(1-e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
	              
                
                x[0] = 1-z1;
                x[1] = z1*(1-z2) + z1*z2*z3*e1;
                
                y[0] = 1-z1 + z1*z2*z3*e1*e2;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*(1-e1*e2);
                
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              //---------------------------------------//
                
                
                x[0] = 1-z1 + z1*z2*z3*e1;
                x[1] = z1*(1-z2)+ z1*z2*z3*(1-e1);
                
                
                y[0] = 1-z1;
                y[1] = z1*(1-z2);
                y[2] = z1*z2*z3*e2;
                
                w    = z1*z1*z1*z1*z2*z2*z2*z3*z3*dw;
                          
	              x_[3].push_back(x);
	              y_[3].push_back(y);
	              w_[3].push_back(w);
	              
	              
	              
	            
            }
          }
        }
      }
    }
    


    R3x3 B3; B3(0,0)=1.;B3(1,1)=1.;B3(0,1)=-1.;B3(0,2)=-1.;B3(2,2)=1.;
    R2x2 B2; B2(0,0)=1.;B2(1,1)=1.;B2(0,1)=-1.;


    // R2 xtemp;
    // R3 ytemp;
    // for(int q=2; q<4; q++){
    //   for(int j=0; j<w_[q].size(); j++){
    //     ytemp = B3*y_[q][j];
    //     y_[q][j] = ytemp;
    //   }
    // }

    // for(int q=2; q<4; q++){
    //   for(int j=0; j<w_[q].size(); j++){

    //     xtemp = B2*x_[q][j];
    //     x_[q][j] = xtemp;
    //   }
    // }


    }

  };







}







#endif
