//===================================================================
//
//  Copyright 2017 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef BEMTOOL_OPERATOR_MAXWELL_HPP
#define BEMTOOL_OPERATOR_MAXWELL_HPP

#include "operator.hpp"

namespace bemtool {


/*====
  EFIE
  ====*/

template <typename PhiX, typename PhiY>
class BIOpKernel<MA,SL_OP,3,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<MA,SL_OP,3,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&      meshx;
  const typename Trait::MeshY&      meshy;
        typename Trait::MatType     inter;
        typename Trait::JacX        dx;
        typename Trait::JacY        dy;
        typename Trait::DivPhiX     div_phix;
        typename Trait::DivPhiY     div_phiy;
                        PhiX        phix;
                        PhiY        phiy;
  const                 Real        kappa, inv_kappa2;
                        R3          x0_y0,x_y,nx,ny,y,x,y0,x0;
                        Real        h,r;
                        Cplx        ker,val,val2;
                        std::function<Real(R3)> f1, f2;

public:
  BIOpKernel<MA,SL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), div_phix(mx), meshy(my), phiy(my), div_phiy(my),
    kappa(k), inv_kappa2( 1./(kappa*kappa) ), f1([](R3 x) {return 1.0;}), 
    f2([](R3 x) {return 1.0;}) {};


  BIOpKernel<MA,SL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k, std::function<Real(R3)> func):
  meshx(mx), phix(mx), div_phix(mx), meshy(my), phiy(my), div_phiy(my),
   kappa(k), inv_kappa2( 1./(kappa*kappa) ), f1(func), f2([](R3 x) {return 1.0;}) {};
   
  BIOpKernel<MA,SL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k, std::function<Real(R3)> func1, std::function<Real(R3)> func2):
  meshx(mx), phix(mx), div_phix(mx), meshy(my), phiy(my), div_phiy(my),
   kappa(k), inv_kappa2( 1./(kappa*kappa) ), f1(func1), f2(func2) {};

  inline void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    phix.Assign(ix);
    phiy.Assign(iy);
    div_phix.Assign(ix);
    div_phiy.Assign(iy);
    h     = DetJac(ex)*DetJac(ey);
    x0_y0 = ex[0]-ey[0];
    x0    = ex[0];
    y0    = ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
  }


  inline const typename Trait::MatType&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);

    if( r < 1e-8){ std::cout << "ALERT! Small number: " << r << std::endl;}
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	val = ( phix(j,tx),phiy(k,ty) )*f2(y)*ker;
	inter(j,k) = val - inv_kappa2*div_phix(j,tx)*f1(y)*div_phiy(k,ty)*ker;

    if(std::isnan(std::abs(inter(j,k)))){std::cout << "WOW! It is not a number!" << std::endl;}
      }
    }
    return inter;
  }


  inline const Cplx&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty,
	     const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    val = ( phix(kx,tx),phiy(ky,ty) )*f2(y)*ker;

    if(std::isnan(std::abs(val))){std::cout << "WOW! It is not a number!" << std::endl;}
    return val2 = val - inv_kappa2*div_phix(kx,tx)*f1(y)*div_phiy(ky,ty)*ker;
  }


};


typedef BIOpKernel<MA,SL_OP,3,RT0_2D,RT0_2D> EFIE_RT0xRT0;


/*====
  MFIE
  ====*/

template <typename PhiX, typename PhiY>
class BIOpKernel<MA,DL_OP,3,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<MA,DL_OP,3,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&      meshx;
  const typename Trait::MeshY&      meshy;
        typename Trait::MatType     inter;
        typename Trait::JacX        dx;
        typename Trait::JacY        dy;
  const          std::vector<R3>&   normaly;
                        PhiX        phix;
                        PhiY        phiy;
  const                 Real        kappa;
                        R3          x0_y0,x_y,nx,ny, y, x,y0,x0;
                        Real        h,r,r3;
                        C3          ker;
                        Cplx        val;

public:
  BIOpKernel<MA,DL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), meshy(my), normaly(NormalTo(my)), phiy(my), kappa(k){};


  inline void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    phix.Assign(ix);
    phiy.Assign(iy);
    ny    = normaly[iy];
    h     = DetJac(ex)*DetJac(ey);
    x0_y0 = ex[0]-ey[0];
    x0    = ex[0];
    y0    = ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
  }


  inline const typename Trait::MatType&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);

    r3  = r*r*r;
    ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = (ker,vprod(phiy(k,ty),phix(j,tx)));

    if(std::isnan(std::abs(inter(j,k)))){std::cout << "WOW! It is not a number!" << std::endl;}
      }
    }
    return inter;
  }


  inline const Cplx&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty,
	     const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);

    r3  = r*r*r;
    ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;
    val = (ker,vprod(phiy(ky,ty),phix(kx,tx)));
    if(std::isnan(std::abs(val))){std::cout << "WOW! It is not a number!" << std::endl;}
    return val;
  }





  inline const typename Trait::MatType&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty, const std::string& crossprod){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);

    r3  = r*r*r;
    ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = (ker,vprod(vprod(phiy(k,ty), ny),phix(j,tx)));

    if(std::isnan(std::abs(inter(j,k)))){std::cout << "WOW! It is not a number!" << std::endl;}
      }
    }
    return inter;
  }


  inline const Cplx&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty,
	     const int& kx, const int& ky, const std::string& crossprod){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);

    r3  = r*r*r;
    ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;
    val = (ker,vprod(vprod(phiy(ky,ty), ny),phix(kx,tx)));
    if(std::isnan(std::abs(val))){std::cout << "WOW! It is not a number!" << std::endl;}
    return val;
  }


};


typedef BIOpKernel<MA,DL_OP,3,RT0_2D,RT0_2D> MFIE_RT0xRT0;






















template <typename PhiX, typename PhiY>
class VIOpKernel<MA,VL_OP,3,PhiX,PhiY>{

public:
  typedef VIOpKernelTraits<MA,VL_OP,3,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&   meshx;
  const typename Trait::MeshY&   meshy;
        typename Trait::MatType  inter;
        typename Trait::JacX     dx;
        typename Trait::JacY     dy;
                        PhiX     phix;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y, y, x,y0,x0;
                        Real     h,r;
                        Cplx     ker;
                        std::function<Real(R3)> f;

public:
  VIOpKernel<MA,VL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k), f([](R3 x) {return 1.0;}) {};


  VIOpKernel<MA,VL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k, std::function<Real(R3)> func):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k), f(func) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    phix.Assign(ix);
    phiy.Assign(iy);
    x0_y0 = ex[0]-ey[0];
    x0    = ex[0];
    y0    = ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*(phix(j,tx), f(y) * phiy(k,ty));

    if(std::isnan(std::abs(inter(j,k)))){
      
      
      std::cout << "WOW! It is not a number!" << std::endl;
      std::cout << "Value of r:" << r << std::endl;
      std::cout << "Value of f:" << f(y) << std::endl;

      std::cout << "Value of Phi(x):" << phix(j,tx) << std::endl;
      std::cout << "Value of Phi(y):" << phiy(k,ty) << std::endl;
      
      }
      }
    }

    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    ker *= (phix(kx,tx), f(y) * phiy(ky,ty));


    if(std::isnan(std::abs(ker))){
      
      
      std::cout << "WOW! It is not a number!" << std::endl;
      std::cout << "Value of r:" << r << std::endl;
      std::cout << "Value of f:" << f(y) << std::endl;

      std::cout << "Value of Phi(x):" << phix(kx,tx) << std::endl;
      std::cout << "Value of Phi(y):" << phiy(ky,ty) << std::endl;
      
      }


    return ker;
  }


};


typedef VIOpKernel<MA,VL_OP,3,NED_3D,NED_3D> MA_VL_NEDxNED;
typedef VIOpKernel<MA,VL_OP,3,NED_3D,Curl_NED_3D> MA_VL_NEDxCurlNED;
typedef VIOpKernel<MA,VL_OP,3,RT0_2D,Curl_NED_3D> MA_VL_RT0xCurlNED;
typedef VIOpKernel<MA,VL_OP,3,RT0_2D,NED_3D> MA_VL_RT0xNED;






template <typename PhiX, typename PhiY>
class VIOpKernel<MA,DVL_OP,3,PhiX,PhiY>{

public:
  typedef VIOpKernelTraits<MA,DVL_OP,3,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&   meshx;
  const typename Trait::MeshY&   meshy;
        typename Trait::MatType  inter;
        typename Trait::JacX     dx;
        typename Trait::JacY     dy;
                        PhiX     phix;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y, y, x,y0,x0;
                        Real     h,r;
                        Cplx     ker;

public:
  VIOpKernel<MA,DVL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    phix.Assign(ix);
    phiy.Assign(iy);
    x0_y0 = ex[0]-ey[0];
    x0    = ex[0];
    y0    = ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*(phix(j,tx), phiy(k,ty));
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    return ker *= (phix(kx,tx), phiy(ky,ty));
  }


};


typedef VIOpKernel<MA,DVL_OP,3,NED_3D,Curl_NED_3D> MA_DVL_3D_NEDxNED;




/*====
  curl N(p curlE)
  ====*/


template <typename PhiX, typename PhiY>
class VIOpKernel<MA,KVL_OP,3,PhiX,PhiY>{

public:
  typedef VIOpKernelTraits<MA,KVL_OP,3,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&   meshx;
  const typename Trait::MeshY&   meshy;
        typename Trait::MatType  inter;
        typename Trait::JacX     dx;
        typename Trait::JacY     dy;
                        PhiX     phix;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y, y, x,y0,x0;
                        Real     h,r,r3;
                        C3       ker;
                        Cplx     val;
                        std::function<Real(R3)> f;

public:

  VIOpKernel<MA,KVL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k), f([](R3 x){return 1.0;}) {};


  VIOpKernel<MA,KVL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k, std::function<Real(R3)> func):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k), f(func) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    phix.Assign(ix);
    phiy.Assign(iy);
    x0_y0 = ex[0]-ey[0];
    x0    = ex[0];
    y0    = ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    r3  = r*r*r;
    ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;
    
    
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = (ker, vprod(phix(j,tx), f(y) * phiy(k,ty)));
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
	  x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    r3  = r*r*r;
    y   = y0 + dy*ty;
    
    
    ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;
    
    
    val = (ker, vprod(phix(kx,tx), f(y) * phiy(ky,ty)));
    
    
    return val;
  }


};


typedef VIOpKernel<MA,KVL_OP,3,NED_3D,Curl_NED_3D> MA_KVL_3D_NEDxNED;

typedef VIOpKernel<MA,KVL_OP,3,Curl_NED_3D,NED_3D> MA_KVL_CurlNEDxNED;

typedef VIOpKernel<MA,KVL_OP,3,NED_3D,RT0_2D> MA_DL_NEDxRT0;

typedef VIOpKernel<MA,KVL_OP,3,RT0_2D,NED_3D> MA_KVL_RT0xNED;


// /*====
//   curl S  ---- MAXWELL DOUBLE LAYER 
//   ====*/


// template <typename PhiX, typename PhiY>
// class VIOpKernel<MA,KDL_OP,3,PhiX,PhiY>{

// public:
//   typedef VIOpKernelTraits<MA,KVL_OP,3,PhiX,PhiY> Trait;

// private:
//   const typename Trait::MeshX&   meshx;
//   const typename Trait::MeshY&   meshy;
//         typename Trait::MatType  inter;
//         typename Trait::JacX     dx;
//         typename Trait::JacY     dy;
//                         PhiX     phix;
//                         PhiY     phiy;
//   const                 Real     kappa;
//                         R3       x0_y0,x_y, y, x,y0,x0;
//                         Real     h,r,r3;
//                         C3       ker;
//                         Cplx     val;
//                         std::function<Real(R3)> f;

// public:

//   VIOpKernel<MA,KVL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
// 				   const typename Trait::MeshY& my,
// 				   const Real& k):
//   meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k), f([](R3 x){return 1.0;}) {};


//   VIOpKernel<MA,KVL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
// 				   const typename Trait::MeshY& my,
// 				   const Real& k, std::function<Real(R3)> func):
//   meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k), f(func) {};

//   void Assign(const int& ix, const int& iy){
//     const typename Trait::EltX& ex=meshx[ix];
//     const typename Trait::EltY& ey=meshy[iy];
//     phix.Assign(ix);
//     phiy.Assign(iy);
//     x0_y0 = ex[0]-ey[0];
//     x0    = ex[0];
//     y0    = ey[0];
//     dx    = MatJac(ex);
//     dy    = MatJac(ey);
//     h     = DetJac(ex)*DetJac(ey);
//   }

//   const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
// 					    const typename Trait::Rdy& ty){
//     x_y = x0_y0 + dx*tx-dy*ty;
//     x   = x0 + dx*tx;
//     y   = y0 + dy*ty;
//     r   = norm2(x_y);
//     r3  = r*r*r;
//     ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;
    
    
//     for(int j=0; j<Trait::nb_dof_x; j++){
//       for(int k=0; k<Trait::nb_dof_y; k++){
// 	inter(j,k) = (ker, vprod(phix(j,tx), f(y) * phiy(k,ty)));
//       }
//     }
//     return inter;
//   }

//   const Cplx& operator()(const typename Trait::Rdx& tx,
// 			 const typename Trait::Rdy& ty,
// 			 const int& kx, const int& ky){
// 	  x_y = x0_y0 + dx*tx-dy*ty;
//     r   = norm2(x_y);
//     r3  = r*r*r;
//     y   = y0 + dy*ty;
    
    
//     ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;
    
    
//     val = (ker, vprod(phix(kx,tx), f(y) * phiy(ky,ty)));
    
    
//     return val;
//   }


// };


// typedef VIOpKernel<MA,KDL_OP,3,NED_3D,RT0_2D> MA_DL_NEDxRT0;





/*====
  Grad N(tau dot E)
  ====*/

template <typename PhiX, typename PhiY>
class VIOpKernel<MA,GVL_OP,3,PhiX,PhiY>{

public:
  typedef VIOpKernelTraits<MA,GVL_OP,3,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&   meshx;
  const typename Trait::MeshY&   meshy;
        typename Trait::MatType  inter;
        typename Trait::JacX     dx;
        typename Trait::JacY     dy;
                        PhiX     phix;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y, y, x,y0,x0;
                        Real     h,r,r3;
                        C3       ker;
                        Cplx     val;
                        std::function<R3(R3)> F;

public:
  VIOpKernel<MA,GVL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k, std::function<R3(R3)> func):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k), F(func) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    phix.Assign(ix);
    phiy.Assign(iy);
    x0_y0 = ex[0]-ey[0];
    x0    = ex[0];
    y0    = ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    r3  = r*r*r;
    ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;
    
    
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = (ker, phix(j,tx)) * (F(y), phiy(k,ty));

    if(std::isnan(std::abs(inter(j,k)))){
      
      
      std::cout << "WOW! It is not a number!" << std::endl;
      std::cout << "Value of r:" << r << std::endl;
      std::cout << "Value of f:" << F(y) << std::endl;

      std::cout << "Value of Phi(x):" << phix(j,tx) << std::endl;
      std::cout << "Value of Phi(y):" << phiy(k,ty) << std::endl;
      
      }



      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    
    
    x_y = x0_y0 + dx*tx-dy*ty;
    r   = norm2(x_y);
    r3  = r*r*r;

    y   = y0 + dy*ty;
    ker = h*(iu*kappa*r-1.)*exp(iu*kappa*r)/(4*pi*r3)*x_y;

    if(std::isnan(std::abs(norm2(ker)))){
      
      
      std::cout << "WOW! It is not a number!" << std::endl;
      std::cout << "Value of r:" << r << std::endl;
      std::cout << "Value of F:" << F(y) << std::endl;

      std::cout << "Value of Phi(x):" << phix(kx,tx) << std::endl;
      std::cout << "Value of Phi(y):" << phiy(ky,ty) << std::endl;
      
      }

    val = (ker, phix(kx,tx)) * (F(y), phiy(ky,ty));  

    return val;
  }


};


typedef VIOpKernel<MA,GVL_OP,3,NED_3D,NED_3D> MA_GVL_NEDxNED;
typedef VIOpKernel<MA,GVL_OP,3,RT0_2D,NED_3D> MA_GVL_RT0xNED;


/*====
  N(grad p x E)
  ====*/


template <typename PhiX, typename PhiY>
class VIOpKernel<MA,TVL_OP,3,PhiX,PhiY>{

public:
  typedef VIOpKernelTraits<MA,TVL_OP,3,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&   meshx;
  const typename Trait::MeshY&   meshy;
        typename Trait::MatType  inter;
        typename Trait::JacX     dx;
        typename Trait::JacY     dy;
                        PhiX     phix;
                        PhiY     phiy;
  const                 Real     kappa;
                        R3       x0_y0,x_y, y, x,y0,x0;
                        Real     h,r;
                        Cplx     ker;
                        std::function<R3(R3)> F;

public:
  VIOpKernel<MA,TVL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k, std::function<R3(R3)> func):
  meshx(mx), phix(mx), meshy(my), phiy(my), kappa(k), F(func) {};

  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    phix.Assign(ix);
    phiy.Assign(iy);
    x0_y0 = ex[0]-ey[0];
    x0    = ex[0];
    y0    = ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
    h     = DetJac(ex)*DetJac(ey);
  }

  const typename Trait::MatType& operator()(const typename Trait::Rdx& tx,
					    const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*(phix(j,tx), vprod(F(y), phiy(k,ty)));
      }
    }
    return inter;
  }

  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    return ker *= (phix(kx,tx), vprod(F(y), phiy(ky,ty)));
  }


};


typedef VIOpKernel<MA,TVL_OP,3,NED_3D,NED_3D> MA_TVL_NEDxNED;
typedef VIOpKernel<MA,TVL_OP,3,RT0_2D,NED_3D> MA_TVL_RT0xNED;





/*====
  MAXWELL SINGLE LAYER POTENTIAL 3D
  ====*/

template <typename PhiX, typename PhiY>
class VIOpKernel<MA,SL_OP,3,PhiX,PhiY>{

public:
  typedef VIOpKernelTraits<MA,SL_OP,3,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&      meshx;
  const typename Trait::MeshY&      meshy;
        typename Trait::MatType     inter;
        typename Trait::JacX        dx;
        typename Trait::JacY        dy;
        typename Trait::DivPhiY     div_phiy;
                        PhiX        phix;
                        PhiY        phiy;
  const                 Real        kappa, inv_kappa2;
                        R3          x0_y0,x_y,nx,ny,y,x,y0,x0;
                        Real        h,r;
                        Cplx        ker,ker1,val,val2;
                        std::function<Real(R3)> f1, f2;

public:
  VIOpKernel<MA,SL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k):
  meshx(mx), phix(mx), meshy(my), phiy(my), div_phiy(my),
    kappa(k), inv_kappa2( 1./(kappa*kappa) ), f1([](R3 x) {return 1.0;}), 
    f2([](R3 x) {return 1.0;}) {};
    
    
    
  VIOpKernel<MA,SL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k, std::function<Real(R3)> func):
  meshx(mx), phix(mx), meshy(my), phiy(my), div_phiy(my),
   kappa(k), inv_kappa2( 1./(kappa*kappa) ), f1(func), f2([](R3 x) {return 1.0;}) {};
   
  VIOpKernel<MA,SL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				   const typename Trait::MeshY& my,
				   const Real& k, std::function<Real(R3)> func1, std::function<Real(R3)> func2):
  meshx(mx), phix(mx), meshy(my), phiy(my), div_phiy(my),
   kappa(k), inv_kappa2( 1./(kappa*kappa) ), f1(func1), f2(func2) {};


  inline void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    phix.Assign(ix);
    phiy.Assign(iy);
    div_phiy.Assign(iy);
    h     = DetJac(ex)*DetJac(ey);
    x0_y0 = ex[0]-ey[0];
    x0    = ex[0];
    y0    = ey[0];
    dx    = MatJac(ex);
    dy    = MatJac(ey);
  }


  inline const typename Trait::MatType&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    ker1 = ker*(iu*kappa*r-1)*inv_kappa2/(r*r);
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	val = (phix(j,tx),phiy(k,ty) )*f2(y)*ker;
	inter(j,k) = val + (phix(j,tx), x_y)*ker1*f1(y)*div_phiy(k,ty);
      }
    }
    return inter;
  }


  inline const Cplx&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty,
	     const int& kx, const int& ky){
    x_y = x0_y0 + dx*tx-dy*ty;
    x   = x0 + dx*tx;
    y   = y0 + dy*ty;
    r   = norm2(x_y);
    ker = h*exp(iu*kappa*r)/(4*pi*r);
    ker1 = ker*(iu*kappa*r-1)*inv_kappa2/(r*r);
    val = ( phix(kx,tx),phiy(ky,ty) )*f2(y)*ker;
    return val2 = val + (phix(kx,tx),x_y)*f1(y)*div_phiy(ky,ty)*ker1;
  }


};


typedef VIOpKernel<MA,SL_OP,3,NED_3D,RT0_2D> MA_SL_NEDxRT0;















}
















#endif
