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
#ifndef BEMTOOL_OPERATOR_OPERATOR_HPP
#define BEMTOOL_OPERATOR_OPERATOR_HPP

#include <fstream>
#include "../quadrature/quad_bem.hpp"
#include "../quadrature/quad_vol.hpp"
#include "../fem/shapefct.hpp"
#include "../equations.hpp"
#include "../miscellaneous/specialfct.hpp"



namespace bemtool {

/*=====================
  Type d'operateurs
  =====================*/

enum BIOpKernelEnum {
  CST_OP,  // noyau constant
   SL_OP,  // trace dirichlet du simple couche
   DL_OP,  // trace dirichlet du double couche
  TDL_OP,  // trace neumann   du simple couche
   HS_OP,   // trace neumann   du double couche
   VL_OP,    // Volume 
   DVL_OP,   
   NVL_OP,   
   KVL_OP,    // curl N (curl E)
   GVL_OP,    // Grad N (tau dot E)
   TVL_OP     // N(grad p x E)
};

template <typename KernelType> class BIOp{

public:

  typedef typename KernelType::Trait         KernelTypeTrait;
  typedef typename KernelTypeTrait::MeshX    MeshX;
  typedef typename KernelTypeTrait::MeshY    MeshY;
  typedef typename KernelTypeTrait::Rdx      RdX;
  typedef typename KernelTypeTrait::Rdy      RdY;
  typedef typename KernelTypeTrait::EltX     EltX;
  typedef typename KernelTypeTrait::EltY     EltY;
  typedef typename KernelTypeTrait::MatType  MatType;

  static  const int dimx = KernelTypeTrait::dimx;
  static  const int dimy = KernelTypeTrait::dimy;
  typedef mat<dimx,dimx,Real>                JacX;
  typedef mat<dimy,dimy,Real>                JacY;
  typedef QuadBEM<dimx,dimy>                 QuadType;

private:

  const MeshX&    meshx;
  const MeshY&    meshy;
  const Geometry& nodex;
  const Geometry& nodey;

  KernelType   ker;
  QuadType     qr;
  MatType      inter;
  int          rule;
  JacX         dx;
  JacY         dy;
  RdX          x0,x,ax[dimx+1];
  RdY          y0,y,ay[dimy+1];
  Cplx         val,val2;

public:
  BIOp<KernelType>(const MeshX& mx, const MeshY& my, const Real& k):
  meshx(mx), nodex(GeometryOf(mx)), ker(mx,my,k),
    meshy(my), nodey(GeometryOf(my)), qr(QR) {
    for(int j=0; j<dimx; j++){ax[j+1][j]=1;}
    for(int j=0; j<dimy; j++){ay[j+1][j]=1;}
  };

  
  template <typename function_t>
  BIOp<KernelType>(const MeshX& mx, const MeshY& my, const Real& k, function_t func):
  meshx(mx), nodex(GeometryOf(mx)), ker(mx,my,k,func),
    meshy(my), nodey(GeometryOf(my)), qr(QR) {
    for(int j=0; j<dimx; j++){ax[j+1][j]=1;}
    for(int j=0; j<dimy; j++){ay[j+1][j]=1;}
  };
  
  template <typename function_t>
  BIOp<KernelType>(const MeshX& mx, const MeshY& my, const Real& k, function_t func1,function_t func2):
  meshx(mx), nodex(GeometryOf(mx)), ker(mx,my,k,func1, func2),
    meshy(my), nodey(GeometryOf(my)), qr(QR) {
    for(int j=0; j<dimx; j++){ax[j+1][j]=1;}
    for(int j=0; j<dimy; j++){ay[j+1][j]=1;}
  };

  void ChooseQuad(const int& jx, const int& jy){
    rule = 0; inter=0.; val=0.;
    EltX ex = meshx[jx];   EltY ey = meshy[jy];
    Elt<dimx,RdX> bx(ax);  Elt<dimy,RdY> by(ay);

    for(int p=0; p<dimx+1; p++){
      for(int q=rule; q<dimy+1; q++){
	if( &ex[p]==&ey[q] ){
	  Swap(ex,rule,p); Swap(bx,rule,p);
	  Swap(ey,rule,q); Swap(by,rule,q);
	  rule++; break;
	}
      }
    }
    dx = MatJac(bx); x0=bx[0];
    dy = MatJac(by); y0=by[0];
  }


  const MatType& operator()(const EltX& ex, const EltY& ey){
    return (*this)(nodex[ex][meshx],nodey[ey][meshy]);}


  const MatType& operator()(const int& jx, const int& jy){
    ChooseQuad(jx,jy);
    if(rule != 0){
    const std::vector<RdX>&  s = qr.x(rule);
    const std::vector<RdY>&  t = qr.y(rule);
    const std::vector<Real>& w = qr.w(rule);
    ker.Assign(jx,jy);
    for(int j=0; j<w.size(); j++){
      x = x0 + dx*s[j];
      y = y0 + dy*t[j];
      inter += w[j]*ker(x,y);
    }
    }
    else{
    inter = 0. * inter;
    }
    return inter;


  }
  

  const Cplx& operator()(const N2& Ix, const N2& Iy){
    ChooseQuad(Ix[0],Iy[0]);
    if(rule == 0){
    const std::vector<RdX>&  s = qr.x(rule);
    const std::vector<RdY>&  t = qr.y(rule);
    const std::vector<Real>& w = qr.w(rule);
    ker.Assign(Ix[0],Iy[0]);
    for(int j=0; j<w.size(); j++){
      x = x0 + dx*s[j]; y = y0 + dy*t[j];
      val += w[j]*ker(x,y,Ix[1],Iy[1]);}
    }
    else{
      val = 0.;
    }
    return val;
  }


  const Cplx& operator()(const std::vector<N2>& vx,
			 const std::vector<N2>& vy){
    val2=0.;
    for(int ix=0; ix<vx.size(); ix++){
      for(int iy=0; iy<vy.size(); iy++){
	val2 += (*this)(vx[ix],vy[iy]);} }
    return val2;
  }












  const MatType& operator()(const EltX& ex, const EltY& ey, const std::string& crossprod){
    return (*this)(nodex[ex][meshx],nodey[ey][meshy], crossprod);}


  const MatType& operator()(const int& jx, const int& jy, const std::string& crossprod){
    ChooseQuad(jx,jy);
    if(rule != 0){
    const std::vector<RdX>&  s = qr.x(rule);
    const std::vector<RdY>&  t = qr.y(rule);
    const std::vector<Real>& w = qr.w(rule);
    ker.Assign(jx,jy);
    for(int j=0; j<w.size(); j++){
      x = x0 + dx*s[j];
      y = y0 + dy*t[j];
      inter += w[j]*ker(x,y,crossprod);
    }
    }
    else{
    inter = 0. * inter;
    }
    return inter;


  }
  

  const Cplx& operator()(const N2& Ix, const N2& Iy,const std::string& crossprod){
    ChooseQuad(Ix[0],Iy[0]);
    if(rule == 0){
    const std::vector<RdX>&  s = qr.x(rule);
    const std::vector<RdY>&  t = qr.y(rule);
    const std::vector<Real>& w = qr.w(rule);
    ker.Assign(Ix[0],Iy[0]);
    for(int j=0; j<w.size(); j++){
      x = x0 + dx*s[j]; y = y0 + dy*t[j];
      val += w[j]*ker(x,y,Ix[1],Iy[1],crossprod);}
    }
    else{
      val = 0.;
    }
    return val;
  }


  const Cplx& operator()(const std::vector<N2>& vx,
			 const std::vector<N2>& vy, const std::string& crossprod){
    val2=0.;
    for(int ix=0; ix<vx.size(); ix++){
      for(int iy=0; iy<vy.size(); iy++){
	val2 += (*this)(vx[ix],vy[iy],crossprod);} }
    return val2;
  }


};


/*==============================
  GENERALITES SUR LES NOYAUX
  =============================*/

template <int,int,int,typename,typename>
class BIOpKernel;

template <int EqT, int OpT, int dim, typename PhiX, typename PhiY>
class BIOpKernelTraits{

public:
  static const int EqType   = EqT;
  static const int OpType   = OpT;
  static const int Dim      = dim;
  static const int dimx     = PhiX::dim;
  static const int dimy     = PhiY::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;

  typedef typename Jac<dimx>::Type     JacX;
  typedef typename Jac<dimy>::Type     JacY;
  typedef typename PhiX::Rd            Rdx;
  typedef typename PhiY::Rd            Rdy;
  typedef Mesh<dimx>                   MeshX;
  typedef Mesh<dimy>                   MeshY;
  typedef Elt<dimx>                    EltX;
  typedef Elt<dimy>                    EltY;
  typedef GradBasisFct<PhiX>           GradPhiX;
  typedef GradBasisFct<PhiY>           GradPhiY;
  typedef DivBasisFct<PhiX>            DivPhiX;
  typedef DivBasisFct<PhiY>            DivPhiY;
  typedef mat<nb_dof_x,nb_dof_y,Cplx>  MatType;

  typedef PhiX ShapeFctX;
  typedef PhiY ShapeFctY;  
  
};


/*=======================
        NOYAU CONSTANT
  =======================*/

template <int D, typename PhiX, typename PhiY>
class BIOpKernel<CST,CST_OP,D,PhiX,PhiY>{

public:
  typedef BIOpKernelTraits<CST,CST_OP,D,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&    meshx;
  const typename Trait::MeshY&    meshy;
        typename Trait::MatType   inter;
                        PhiX      phix;
                        PhiY      phiy;
                        Real      ker;
                        Cplx      val;

public:
  BIOpKernel<CST,CST_OP,D,PhiX,PhiY>(const typename Trait::MeshX& mx,
				 const typename Trait::MeshY& my,
				 const Real& k=0.):
  meshx(mx), phix(mx), meshy(my), phiy(my) {};


  void Assign(const int& ix, const int& iy){
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    ker = DetJac(ex)*DetJac(ey);
  }


  const typename Trait::MatType&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty){
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
	inter(j,k) = ker*phix(j,tx)*phiy(k,ty);
      }
    }
    return inter;
  }


  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    return val = ker*phix(kx,tx)*phiy(ky,ty);
  }


};

typedef BIOpKernel<CST,CST_OP,1,P1_1D,P1_1D> CST_1D_P1xP1;
typedef BIOpKernel<CST,CST_OP,1,P0_1D,P0_1D> CST_1D_P0xP0;
typedef BIOpKernel<CST,CST_OP,1,P0_1D,P1_1D> CST_1D_P0xP1;
typedef BIOpKernel<CST,CST_OP,1,P1_1D,P0_1D> CST_1D_P1xP0;
typedef BIOpKernel<CST,CST_OP,2,P1_2D,P1_2D> CST_2D_P1xP1;
typedef BIOpKernel<CST,CST_OP,2,P0_2D,P0_2D> CST_2D_P0xP0;
typedef BIOpKernel<CST,CST_OP,2,P0_2D,P1_2D> CST_2D_P0xP1;
typedef BIOpKernel<CST,CST_OP,2,P1_2D,P0_2D> CST_2D_P1xP0;
typedef BIOpKernel<CST,CST_OP,1,P2_1D,P2_1D> CST_1D_P2xP2;
typedef BIOpKernel<CST,CST_OP,2,P2_2D,P2_2D> CST_2D_P2xP2;














template <typename KernelType> class VIOp{

public:

  typedef typename KernelType::Trait         KernelTypeTrait;
  typedef typename KernelTypeTrait::MeshX    MeshX;
  typedef typename KernelTypeTrait::MeshY    MeshY;
  typedef typename KernelTypeTrait::Rdx      RdX;
  typedef typename KernelTypeTrait::Rdy      RdY;
  typedef typename KernelTypeTrait::EltX     EltX;
  typedef typename KernelTypeTrait::EltY     EltY;
  typedef typename KernelTypeTrait::MatType  MatType;

  static  const int dimx = KernelTypeTrait::dimx;
  static  const int dimy = KernelTypeTrait::dimy;
  typedef mat<dimx,dimx,Real>                JacX;
  typedef mat<dimy,dimy,Real>                JacY;
  typedef QuadVol<dimx,dimy>                 QuadType;

private:

  const MeshX&    meshx;
  const MeshY&    meshy;
  const Geometry& nodex;
  const Geometry& nodey;

  KernelType   ker;
  QuadType     qr;
  MatType      inter;
  int          rule;
  JacX         dx;
  JacY         dy;
  RdX          x0,x,ax[dimx+1];
  RdY          y0,y,ay[dimy+1];
  Cplx         val,val2;

public:
  VIOp<KernelType>(const MeshX& mx, const MeshY& my, const Real& k):
  meshx(mx), nodex(GeometryOf(mx)), ker(mx,my,k),
    meshy(my), nodey(GeometryOf(my)), qr(QR2) {
    for(int j=0; j<dimx; j++){ax[j+1][j]=1;}
    for(int j=0; j<dimy; j++){ay[j+1][j]=1;}
  };
  
  template <typename function_t>
  VIOp<KernelType>(const MeshX& mx, const MeshY& my, const Real& k, function_t func):
  meshx(mx), nodex(GeometryOf(mx)), ker(mx,my,k,func),
    meshy(my), nodey(GeometryOf(my)), qr(QR2) {
    for(int j=0; j<dimx; j++){ax[j+1][j]=1;}
    for(int j=0; j<dimy; j++){ay[j+1][j]=1;}
  };
  
  template <typename function_t>
  VIOp<KernelType>(const MeshX& mx, const MeshY& my, const Real& k, function_t func1, function_t func2):
  meshx(mx), nodex(GeometryOf(mx)), ker(mx,my,k,func1,func2),
    meshy(my), nodey(GeometryOf(my)), qr(QR2) {
    for(int j=0; j<dimx; j++){ax[j+1][j]=1;}
    for(int j=0; j<dimy; j++){ay[j+1][j]=1;}
  };
  

  void ChooseQuad(const int& jx, const int& jy){
    rule = 0; inter=0.; val=0.;
    EltX ex = meshx[jx];   EltY ey = meshy[jy];
    Elt<dimx,RdX> bx(ax);  Elt<dimy,RdY> by(ay);

    for(int p=0; p<dimx+1; p++){
      for(int q=rule; q<dimy+1; q++){
	if( &ex[p]==&ey[q] ){
	  Swap(ex,rule,p); Swap(bx,rule,p);
	  Swap(ey,rule,q); Swap(by,rule,q);
	  rule++; break;
	}
      }
    }
    
    //std::cout << rule << std::endl;
    dx = MatJac(bx); x0=bx[0];
    dy = MatJac(by); y0=by[0];
    
    
    
    //std::cout << det(dx) << "\t" << det(dy) << std::endl;
    //std::cout << dx << std::endl;
  }


  const MatType& operator()(const EltX& ex, const EltY& ey){
    return (*this)(nodex[ex][meshx],nodey[ey][meshy]);}


  const MatType& operator()(const int& jx, const int& jy){
    ChooseQuad(jx,jy);
    
    //std::cout << "Call Operator" << std::endl;
    if(rule != 0){
    const std::vector<RdX>&  s = qr.x(rule);
    const std::vector<RdY>&  t = qr.y(rule);
    const std::vector<Real>& w = qr.w(rule);
    ker.Assign(jx,jy);
    
    //std::cout << inter << std::endl;
    
    int size_w = w.size();
    
    //std::cout << dx << std::endl;
    
    //std::cout << dy << std::endl;
    
    double value = 0.;
    for(int j=0; j<size_w; j++){
      value += w[j];
      x = x0 + dx*s[j];
      y = y0 + dy*t[j];
      //std::cout << "Before" << std::endl;
      inter += w[j]*ker(x,y);
      
      //std::cout << "After" << std::endl;
    }
      
    
    }
    else{
    inter = 0. * inter;
    }
    return inter;
    
    
    
    
  }
  

  const Cplx& operator()(const N2& Ix, const N2& Iy){
    ChooseQuad(Ix[0],Iy[0]);
    if(rule == 0){
    const std::vector<RdX>&  s = qr.x(rule);
    const std::vector<RdY>&  t = qr.y(rule);
    const std::vector<Real>& w = qr.w(rule);
    ker.Assign(Ix[0],Iy[0]);
    for(int j=0; j<w.size(); j++){
      x = x0 + dx*s[j]; y = y0 + dy*t[j];
      val += w[j]*ker(x,y,Ix[1],Iy[1]);}
    }
    else{
      val = 0.;
    }
    
    return val;
    
  }


  const Cplx& operator()(const std::vector<N2>& vx,
			 const std::vector<N2>& vy){
    val2=0.;
    for(int ix=0; ix<vx.size(); ix++){
      for(int iy=0; iy<vy.size(); iy++){
	val2 += (*this)(vx[ix],vy[iy]);} }
    return val2;
  }


};






/*==============================
  GENERALITES SUR LES NOYAUX
  =============================*/

template <int,int,int,typename,typename>
class VIOpKernel;

template <int EqT, int OpT, int dim, typename PhiX, typename PhiY>
class VIOpKernelTraits{

public:
  static const int EqType   = EqT;
  static const int OpType   = OpT;
  static const int Dim      = dim;
  static const int dimx     = PhiX::dim;
  static const int dimy     = PhiY::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;

  typedef typename Jac<dimx>::Type     JacX;
  typedef typename Jac<dimy>::Type     JacY;
  typedef typename PhiX::Rd            Rdx;
  typedef typename PhiY::Rd            Rdy;
  typedef Mesh<dimx>                   MeshX;
  typedef Mesh<dimy>                   MeshY;
  typedef Elt<dimx>                    EltX;
  typedef Elt<dimy>                    EltY;
  typedef GradBasisFct<PhiX>           GradPhiX;
  typedef GradBasisFct<PhiY>           GradPhiY;
  typedef DivBasisFct<PhiX>            DivPhiX;
  typedef DivBasisFct<PhiY>            DivPhiY;
  typedef mat<nb_dof_x,nb_dof_y,Cplx>  MatType;

  typedef PhiX ShapeFctX;
  typedef PhiY ShapeFctY;  
  
};


/*=======================
        NOYAU CONSTANT
  =======================*/

template <int D, typename PhiX, typename PhiY>
class VIOpKernel<CST,CST_OP,D,PhiX,PhiY>{

public:
  typedef VIOpKernelTraits<CST,CST_OP,D,PhiX,PhiY> Trait;

private:
  const typename Trait::MeshX&    meshx;
  const typename Trait::MeshY&    meshy;
        typename Trait::MatType   inter;
                        PhiX      phix;
                        PhiY      phiy;
                        Real      ker;
                        Cplx      val;

public:
  VIOpKernel<CST,CST_OP,D,PhiX,PhiY>(const typename Trait::MeshX& mx,
				 const typename Trait::MeshY& my,
				 const Real& k=0.):
  meshx(mx), phix(mx), meshy(my), phiy(my) {};


  void Assign(const int& ix, const int& iy){
    phix.Assign(ix);
    phiy.Assign(iy);
    const typename Trait::EltX& ex=meshx[ix];
    const typename Trait::EltY& ey=meshy[iy];
    ker = DetJac(ex)*DetJac(ey);
  }


  const typename Trait::MatType&
  operator()(const typename Trait::Rdx& tx,
	     const typename Trait::Rdy& ty){
    for(int j=0; j<Trait::nb_dof_x; j++){
      for(int k=0; k<Trait::nb_dof_y; k++){
      
      /*
      std::cout << "----------------------" << std::endl;
      
      std::cout << j << " " << k << std::endl;
      std::cout << phix(j, tx) << std::endl;
      std::cout << phiy(k, ty) << std::endl;
      
      std::cout << tx << std::endl;
      std::cout << ty << std::endl;
      */
      
	inter(j,k) = ker*(phix(j,tx), phiy(k,ty));
	
      //std::cout << inter(j,k) << std::endl;
      }
    }
    return inter;
  }


  const Cplx& operator()(const typename Trait::Rdx& tx,
			 const typename Trait::Rdy& ty,
			 const int& kx, const int& ky){
    return val = ker*(phix(kx,tx),phiy(ky,ty));
  }


};



typedef VIOpKernel<CST,CST_OP,3,P1_3D,P1_3D>             CST_3D_P1xP1;

typedef VIOpKernel<CST,CST_OP,3,NED_3D,NED_3D>           CST_3D_NEDxNED;
typedef VIOpKernel<CST,CST_OP,3,Curl_NED_3D,Curl_NED_3D> CST_3D_cNEDxcNED;

typedef VIOpKernel<CST,CST_OP,3,P1_2D,P1_3D>             CST_2Dx3D_P1xP1;
typedef VIOpKernel<CST,CST_OP,3,P1_3D,P1_2D>             CST_3Dx2D_P1xP1;


}

#endif
