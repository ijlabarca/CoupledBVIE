//===================================================================
//
//  Copyright 2017 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  (see ../LICENSE.txt)
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef BEMTOOL_FEM_FEM_HPP
#define BEMTOOL_FEM_FEM_HPP
#include "../quadrature/quad.hpp"
#include "../mesh/mesh.hpp"
#include "shapefct.hpp"


namespace bemtool {

template <typename T>
Real ScalarProd(const T&, const T&);


template <>
Real ScalarProd(const Real& u, const Real& v){
  return u*v;}

template <int D>
Real ScalarProd(const array<D,Real>& u, const array<D,Real>& v){
  return (u,v);}
    




template <typename PhiX,typename PhiY = PhiX>
class LocalMatrix{

  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;
  typedef mat<nb_dof_x,nb_dof_y,Real>  ReturnType;
  typedef Mesh<dim>                    MeshType;
  typedef mat<3,dim,Real>              JacX;
  typedef typename PhiX::Rd            Rd;

private:
  const MeshType&     mesh;
  ReturnType          inter;
  PhiX                phix;
  PhiY                phiy;
  Real                h;
  JacX                dx;
  R3                  x0, x;
  Cplx                val, val2;
  QuadFEM<dim>        qr;
  const std::vector<Rd>&   t;
  const std::vector<Real>& dt;

public:
  LocalMatrix<PhiX,PhiY>(const MeshType& m):
  mesh(m), phix(m), phiy(m),
    qr(10), t(points(qr)), dt(weights(qr)) {};



  const ReturnType& operator()(const int& l) {
    const Elt<dim>& e = mesh[l];
    phix.Assign(l);
    phiy.Assign(l);
    inter = 0.;
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
      for(int kx = 0; kx<nb_dof_x; kx++){
	      for(int ky = 0; ky<nb_dof_y; ky++){
	        inter(kx,ky)  += (phix(kx,t[j]),phiy(ky,t[j]))*h*dt[j];
	      }
	
      }
      
    }
    
    return inter;
  }



  template <typename fct_t>
  const ReturnType& operator()(const fct_t& f, const int& l) {
    const Elt<dim>& e = mesh[l];
    phix.Assign(l);
    phiy.Assign(l);
    inter = 0.;
    h     = DetJac(e);
    dx = MatJac(e); x0=e[0];
    for(int j=0; j<t.size(); j++){
    
      x     = x0 + dx*t[j];
    
      for(int kx = 0; kx<nb_dof_x; kx++){
      
	for(int ky = 0; ky<nb_dof_y; ky++){
	
	   inter(kx,ky)  += (phix(kx,t[j]),phiy(ky,t[j]))*f(x)*h*dt[j];
	}
	
      }
      
    }
    
    return inter;
  }


 

  const Cplx& operator()(const N2& Ix, const N2& Iy){
    
    const Elt<dim>& e = mesh[Ix[0]];
    phix.Assign(Ix[0]);
    phiy.Assign(Iy[0]);
    val = 0.;
    h     = DetJac(e);
    dx = MatJac(e); x0=e[0];
    for(int j=0; j<t.size(); j++){
      x = x0 + dx*t[j];
      val  += (phix(Ix[1],t[j]),phiy(Iy[1],t[j]))*h*dt[j];
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









template <>
class LocalMatrix<RT0_2D>{

  typedef RT0_2D                      PhiX;
  typedef RT0_2D                      PhiY;
  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;
  typedef mat<nb_dof_x,nb_dof_y,Real>  ReturnType;
  typedef Mesh<dim>                    MeshType;
  typedef mat<3,dim,Real>              JacX;
  typedef typename PhiX::Rd            Rd;

private:
  const MeshType&     mesh;
  ReturnType          inter;
  PhiX                phix;
  PhiY                phiy;
  const std::vector<R3>&   normal;
  R3                       nx;
  Real                h;
  JacX                dx;
  R3                  x0, x;
  Cplx                val, val2;
  QuadFEM<dim>        qr;
  const std::vector<Rd>&   t;
  const std::vector<Real>& dt;

public:


  LocalMatrix<RT0_2D,RT0_2D>(const MeshType& m):
  mesh(m), normal(NormalTo(m)), phix(m), phiy(m),
    qr(10), t(points(qr)), dt(weights(qr)) {};

  const ReturnType& operator()(const int& l) {
    const Elt<dim>& e = mesh[l];
    phix.Assign(l);
    phiy.Assign(l);
    inter = 0.;
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
      for(int kx = 0; kx<nb_dof_x; kx++){
	      for(int ky = 0; ky<nb_dof_y; ky++){
	
	        inter(kx,ky)  += (phix(kx,t[j]),phiy(ky,t[j]))*h*dt[j];
	      }
	
      }
      
    }
    
    return inter;
  }


  const ReturnType& operator()(const int& l, const std::string& crossprod) {
    const Elt<dim>& e = mesh[l];
    nx = normal[l];
    phix.Assign(l);
    phiy.Assign(l);
    inter = 0.;
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
    
    
      for(int kx = 0; kx<nb_dof_x; kx++){
      
	for(int ky = 0; ky<nb_dof_y; ky++){
	
	   inter(kx,ky)  += (phix(kx,t[j]),vprod(phiy(ky,t[j]), nx ))*h*dt[j];
	}
	
      }
      
    }
    
    return inter;
  }

  template <typename fct_t>
  const ReturnType& operator()(const fct_t& f, const int& l) {
    const Elt<dim>& e = mesh[l];
    phix.Assign(l);
    phiy.Assign(l);
    inter = 0.;
    h     = DetJac(e);
    dx = MatJac(e); x0=e[0];
    for(int j=0; j<t.size(); j++){
    
      x     = x0 + dx*t[j];
    
      for(int kx = 0; kx<nb_dof_x; kx++){
      
	for(int ky = 0; ky<nb_dof_y; ky++){
	
	   inter(kx,ky)  += (phix(kx,t[j]),phiy(ky,t[j]))*f(x)*h*dt[j];
	}
	
      }
      
    }
    
    return inter;
  }


 

  const Cplx& operator()(const N2& Ix, const N2& Iy){
    
    const Elt<dim>& e = mesh[Ix[0]];
    phix.Assign(Ix[0]);
    phiy.Assign(Iy[0]);
    val = 0.;
    h     = DetJac(e);
    dx = MatJac(e); x0=e[0];
    for(int j=0; j<t.size(); j++){
        
    x = x0 + dx*t[j];
    val  += (phix(Ix[1],t[j]),phiy(Iy[1],t[j]))*h*dt[j];
    
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



template <typename PhiX>
class LocalVec{

  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  typedef array<nb_dof_x,Cplx>         ReturnType;
  typedef Mesh<dim>                    MeshType;
  typedef typename PhiX::Rd            Rd;
  typedef mat<3,dim,Real>              JacX;

private:
  const MeshType&     mesh;
  const std::vector<R3>&   normal;
  R3                       nx;
  ReturnType          inter;
  PhiX                phix;
  Real                h;
  QuadFEM<dim>        qr;
  JacX                dx;
  R3                  x0, x;
  Cplx                val, val2;
  const std::vector<Rd>&   t;
  const std::vector<Real>& dt;
        

public:
  LocalVec<PhiX>(const MeshType& m):
  mesh(m), normal(NormalTo(m)), phix(m),
    qr(10), t(points(qr)), dt(weights(qr)) {};

  
  template <typename fct_t>
  const ReturnType& operator()(const fct_t& f, const int& l){
    const Elt<dim>& e = mesh[l];
    
    phix.Assign(l);
    dx = MatJac(e); x0=e[0];
    for(int kx = 0; kx<nb_dof_x; kx++){
      inter[kx] = 0.;
    }
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
      x     = x0 + dx*t[j];
      for(int kx = 0; kx<nb_dof_x; kx++){
	      inter[kx]+= phix(kx,t[j])*f(x)*h*dt[j];
      }
      
    }
    return inter;
  }



  template <typename fct_t>
  const ReturnType& operator()(const fct_t& fx,const fct_t& fy,const fct_t& fz, const int& l) {
    const Elt<dim>& e = mesh[l];
    nx = normal[l];
    phix.Assign(l);
    dx = MatJac(e); x0=e[0];
    
    for(int kx = 0; kx<nb_dof_x; kx++){
      inter[kx] = 0.;
    }
    
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
    
      
      x     = x0 + dx*t[j];
      for(int kx = 0; kx<nb_dof_x; kx++){
        Cplx dfn = fx( x ) * nx[0] +
              fy( x ) * nx[1] +
              fz( x ) * nx[2];
	      inter[kx]+= phix(kx,t[j])* dfn*h*dt[j];
      }
    }
    return inter;
  }

};






template <>
class LocalVec<RT0_2D>{
  
  typedef RT0_2D                       PhiX;
  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  typedef array<nb_dof_x,Cplx>         ReturnType;
  typedef Mesh<dim>                    MeshType;
  typedef typename RT0_2D::Rd            Rd;
  typedef mat<3,dim,Real>              JacX;

private:
  const MeshType&     mesh;
  const std::vector<R3>&   normal;
  R3                       nx;
  ReturnType          inter;
  PhiX                phix;
  Real                h;
  QuadFEM<dim>        qr;
  JacX                dx;
  R3                  x0, x;
  const std::vector<Rd>&   t;
  const std::vector<Real>& dt;
        

public:
  LocalVec<RT0_2D>(const MeshType& m):
  mesh(m), normal(NormalTo(m)), phix(m),
    qr(10), t(points(qr)), dt(weights(qr)) {};


  template <typename fct_t>
  const ReturnType& operator()(const fct_t& fx,const fct_t& fy,const fct_t& fz, const int& l) {
    const Elt<dim>& e = mesh[l];
    nx = normal[l];
    phix.Assign(l);
    dx = MatJac(e); x0=e[0];
    
    for(int kx = 0; kx<nb_dof_x; kx++){
      inter[kx] = 0.;
    }
    
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
    
      
      x     = x0 + dx*t[j];
      for(int kx = 0; kx<nb_dof_x; kx++){
        C3 F;
        F[0] = fx( x );
        F[1] = fy( x );
        F[2] = fz( x );
	      inter[kx]+= (phix(kx,t[j]),F)*h*dt[j];
      }
    }
    return inter;
  }

};



template <>
class LocalVec<NED_3D>{
  
  typedef NED_3D                       PhiX;
  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  typedef array<nb_dof_x,Cplx>         ReturnType;
  typedef Mesh<dim>                    MeshType;
  typedef typename PhiX::Rd            Rd;
  typedef mat<3,dim,Real>              JacX;

private:
  const MeshType&     mesh;
  ReturnType          inter;
  PhiX                phix;
  Real                h;
  QuadFEM<dim>        qr;
  JacX                dx;
  R3                  x0, x;
  const std::vector<Rd>&   t;
  const std::vector<Real>& dt;
        

public:
  LocalVec<NED_3D>(const MeshType& m):
  mesh(m), phix(m), qr(15), t(points(qr)), dt(weights(qr)) {};


  template <typename fct_t>
  const ReturnType& operator()(const fct_t& fx,const fct_t& fy,const fct_t& fz, const int& l) {
    const Elt<dim>& e = mesh[l];
    phix.Assign(l);
    dx = MatJac(e); x0=e[0];
    
    for(int kx = 0; kx<nb_dof_x; kx++){
      inter[kx] = 0.;
    }
    
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
    
      
      x     = x0 + dx*t[j];
      for(int kx = 0; kx<nb_dof_x; kx++){
        C3 F;
        F[0] = fx( x );
        F[1] = fy( x );
        F[2] = fz( x );
	      inter[kx]+= (phix(kx,t[j]),F)*h*dt[j];
      }
    }
    return inter;
  }

};



template <>
class LocalVec<Curl_NED_3D>{
  
  typedef NED_3D                       PhiX;
  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  typedef array<nb_dof_x,Cplx>         ReturnType;
  typedef Mesh<dim>                    MeshType;
  typedef typename PhiX::Rd            Rd;
  typedef mat<3,dim,Real>              JacX;

private:
  const MeshType&     mesh;
  ReturnType          inter;
  PhiX                phix;
  Real                h;
  QuadFEM<dim>        qr;
  JacX                dx;
  R3                  x0, x;
  const std::vector<Rd>&   t;
  const std::vector<Real>& dt;
        

public:
  LocalVec<Curl_NED_3D>(const MeshType& m):
  mesh(m), phix(m), qr(15), t(points(qr)), dt(weights(qr)) {};


  template <typename fct_t>
  const ReturnType& operator()(const fct_t& fx,const fct_t& fy,const fct_t& fz, const int& l) {
    const Elt<dim>& e = mesh[l];
    phix.Assign(l);
    dx = MatJac(e); x0=e[0];
    
    for(int kx = 0; kx<nb_dof_x; kx++){
      inter[kx] = 0.;
    }
    
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
    
      
      x     = x0 + dx*t[j];
      for(int kx = 0; kx<nb_dof_x; kx++){
        C3 F;
        F[0] = fx( x );
        F[1] = fy( x );
        F[2] = fz( x );
	      inter[kx]+= (phix(kx,t[j]),F)*h*dt[j];
      }
    }
    return inter;
  }

};






// P0
template <int dimX,int dimY>
class LocalMatrix<BasisFct<0,dimX>, BasisFct<0,dimY>> {
  typedef BasisFct<0,dimX> PhiX;
  typedef BasisFct<0,dimY> PhiY;
  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;
  typedef mat<nb_dof_x,nb_dof_y,Real>  ReturnType;
  typedef Mesh<dim>                    MeshType;

private:
  const MeshType&     mesh;
  ReturnType          inter;

public:
  LocalMatrix<PhiX,PhiY>(const MeshType& m): mesh(m) {};

  const ReturnType& operator()(const int& l) {
    const Elt<dim>& e = mesh[l];
    inter(0,0) = Vol(e);
    return inter;
  }

};

/*
// P1
template <int dimX,int dimY>
class LocalMatrix<BasisFct<1,dimX>, BasisFct<1,dimY>> {
  typedef BasisFct<1,dimX> PhiX;
  typedef BasisFct<1,dimY> PhiY;
  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  static const int nb_dof_y = PhiY::nb_dof_loc;
  typedef mat<nb_dof_x,nb_dof_y,Real>  ReturnType;
  typedef Mesh<dim>                    MeshType;

private:
  const MeshType&     mesh;
  ReturnType          inter;

public:
  LocalMatrix<PhiX,PhiY>(const MeshType& m): mesh(m) {};

  const ReturnType& operator()(const int& l) {
    const Elt<dim>& e = mesh[l];
    inter = MassP1(e);
    return inter;
  }

};


*/






/*

template <typename PhiX>
class Projection{

  static const int dim      = PhiX::dim;
  static const int nb_dof_x = PhiX::nb_dof_loc;
  typedef mat<nb_dof_x,1,Cplx>  ReturnType;
  typedef Mesh<dim>                    MeshType;
  typedef typename PhiX::Rd            Rd;

private:
  const MeshType&     mesh;
  ReturnType          inter;
  PhiX                phix;
  Real                h;
  QuadFEM<dim>        qr;
  const std::vector<Rd>&   t;
  const std::vector<Real>& dt;

public:
  Projection<PhiX>(const MeshType& m):
  mesh(m), phix(m),
    qr(5), t(points(qr)), dt(weights(qr)) {};

  template <typename function_t>
  const ReturnType& operator()(const function_t& f, const int& l) {
    const Elt<dim>& e = mesh[l];
    phix.Assign(l);
    inter = 0.;
    h     = DetJac(e);
    for(int j=0; j<t.size(); j++){
      for(int kx = 0; kx<nb_dof_x; kx++){
	  inter(kx,1)+= ScalarProd( phix(kx,t[j]),f(t[j]) )*h*dt[j];
      }
    }
    return inter;
  }

};

*/





}

#endif
