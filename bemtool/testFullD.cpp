// This file computes the block matrix for the combined Boundary-Volume formulation for Maxwell Transmission Problems.
// The block matrix corresponds to 4 blocks:
// - BIE b  : STF (or PMCHWT) Formulation with weighted operators.
// - TRACES : Traces of Volume Integral Operators
// - POT    : Layer Potentials
// - VOL    : Volume Integrla Operators (I - A)u + Kv 

// The unkwnown are assumed to be:
// alpha := (Rotated) Tangential trace of the Electric Field,                         \gammat u
// beta  := (scaled) (Rotated) Tangential trace of the Magnetic Field, i\omega \mu_1 \gammatau v
// u     := Electric Field
// v     := Magnetic Field 

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include "bemtool/tools.hpp"
#include "bemtool/parameters.hpp"
#include "castor/matrix.hpp"
#include "castor/linalg.hpp"
#include "castor/smatrix.hpp"
#include "castor/hmatrix.hpp"
#include <omp.h>

using namespace bemtool;
using namespace castor;


template <typename OperatorType>
struct Test{
  static inline void Launch(std::string filename, std::string extension, int ntest){

    
    
    const int NUM = 24;
    
    double    tol = 1e-6;
    
    Geometry node(filename + std::to_string(ntest) + extension);
    Mesh3D mesh; mesh.Load(node,1);
    
    Mesh2D mesh2D; mesh2D.Load(node,1); Orienting(mesh2D);
    int nb_elt = NbElt(mesh);
    int nb_elt2D = NbElt(mesh2D);
    
    std::cout << "nb_elt:\t" << nb_elt << std::endl;

    bool isConstant = (testEps(mesh) && testMu(mesh));
    // bool isConstant = false;
//---------------------------------------------------------------------------------//
    VIOp<MA_VL_NEDxCurlNED>   K1e(mesh, mesh,kappa1,p_e);

    VIOp<MA_TVL_NEDxNED>      K2e(mesh, mesh,kappa1,gradp_e);

    VIOp<MA_VL_NEDxNED>       Ve(mesh, mesh,kappa1,p_e);
    
    VIOp<MA_GVL_NEDxNED>      Ge(mesh, mesh,kappa1,tau_e);
//---------------------------------------------------------------------------------//
    
    VIOp<MA_SL_NEDxRT0>       SLe(mesh, mesh2D,kappa1,invepstilde,mutilde);
    VIOp<MA_SL_NEDxRT0>       SLm(mesh, mesh2D,kappa1,invmutilde,epstilde);
    
    VIOp<MA_DL_NEDxRT0>       DL(mesh, mesh2D,kappa1);
    
    
    VIOp<MA_VL_RT0xNED>       tVe(mesh2D,mesh,kappa1,p_e);
    
    VIOp<MA_GVL_RT0xNED>      tGe(mesh2D,mesh,kappa1,tau_e);
    
    VIOp<MA_VL_RT0xCurlNED>   tK1e(mesh2D,mesh,kappa1,p_e);

    VIOp<MA_TVL_RT0xNED>      tK2e(mesh2D,mesh,kappa1,gradp_e);
    
    
//---------------------------------------------------------------------------------//
    
    BIOp<EFIE_RT0xRT0>        V0(mesh2D,mesh2D,kappa0);
    BIOp<MFIE_RT0xRT0>        K0(mesh2D,mesh2D,kappa0);
    BIOp<MFIE_RT0xRT0>        Kp0(mesh2D,mesh2D,kappa0);
    BIOp<EFIE_RT0xRT0>        W0(mesh2D,mesh2D,kappa0);
    
    BIOp<EFIE_RT0xRT0>        V1(mesh2D,mesh2D,kappa1,invepstilde,mutilde);
    BIOp<MFIE_RT0xRT0>        K1(mesh2D,mesh2D,kappa1);
    BIOp<MFIE_RT0xRT0>        Kp1(mesh2D,mesh2D,kappa1);
    BIOp<EFIE_RT0xRT0>        W1(mesh2D,mesh2D,kappa1,invmutilde,epstilde);

//---------------------------------------------------------------------------------//
    
    
    
    
    LocalMatrix<NED_3D, Curl_NED_3D>      LocalX(mesh);
    LocalMatrix<NED_3D, NED_3D>           Local(mesh);
    LocalMatrix<Curl_NED_3D, Curl_NED_3D> LocalC(mesh);
    
    LocalVec<RT0_2D> Local0(mesh2D);

    LocalVec<NED_3D> Local00(mesh);
    
    
    Dof<NED_3D> dof(mesh);
    Dof<RT0_2D> dofRT(mesh2D);
    
    
    int nb_dof   = NbDof(dof);
    int nb_dof2D = NbDof(dofRT); 
    
    int nb_dof_loc = NED_3D::nb_dof_loc;
    
    int nb_dof_loc2D = RT0_2D::nb_dof_loc;
    
    std::cout << "nb_dof:\t" << nb_dof << std::endl;
    
    
    std::vector<double>  sol_real(nb_dof, 0.);
    
      
    // -----------------------------------------------------------------------//
    // --------------- NEAR FIELD MATRIX (SPARSE) ----------------------------//
    // -----------------------------------------------------------------------//
    std::vector<N2> Close, Close_Gamma, Close_2Dx3D, Close_3Dx2D;
    
    
    // -----------------------------------------------------------------------//
    // 3D
    for(int j = 0; j < nb_elt; j++){
      for(int k = 0; k < nb_elt; k++){
        if(CheckInteraction(j, k, mesh)){
          Close.push_back(N2_(j, k));
        }
      }
    }
    
    // -----------------------------------------------------------------------//
    // 2D
    for(int j = 0; j < nb_elt2D; j++){
      for(int k = 0; k < nb_elt2D; k++){
        if(CheckInteraction(j, k, mesh2D)){
          Close_Gamma.push_back(N2_(j, k));
        }
      }
    }
    
    // -----------------------------------------------------------------------//
    // 3D x 2D
    for(int j = 0; j < nb_elt; j++){
      for(int k = 0; k < nb_elt2D; k++){
        if(CheckInteraction(j, k, mesh, mesh2D)){
          Close_3Dx2D.push_back(N2_(j, k));
        }
      }
    }
    
    // -----------------------------------------------------------------------//
    // 2D x 3D
    for(int j = 0; j < nb_elt2D; j++){
      for(int k = 0; k < nb_elt; k++){
        if(CheckInteraction(j, k, mesh2D, mesh)){
          Close_2Dx3D.push_back(N2_(j, k));
        }
      }
    }
    
    // -----------------------------------------------------------------------//
    const int nb_close        = Close.size();
    const int nb_close2D      = Close_Gamma.size();
    const int nb_close3Dx2D   = Close_3Dx2D.size();
    const int nb_close2Dx3D   = Close_2Dx3D.size();

    std::cout << nb_close << "\t" << nb_close2D << "\t" << nb_close2Dx3D << "\t" << nb_close3Dx2D << std::endl;

    omp_set_num_threads(NUM);
    

    progress bar("VIE Matrix Assembly\t",nb_close/NUM);
    
    std::vector<std::size_t>  rowA;
    std::vector<std::size_t>  colA;
    std::vector<Cplx>         valAe;
    
    std::vector<Cplx>         valKe;
    

    
    if(!isConstant){
    

    
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      
      std::vector<std::size_t>  rowA_priv;
      std::vector<std::size_t>  colA_priv;
      std::vector<Cplx>         valAe_priv;
      
      std::vector<Cplx>         valKe_priv;
      
      
      C6x6 temp_e;

      #pragma omp for firstprivate(Ve, Ge, dof, Close)
      for(int l=0; l<nb_close; l++){
      
         int j = Close[l][0];
         int k = Close[l][1];
             
         temp_e = -kappa1 * kappa1 * Ve(j, k) + Ge(j, k);
          
         
          
         for(int nx = 0; nx < nb_dof_loc; nx++){
           for(int ny = 0; ny < nb_dof_loc; ny++){
             
	     if(std::abs(temp_e(nx,ny))>1e-10){	   
             rowA_priv.push_back(dof[j][nx]);
             colA_priv.push_back(dof[k][ny]);
              
             valAe_priv.push_back(temp_e(nx, ny));
	     }
            
            
           } 
          }
	         
        if(tid == 0){bar++;} 
      }
      
      
      #pragma omp critical
      {
       rowA.insert(rowA.end(), rowA_priv.begin(), rowA_priv.end());
       colA.insert(colA.end(), colA_priv.begin(), colA_priv.end());
       
       valAe.insert(valAe.end(), valAe_priv.begin(), valAe_priv.end());
       
       
      }
      
      
      if(tid == 0){bar.end();}
    
    }
    
    
    }
    else{
      bar.end();
    }
    
    
    smatrix<Cplx> Ae_sparse(nb_dof,nb_dof, rowA, colA, valAe);


    
    
    
    
    // -----------------------------------------------------------------------//
    
    progress bar00("BEM Matrix Assembly\t",nb_close2D/NUM);
    
    std::vector<std::size_t>  rowBEM;
    std::vector<std::size_t>  colBEM;
    
    std::vector<Cplx>         valVk;
    std::vector<Cplx>         valWk;
    std::vector<Cplx>         valKk;
    std::vector<Cplx>         valKpk;
    

    
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      
    std::vector<std::size_t>  rowBEM_priv;
    std::vector<std::size_t>  colBEM_priv;
    
    std::vector<Cplx>         valVk_priv;
    std::vector<Cplx>         valWk_priv;
    std::vector<Cplx>         valKk_priv;
    std::vector<Cplx>         valKpk_priv;
    
      
      
      C3x3 temp1,temp2,temp3,temp4;

      #pragma omp for firstprivate(V0,V1,W0,W1,K0,K1,Kp0,Kp1, dofRT, Close_Gamma)
      for(int l=0; l<nb_close2D; l++){
      
         int j = Close_Gamma[l][0];
         int k = Close_Gamma[l][1];
             
         temp1 = (mu0 / mu1)  * V0(j, k) + V1(j, k);
         temp2 = K0(j, k) + K1(j, k);
         temp3 = Kp0(j, k)+ Kp1(j, k); 
         temp4 = (mu1 / mu0) * kappa0 * kappa0*W0(j, k) + kappa1 * kappa1 *W1(j, k);
          
          
         for(int nx = 0; nx < nb_dof_loc2D; nx++){
           for(int ny = 0; ny < nb_dof_loc2D; ny++){ 
            
             rowBEM_priv.push_back(dofRT[j][nx]);
             colBEM_priv.push_back(dofRT[k][ny]);
             
             valVk_priv.push_back(temp1(nx, ny));
             valKk_priv.push_back(temp2(nx, ny));
             valKpk_priv.push_back(temp3(nx, ny));
             valWk_priv.push_back(temp4(nx, ny));
            
           } 
          }
	         
        if(tid == 0){bar00++;} 
      }
      
      
      #pragma omp critical
      {
       rowBEM.insert(rowBEM.end(), rowBEM_priv.begin(), rowBEM_priv.end());
       colBEM.insert(colBEM.end(), colBEM_priv.begin(), colBEM_priv.end());
       
       valVk.insert(valVk.end(), valVk_priv.begin(), valVk_priv.end());
       valKk.insert(valKk.end(), valKk_priv.begin(), valKk_priv.end());
       valKpk.insert(valKpk.end(), valKpk_priv.begin(), valKpk_priv.end());
       valWk.insert(valWk.end(), valWk_priv.begin(), valWk_priv.end());
      }
      
      
      if(tid == 0){bar00.end();}
    
    }
    
    
    smatrix<Cplx> V_sparse(nb_dof2D,nb_dof2D, rowBEM, colBEM, valVk);
    smatrix<Cplx> K_sparse(nb_dof2D,nb_dof2D, rowBEM, colBEM, valKk);
    smatrix<Cplx> Kp_sparse(nb_dof2D,nb_dof2D, rowBEM, colBEM, valKpk);
    smatrix<Cplx> W_sparse(nb_dof2D,nb_dof2D, rowBEM, colBEM, valWk);
    
    
    
    // -----------------------------------------------------------------------//
    
    progress barL("Boundary-Volume Matrix Assembly\t",nb_close2Dx3D/NUM);
    omp_set_num_threads(NUM);
    
    std::vector<std::size_t>  rowx1;
    std::vector<std::size_t>  colx1;

    std::vector<std::size_t>  rowx2;
    std::vector<std::size_t>  colx2;
   
    std::vector<std::size_t>  rowxx1;
    std::vector<std::size_t>  colxx1;
    
    std::vector<std::size_t>  rowxx2;
    std::vector<std::size_t>  colxx2;
   
    std::vector<Cplx>         valAde;
    std::vector<Cplx>         valAne;
    std::vector<Cplx>         valSLe;
    std::vector<Cplx>         valDLe;
    
    
    
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      
      std::vector<std::size_t>  rowx1_priv;
      std::vector<std::size_t>  colx1_priv;

      std::vector<std::size_t>  rowx2_priv;
      std::vector<std::size_t>  colx2_priv;

      std::vector<std::size_t>  rowxx1_priv;
      std::vector<std::size_t>  colxx1_priv;
    
      std::vector<std::size_t>  rowxx2_priv;
      std::vector<std::size_t>  colxx2_priv;
     
      std::vector<Cplx>         valAde_priv;
      std::vector<Cplx>         valAne_priv;
      std::vector<Cplx>         valSLe_priv;
      std::vector<Cplx>         valDLe_priv;
      
      
      
      
      C3x6 temp1e,temp2e;
      C6x3 temp3e,temp4;
      
      #pragma omp for firstprivate(tGe,tVe,tK1e,tK2e,SLe,DL,dof,dofRT)
      for(int l=0; l<nb_close2Dx3D; l++){


         int j = Close_2Dx3D[l][0];
         int k = Close_2Dx3D[l][1];
             
         if(!isConstant){
         temp1e = kappa1 * kappa1 * tVe(j,k) - tGe(j,k);
         
         temp2e = kappa1 * kappa1 * tK1e(j,k) + kappa1 * kappa1 * tK2e(j,k);
         }
         temp3e = 0. - SLe(k,j);

         temp4  = DL(k,j); 
         
         for(int nx = 0; nx < nb_dof_loc2D; nx++){
           for(int ny = 0; ny < nb_dof_loc; ny++){
            
             
            
             if(!isConstant){

             if(std::abs(temp1e(nx,ny)) > 1e-16){
             rowxx1_priv.push_back(dofRT[j][nx]);
             colxx1_priv.push_back(dof[k][ny]);
             valAde_priv.push_back(temp1e(nx, ny));}
             
             if(std::abs(temp2e(nx,ny)) > 1e-16){
             rowxx2_priv.push_back(dofRT[j][nx]);
             colxx2_priv.push_back(dof[k][ny]);
             valAne_priv.push_back(temp2e(nx, ny));}


             }


             if(std::abs(temp3e(ny,nx)) > 1e-16){
             rowx1_priv.push_back(dofRT[j][nx]);
             colx1_priv.push_back(dof[k][ny]);
             valSLe_priv.push_back(temp3e(ny, nx));}
            

             if(std::abs(temp4(ny,nx)) > 1e-16){
             rowx2_priv.push_back(dofRT[j][nx]);
             colx2_priv.push_back(dof[k][ny]);
             valDLe_priv.push_back(temp4(ny, nx));}
            
           } 
          }

  
        if(tid == 0){barL++;}
      }
      
      
      
      #pragma omp critical
      {
      
       rowx1.insert(rowx1.end(), rowx1_priv.begin(), rowx1_priv.end());
       colx1.insert(colx1.end(), colx1_priv.begin(), colx1_priv.end());

       rowx2.insert(rowx2.end(), rowx2_priv.begin(), rowx2_priv.end());
       colx2.insert(colx2.end(), colx2_priv.begin(), colx2_priv.end());

       rowxx1.insert(rowxx1.end(), rowxx1_priv.begin(), rowxx1_priv.end());
       colxx1.insert(colxx1.end(), colxx1_priv.begin(), colxx1_priv.end());
       
       rowxx2.insert(rowxx2.end(), rowxx2_priv.begin(), rowxx2_priv.end());
       colxx2.insert(colxx2.end(), colxx2_priv.begin(), colxx2_priv.end());
      

       valAde.insert(valAde.end(), valAde_priv.begin(), valAde_priv.end());
       valAne.insert(valAne.end(), valAne_priv.begin(), valAne_priv.end());
       valSLe.insert(valSLe.end(), valSLe_priv.begin(), valSLe_priv.end());
       valDLe.insert(valDLe.end(), valDLe_priv.begin(), valDLe_priv.end());
       
       
      }
    
      if(tid == 0){barL.end();}
    
    }
    
    
    smatrix<Cplx> tAe_sparse(nb_dof2D,nb_dof, rowxx1, colxx1, valAde);
    smatrix<Cplx> tKe_sparse(nb_dof2D,nb_dof, rowxx2, colxx2, valAne);
    smatrix<Cplx> SLe_sparse(nb_dof,nb_dof2D, colx1, rowx1, valSLe);
    smatrix<Cplx> DLe_sparse(nb_dof,nb_dof2D, colx2, rowx2, valDLe);
    
    
    
    // -----------------------------------------------------------------------//
    // --------------- END NEAR FIELD MATRIX ---------------------------------//
    // -----------------------------------------------------------------------//
    
    
    
    
    
    // -----------------------------------------------------------------------//
    // --------------- MASS MATRICES -----------------------------------------//
    // -----------------------------------------------------------------------//
    std::vector<std::size_t> row, col; 
    
    std::vector<Cplx> valX, val, valC; 
    
    for(int j=0; j<nb_elt; j++){
    
      C6x6 tempX = LocalX(j);
      C6x6 temp = Local(j);
      C6x6 tempC = LocalC(j);
      
    
      //std::cout << Local(j) << std::endl;
      for(int nx = 0; nx < nb_dof_loc; nx++){
        for(int ny = 0; ny < nb_dof_loc; ny++){
          row.push_back(dof[j][nx]);
          col.push_back(dof[j][ny]);
          
          valX.push_back(tempX(nx, ny));
          val.push_back(temp(nx, ny));
          valC.push_back(tempC(nx, ny));   
        }
      }
    }
    
    
    
    smatrix<Cplx> M0(nb_dof,nb_dof, row, col, val);
    smatrix<Cplx> Mx(nb_dof,nb_dof, row, col, valX);
    smatrix<Cplx> M1(nb_dof,nb_dof, row, col, valC);

    smatrix<Cplx> T = M0 + M1;


    // -----------------------------------------------------------------------//
    // ------------------ Boundary-to-Volume DOFs-----------------------------//
    // -----------------------------------------------------------------------//
    // std::vector<std::size_t> Pi, Pj;
    // std::size_t countt = 0;

    // for(int ii = 0; ii < nb_dof2D; ii++){
    
    //   for(int jj = 0; jj < nb_dof; jj++){

    //     if(norm(eval(Xb(ii, colon(0, 2))) - eval(X(jj, colon(0, 2))), "2") < 1e-10){

    //       Pi.push_back(ii);

    //       Pj.push_back(jj);

    //       countt++;
    //     }
    //   }
    // }

    // std::vector<Cplx> vec_P(nb_dof2D, 1.);
    // smatrix<Cplx> P_Gamma(nb_dof2D, nb_dof, Pi, Pj, vec_P);
    
    
    // -----------------------------------------------------------------------//
    // --------------- END MASS MATRICES -------------------------------------//
    // -----------------------------------------------------------------------//
    
    
    
    // -----------------------------------------------------------------------//
    // --------------------- H-MATRIX BEM -------------------------------//
    // -----------------------------------------------------------------------//
    matrix<double> Xb(nb_dof2D, 3);
    matrix<double> Yb(nb_dof2D, 3);
    
    for(int j = 0; j< nb_dof2D; j++){
    
    Xb(j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Xb(j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Xb(j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];
    
    
    Yb(j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Yb(j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Yb(j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];
    
    }

    matrix<double> Xbb = vertcat(Xb, Xb);
    matrix<double> Ybb = vertcat(Yb, Yb);
    
    castor::tic();
    
    // -----------------------------------------------------------------------//
    auto fctV = [dofRT, &V_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> V(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
        BIOp<EFIE_RT0xRT0> K0(MeshOf(dofRT),MeshOf(dofRT),kappa0); 
        BIOp<EFIE_RT0xRT0> K1(MeshOf(dofRT),MeshOf(dofRT),kappa1,invepstilde,mutilde); 
       
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = (mu0 / mu1) * K0(dofRT.ToElt(ix),dofRT.ToElt(iy)) + K1(dofRT.ToElt(ix),dofRT.ToElt(iy));
            }
        }
        
        V = full(eval(V_sparse(Ix, Iy)));
        M += V;
        return M;};
    
    // -----------------------------------------------------------------------//
    auto fctK = [dofRT, &K_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> V(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
        BIOp<MFIE_RT0xRT0> K0(MeshOf(dofRT),MeshOf(dofRT),kappa0); 
        BIOp<MFIE_RT0xRT0> K1(MeshOf(dofRT),MeshOf(dofRT),kappa1); 
       
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = K0(dofRT.ToElt(ix),dofRT.ToElt(iy)) + K1(dofRT.ToElt(ix),dofRT.ToElt(iy));
            }
        }
        
        V = full(eval(K_sparse(Ix, Iy)));
        M += V;
        return M;};
        
        
    // -----------------------------------------------------------------------//
    auto fctKp = [dofRT, &Kp_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> V(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
        BIOp<MFIE_RT0xRT0> K0(MeshOf(dofRT),MeshOf(dofRT),kappa0); 
        BIOp<MFIE_RT0xRT0> K1(MeshOf(dofRT),MeshOf(dofRT),kappa1); 
       
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = K0(dofRT.ToElt(ix),dofRT.ToElt(iy)) + K1(dofRT.ToElt(ix),dofRT.ToElt(iy));
            }
        }
        
        V = full(eval(Kp_sparse(Ix, Iy)));
        M += V;
        return M;};
        
    // -----------------------------------------------------------------------//
    auto fctW = [dofRT, &W_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> V(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
        BIOp<EFIE_RT0xRT0> K0(MeshOf(dofRT),MeshOf(dofRT),kappa0); 
        BIOp<EFIE_RT0xRT0> K1(MeshOf(dofRT),MeshOf(dofRT),kappa1,invmutilde,epstilde); 
       
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) =     (mu1/mu0) *    kappa0*kappa0* K0(dofRT.ToElt(ix),dofRT.ToElt(iy)) 
                           +                kappa1*kappa1* K1(dofRT.ToElt(ix),dofRT.ToElt(iy));
            }
        }
        
        V = full(eval(W_sparse(Ix, Iy)));
        M += V;
        return M;};
    // -----------------------------------------------------------------------//
    
    
    hmatrix<Cplx> Bh(Xb,Yb,tol,fctV);//disp(Bh);
    hmatrix<Cplx> Ah(Xb,Yb,tol,fctK);Ah *= (-1);  //disp(Ah);
    hmatrix<Cplx> Dh(Xb,Yb,tol,fctKp);Dh *= (-1); //disp(Dh);
    hmatrix<Cplx> Ch(Xb,Yb,tol,fctW); //disp(Ch);
    
    hmatrix<Cplx> BEMh(Ch, Dh, Ah, Bh);


    // hmatrix<Cplx> BEM0(Xb, Ybb, Ch, Dh);
    // BEM0.htranspose();

    // hmatrix<Cplx> BEM1(Xb, Ybb, Ah, Bh);
    // BEM1.htranspose();

    // hmatrix<Cplx> BEMh(Xbb, Ybb, BEM0, BEM1);
    // BEMh.htranspose();
    
    
    castor::toc();
    // -----------------------------------------------------------------------//
    // --------------------- END H-MATRIX ------------------------------------//
    // -----------------------------------------------------------------------//
    
    
    
    // -----------------------------------------------------------------------//
    // --------------------- H-MATRIX ELECTRIC  ------------------------------//
    // -----------------------------------------------------------------------//
    matrix<double> X(nb_dof, 3);
    matrix<double> Y(nb_dof, 3);
    
    for(int j = 0; j< nb_dof; j++){
    
    X(j, 0) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][0];
    X(j, 1) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][1];
    X(j, 2) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][2];
    
    
    Y(j, 0) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][0];
    Y(j, 1) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][1];
    Y(j, 2) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][2];
    
    }

    
    castor::tic();
    auto fctE = [isConstant, dof, &Ae_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> A1(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
        
        VIOp<MA_VL_NEDxNED>  Ve(MeshOf(dof),MeshOf(dof),kappa1,p_e);
        VIOp<MA_GVL_NEDxNED> Ge(MeshOf(dof),MeshOf(dof),kappa1,tau_e); 
        
        
        if(!isConstant){
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = kappa1*kappa1*Ve(dof.ToElt(ix),dof.ToElt(iy)) - Ge(dof.ToElt(ix),dof.ToElt(iy));

                if(std::isnan(std::abs(M(i,j)))){

                  std::cout << "Oops! Something is wrong" << std::endl;
                }
                
                
            }
        }
        }
        
        A1 = full(eval(Ae_sparse(Ix, Iy)));
        
        M += A1;
        
        return M;};
    
    
    // -----------------------------------------------------------------------//
    
    hmatrix<Cplx> Ie(X,Y,tol,M0);  

    hmatrix<Cplx> M0h(X,Y,tol,M0);  
    hmatrix<Cplx> VOLh = Ie;


    if(!isConstant){
    
    hmatrix<Cplx> VOLhh(X,Y,tol,fctE);

    VOLh += VOLhh;
  
    }
  


    //hmatrix<Cplx> VOLh = VOLh_Ae;
    
    castor::toc();

    // -----------------------------------------------------------------------//
    // --------------------- END H-MATRIX ------------------------------------//
    // -----------------------------------------------------------------------//
    
    
    
    
    
    
    // -----------------------------------------------------------------------//
    // --------------------- H-MATRIX TRACES E/M -----------------------------//
    // -----------------------------------------------------------------------//

    
    castor::tic();
    auto fcttE = [isConstant, dof, dofRT, &tAe_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> A1(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
        
        VIOp<MA_VL_RT0xNED>  tVe(MeshOf(dofRT),MeshOf(dof),kappa1,p_e);
        VIOp<MA_GVL_RT0xNED> tGe(MeshOf(dofRT),MeshOf(dof),kappa1,tau_e); 
        
        
        if(!isConstant){
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = kappa1*kappa1*tVe(dofRT.ToElt(ix),dof.ToElt(iy)) 
                       -               tGe(dofRT.ToElt(ix),dof.ToElt(iy));
                
            }
        }
        }
        A1 = full(eval(tAe_sparse(Ix, Iy)));
        
        M += A1;
        
        return M;};
    
    
    // -----------------------------------------------------------------------//
    
    
    
    auto fcttKe = [isConstant, dof, dofRT, &tKe_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> A1(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
      
        
        VIOp<MA_VL_RT0xCurlNED>   tK1e(MeshOf(dofRT),MeshOf(dof),kappa1,p_e);
        VIOp<MA_TVL_RT0xNED>      tK2e(MeshOf(dofRT),MeshOf(dof),kappa1,gradp_e);
        
        
        if(!isConstant){
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = kappa1 * kappa1 * tK1e(dofRT.ToElt(ix),dof.ToElt(iy)) 
                       + kappa1 * kappa1 * tK2e(dofRT.ToElt(ix),dof.ToElt(iy));
              
                
            }
        }
        }
        
        A1 = full(eval(tKe_sparse(Ix, Iy)));
        M += A1;
        
        return M;};
        
        
    // -----------------------------------------------------------------------//
    
    
    
    hmatrix<Cplx> TRACESh = hzeros<Cplx>(2*nb_dof2D, nb_dof);


    if(!isConstant){

    hmatrix<Cplx> TRACESh_Ae(Xb,Y,tol,fcttE);   TRACESh_Ae.htranspose();disp(TRACESh_Ae);
    hmatrix<Cplx> TRACESh_Ke(Xb,Y,tol,fcttKe);  TRACESh_Ke.htranspose();disp(TRACESh_Ke);
    
    hmatrix<Cplx> TRACEShh(X, Ybb, TRACESh_Ke, TRACESh_Ae);
    TRACEShh.htranspose();

    TRACESh = TRACEShh;

    }
    

    // TRACESh *= (-1.);
    //TRACESh_Ke.clear();
    //TRACESh_Ae.clear();

    // hmatrix<Cplx> TRACESh(TRACESh_Ke, TRACESh_Am, TRACESh_Ae, TRACESh_Km);
    
    castor::toc();
    
    disp(TRACESh);
    // -----------------------------------------------------------------------//
    
    
    
    
    
    // -----------------------------------------------------------------------//
    // --------------------- H-MATRIX POTENTIALS E/M -------------------------//
    // -----------------------------------------------------------------------//

    
    castor::tic();
    auto fctSLe = [dof, dofRT, &SLe_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> A1(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
        
        VIOp<MA_SL_NEDxRT0>  hSLe(MeshOf(dof),MeshOf(dofRT),kappa1,invepstilde,mutilde);
        
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = - hSLe(dof.ToElt(ix),dofRT.ToElt(iy));
                
                
            }
        }
        A1 = full(eval(SLe_sparse(Ix, Iy)));
        
        M += A1;
        
        return M;};
    
    
    // -----------------------------------------------------------------------//
    
    
    
    auto fctDLe = [dof, dofRT, &DLe_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> A1(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
      
        VIOp<MA_DL_NEDxRT0>  hDL(MeshOf(dof),MeshOf(dofRT),kappa1);
        
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = hDL(dof.ToElt(ix),dofRT.ToElt(iy));
                
                
            }
        }
        
        A1 = full(eval(DLe_sparse(Ix, Iy)));
        M += A1;
        
        return M;};
        
   
    // -----------------------------------------------------------------------//
    
    
    
    
    
    

    hmatrix<Cplx> SLeh(X,Yb,tol,fctSLe); SLeh *= -1.;
    hmatrix<Cplx> DLeh(X,Yb,tol,fctDLe);
    
    // disp(full(DLeh));
    
    hmatrix<Cplx> POTh(X, Ybb, DLeh, SLeh);


    Cplx kkappa1 = kappa1 * kappa1;
    hmatrix<Cplx> POThm(X, Ybb, kkappa1*SLeh, DLeh);
    


   // DLeh.clear();
   //SLeh.clear();

    // hmatrix<Cplx> POTh(SLeh, DLeh, DLmh, SLmh);
    
    castor::toc();
    
    // -----------------------------------------------------------------------//
    
    std::cout << "\n\n" << std::endl;    

    disp(BEMh);
    disp(TRACESh);
    disp(POTh);
    disp(VOLh);

    // disp(full(POTh));

    // disp("\n\n");
    // disp(full(VOLh));
    
    
    // castor::tic();
    
    // matrix<double> XXbb = vertcat(Xbb, X);
    // matrix<double> YYbb = vertcat(Ybb, Y);


    // hmatrix<Cplx> Matrix0(Xbb, YYbb, BEMh, TRACESh);
    // Matrix0.htranspose();


    // BEMh.clear();
    // TRACESh.clear();

    // hmatrix<Cplx> Matrix1(X, YYbb, POTh, VOLh);
    // Matrix1.htranspose();


    // POTh.clear();
    // VOLh.clear();


    // hmatrix<Cplx> Matrix(XXbb, YYbb, Matrix0, Matrix1);
    // Matrix.htranspose();
    // // hmatrix<Cplx> Matrix(BEMh, TRACESh, POTh, VOLh);
    

    // Matrix0.clear();
    // Matrix1.clear();

    // // std::cout << "Concatenation of Matrices: " << std::endl;

    // // matrix<Cplx> Test0 = horzcat(full(BEMh), full(TRACESh));disp("Test0");

    // // matrix<Cplx> Test1 = horzcat(full(POTh), full(VOLh));disp("Test1");

    // // matrix<Cplx> Test  = vertcat(Test0, Test1);disp("Test");

    // // std::cout << "Difference Matrix and Test: " << std::endl; 
    // // disp( norm(full(Matrix)-Test,"inf")/norm(Test,"inf") );




    // castor::toc();

    // disp(Matrix);
    

    // matrix<Cplx> Mat = full(Matrix);



    // std::ofstream fileM("matrix/Full-" + std::to_string(ntest) + ".txt");
    
    // for(int i=0; i<(nb_dof +2* nb_dof2D); i++){
    
    // for(int j=0; j<(nb_dof + 2*nb_dof2D); j++){
      
    
    //   if (fileM.is_open()){
    //     fileM << Mat(i, j) << '\n';
    //   }
      
    // }
    // fileM << "-----------------------------------------" << '\n';
    // }
    
    // fileM.close();

    std::cout << "Computing inverse ..." << std::endl;

    // castor::tic();
    // Matrix.hinv();
    // castor::toc();

    // disp(Matrix);
    // disp(inv(Matrix));
    // --------------------- END H-MATRIX ------------------------------------//
    // -----------------------------------------------------------------------//
    
    matrix<Cplx> b0(2*nb_dof2D, 1);

    matrix<Cplx> b1(2*nb_dof2D, 1);
    for(int j=0; j<nb_elt2D; j++){
      
      C3 vect0 = Local0(Hx, Hy, Hz, j); vect0 = (mu1/mu0) * iu * kappa0 * vect0; 
      
      C3 vect1 = Local0(Ex, Ey, Ez, j);// vect1 = (1./ iu) * vect1;
      
      
      for(int kj = 0; kj < nb_dof_loc2D; kj++){
      
        b0(dofRT[j][kj])              += vect0[kj];
        b0(nb_dof2D + dofRT[j][kj])   += vect1[kj]; 

        
        b1(dofRT[j][kj])              += vect1[kj];
        b1(nb_dof2D + dofRT[j][kj])   += vect0[kj];  
      }
    }
    

    castor::tic();    

    hmatrix<Cplx> BEMhm1 = inv(BEMh);



    std::function<matrix<Cplx>(matrix<Cplx>const&)> Schurfct, Am1fct;
    Am1fct = [&BEMhm1](matrix<Cplx> const B) {
        
      matrix<Cplx> X = mtimes(BEMhm1, B);

      return X;

    };
    
    Schurfct = [&M0, &BEMh, &VOLh, &POTh, &TRACESh, &isConstant](matrix<Cplx> const B) {
        
        // disp(norm(B, "2"));
        matrix<Cplx> Xa = mtimes(BEMh, B);

        // disp(norm(Xa, "2"));
        if(!isConstant){
        // disp(norm(B, "2"));
        matrix<Cplx> Xb = mtimes(POTh, B); //disp(norm(Xb, "2"));

        //disp("Solving Volume Integral Equation");
        if(norm(Xb, "2") > 1e-16){
        matrix<Cplx> Xbb= gmres(VOLh, Xb, 1e-6,500);
                     Xbb= mtimes(TRACESh, Xbb);
                     Xa -= Xbb;}                   
        }


        return Xa;

    };


    disp("Solving Schur Complement");
    matrix<Cplx> sol_boundary = gmres(Schurfct, b0, 1e-6,500,Am1fct,matrix<Cplx>(), 1);
    matrix<Cplx> sol0         = mtimes(POTh,  sol_boundary);
    matrix<Cplx> sol1         = mtimes(POThm, sol_boundary);

    disp("Computing Final Solution");
    matrix<Cplx> sol     =  gmres(VOLh, sol0, 1e-6,500);
    matrix<Cplx> solCurl =  gmres(VOLh, sol1, 1e-6,500); // WRONG!

    // matrix<Cplx> sol =  gmres(M0, sol0, 1e-6,500);
    
    castor::toc();
    



    std::cout << "L2 Norm of U: \t" << mtimes(transpose(real(sol)), mtimes(M0, real(sol))) << std::endl;

    std::cout << "Hcurl Norm U: \t" << mtimes(transpose(real(sol)), mtimes(T, real(sol))) << std::endl;



    std::cout << "L2 Norm of curl(U): \t" << mtimes(transpose(real(solCurl)), mtimes(M0, real(solCurl))) << std::endl;

    std::cout << "Hcurl Norm of curl(U): \t" << mtimes(transpose(real(solCurl)), mtimes(T, real(solCurl))) << std::endl;



    matrix<Cplx> UDir(nb_dof, 1);


    // disp("Maybe here?");

    for(int j=0; j<nb_elt; j++){
      

      C6 vect1 = Local00(Ex, Ey, Ez, j);

      // C6 vect0 = Local11(Hx, Hy, Hz, j); vect0 = iu * kappa0 * vect0;
      
      
      for(int kj = 0; kj < nb_dof_loc; kj++){
      
        UDir(dof[j][kj])   += vect1[kj];  

 
      }
    }


    // matrix<Cplx> solref = mtimes(inv(M0h), UDir); solref = real(solref);

    
    std::ofstream file("solutions/STF-VIE-Alpha-" + std::to_string(ntest) + ".txt");
    
    for(int n=0; n<nb_dof2D; n++){
      
      sol_real[n] = -sol_boundary(n).real();
    
      if (file.is_open()){
        file << sol_real[n] << '\n';
      }
      
    }
    
    file.close();



    std::vector<double> sol_real2(nb_dof2D, 0.);
    std::ofstream file2("solutions/STF-VIE-Beta-" + std::to_string(ntest) + ".txt");
    
    for(int n=0; n<nb_dof2D; n++){
      
      sol_real2[n] = -sol_boundary(nb_dof2D + n).real();
    
      if (file2.is_open()){
        file2 << sol_real2[n] << '\n';
      }
      
    }
    
    file2.close();
    

    std::vector<double> solU(nb_dof, 0.);
    
    std::ofstream fileU("solutions/STF-VIE-U-" + std::to_string(ntest) + ".txt");
    
    for(int n=0; n<nb_dof; n++){
      
      solU[n] = sol(n).real();
    
      if (fileU.is_open()){
        fileU << solU[n] << '\n';
      }
      
    }
    
    fileU.close();
    
    sol = real(sol);



    std::vector<double> sol_curlU(nb_dof, 0.);
    
    std::ofstream filecurlU("solutions/STF-VIE-curlU-" + std::to_string(ntest) + ".txt");
    
    for(int n=0; n<nb_dof; n++){
      
      sol_curlU[n] = solCurl(n).real();
    
      if (filecurlU.is_open()){
        filecurlU << sol_curlU[n] << '\n';
      }
      
    }
    
    filecurlU.close();
    
    solCurl = real(solCurl);

  // -------------------------------------------------------------------------------------------- //

    // double norms    = 0.;
    // double ref_norm = 0.;
    // double error    = 0.;

    // double norms0    = 0.;
    // double ref_norm0 = 0.;
    // double error0    = 0.;
    
    // matrix<Cplx> diff(nb_dof, 1);
    // for(int n = 0; n < nb_dof; n++){
    //   diff(n) = solU[n] - solref(n);
    // }
    
    // matrix<Cplx> refval(nb_dof, 1);
    // for(int n = 0; n < nb_dof; n++){
    //   refval(n) = solref(n);
    // }
    

    
    
    // matrix<Cplx> temp0 = mtimes(T, diff);
    
    // matrix<Cplx> temp00 = mtimes(M0, diff);
    
    
    // for(int n = 0; n < nb_dof; n++){
    //   // ref_norm += std::abs(sol_reference(n) * temp1(n));
    //   // norm += std::abs(sol_interp[n] * temp1(n));
    //   error += std::abs(diff(n) * temp0(n));

    //   // ref_norm0 += std::abs(sol_reference(n) * temp01(n));
    //   // norms0 += std::abs(sol_interp[n] * temp01(n));
    //   error0 += std::abs(diff(n) * temp00(n));
    // }
    
    // // std::cout << std::sqrt(error0) / std::sqrt(ref_norm0) << "\t" <<  std::sqrt(error) / std::sqrt(ref_norm) << std::endl;
    // std::cout << std::sqrt(error0) << "\t" <<  std::sqrt(error) << std::endl;

    // std::cout << castor::sqrt(mtimes(transpose(diff), mtimes(M0, diff))) << "\t" << castor::sqrt(mtimes(transpose(diff), mtimes(T, diff))) << std::endl;
    // std::cout << castor::sqrt(mtimes(transpose(refval), mtimes(M0, refval))) << "\t" << castor::sqrt(mtimes(transpose(refval), mtimes(T, refval))) << std::endl;
    // std::cout << castor::sqrt(mtimes(transpose(sol), mtimes(M0, sol))) << "\t" << castor::sqrt(mtimes(transpose(sol), mtimes(T, sol))) << std::endl;


    // std::cout << ref_norm0 << "\t" <<  ref_norm << std::endl;
  
  
  // -------------------------------------------------------------------------------------------- //



    
       
    std::vector<N5> boundary; 
    
    
    for(int j = 0; j < nb_elt2D; j++){
      Elt2D e = mesh2D[j]; Order(e);
      for(int k=0; k< nb_elt; k++){
        Elt3D E = mesh[k];
        
        array<4, Elt2D> faces = FacesOf(E);
        for(int l = 0; l < 4; l++){
          Elt2D face = faces[l]; Order(face);
          if(e == face){
            array<6, Elt1D> edges   = EdgesOf(E);
            array<3, Elt1D> edges_b = EdgesOf(e);
            N3 Idx;
            int count = 0;
            for(int elt = 0; elt<6; elt++){
              Elt1D edgeA = edges[elt]; Order(edgeA);
                
              for(int eltB = 0; eltB<3; eltB++){ 
                Elt1D edgeB = edges_b[eltB]; Order(edgeB);
                if(edgeA == edgeB){
                  Idx[count] = elt;  
                  
                  //std::cout << edgeA << std::endl;
                  //std::cout << "----------------------" << std::endl;
                  
                
                  count++; 
                }
              } 
              
            }
            
            
            boundary.push_back(N5_(j, k, Idx[0],Idx[1],Idx[2]));
            //std::cout<< boundary[j] << std::endl;
            
            //std::cout << mesh[k] << std::endl;
            //std::cout << "=========================" << std::endl;
          }
        }
        
      }
    }


    
    
    WriteEltVectGmsh(dofRT, "outputs/valsSTF-VIE-Alpha-" + std::to_string(ntest) + ".txt", sol_real);
    WriteEltVectGmsh(dofRT, "outputs/valsSTF-VIE-Beta-" + std::to_string(ntest) + ".txt", sol_real2);
    WriteEltVectGmsh(dof, mesh2D, boundary, "outputs/valsSTF-VIE-U-" + std::to_string(ntest) + ".txt", solU);
    WriteEltVectGmsh(dof, mesh2D, boundary, "outputs/valsSTF-VIE-curlU-" + std::to_string(ntest) + ".txt", sol_curlU);
   
    
   
    // int rule = 0;

// for(int nb = 0; nb < nb_close3Dx2D; nb++){
// int jj = Close_3Dx2D[nb][0];
// int kk = Close_3Dx2D[nb][1];

// disp(CheckInteraction(jj, kk, mesh, mesh2D));
// std::cout << std::setprecision(16) << std::endl;
// for(int q = 1; q < 12; q++){
//     QuadVol<3,2> quad_rule(q);
//     for(int rule = 0; rule < 4; rule++){
//     double sum = 0;
    
//     const std::vector<Real>& weights = quad_rule.w(rule);
//     const std::vector<R3>&  s = quad_rule.x(rule);  
//     const std::vector<R2>&  t = quad_rule.y(rule); 
//     for(int n = 0; n < weights.size(); n++){
//       sum += 1./norm2(MatJac(mesh[jj]) * s[n] - MatJac(mesh2D[kk]) * t[n]) * weights[n];

//     }
//     std::cout << sum << "\t";
//     }
//     std::cout << std::endl;
// }
// std::cout << std::endl;
//   }



    // std::cout << std::setprecision(16) << std::endl;
    // for(int n = 0; n < 5; n++){
    // std::cout << valSLe[n] << "\t" << valDLe[n] << "\t" << valVk[n] << std::endl;
    // }

    std::cout << "---------------------------------------" << std::endl;
  }
  
  
  
  
  
  
  static void Errors(std::string filename, std::string extension, const int ntest, const int nmax){
  
    
    Geometry node0(filename + std::to_string(ntest) + extension);
    Mesh3D mesh0; mesh0.Load(node0,1);
    
    Geometry node(filename + std::to_string(nmax) + extension);
    Mesh3D mesh; mesh.Load(node,1);
    
    
    
    int nb_elt  = NbElt(mesh);
    int nb_elt0 = NbElt(mesh0);
    

    Dof<NED_3D> dof(mesh);
    Dof<NED_3D> dof0(mesh0); NED_3D phi(mesh0);
    
    
    int nb_dof = NbDof(dof);
    int nb_dof0 = NbDof(dof0);
    int nb_dof_loc = NED_3D::nb_dof_loc;
    
    std::vector<double>  sol_real(nb_dof0, 0.);
    
    std::vector<double>  sol_interp(nb_dof, 0.);
    matrix<Cplx>         sol_reference(nb_dof, 1);
    
    // ---------------------------------------------------------------//
    // ----------- LOAD SOLUTIONS INTO VECTORS -----------------------//
    // ---------------------------------------------------------------//
    
    std::ifstream file0("solutions0/STF-VIE-U-" + std::to_string(ntest) + ".txt");
    double value;
    
    for(int n = 0; n < nb_dof0; n++){
      file0 >> value;
      sol_real[n] = value;
    }
    
    file0.close();
    
    
    
    
    std::ifstream file("solutions0/STF-VIE-U-" + std::to_string(nmax) + ".txt");

    // std::ifstream file("solutions1/FEM-BEM-U-" + std::to_string(nmax) + ".txt");

    
    for(int n = 0; n < nb_dof; n++){
      file >> value;
      sol_reference(n) = value;
    }
    
    file.close();
    
    
    
    // ---------------------------------------------------------------//
    // ----------- FINITE ELEMENT MATRIX ASSEMBLY --------------------//
    // ---------------------------------------------------------------//
    
    LocalMatrix<NED_3D, NED_3D>           Local(mesh);
    LocalMatrix<Curl_NED_3D, Curl_NED_3D> LocalC(mesh);
    
    // ---------------------------------------------------------------//
    // ---------------------------------------------------------------//
    
    
    std::vector<std::size_t> row0, col0, row1, col1; 
    
    std::vector<Cplx> val, valC; 
    
    for(int j=0; j<nb_elt; j++){
    
      C6x6 temp = Local(j);
      C6x6 tempC = LocalC(j);
      
    
      //std::cout << Local(j) << std::endl;
      for(int nx = 0; nx < nb_dof_loc; nx++){
        
      
        for(int ny = 0; ny < nb_dof_loc; ny++){
          
          if(std::abs(temp(nx,ny)) > 1e-8){
          row0.push_back(dof[j][nx]);
          col0.push_back(dof[j][ny]);
          
          val.push_back(temp(nx, ny));
          }

          if(std::abs(tempC(nx,ny)) > 1e-8){
          row1.push_back(dof[j][nx]);
          col1.push_back(dof[j][ny]);
          valC.push_back(tempC(nx, ny));
          }
        }
          
        
      
        
      }
      
      
    }
    
    
    
    smatrix<Cplx> M0(nb_dof,nb_dof, row0, col0, val);
    smatrix<Cplx> M1(nb_dof,nb_dof, row1, col1, valC);
    
    smatrix<Cplx> T = M0 + M1;
    
    
    
    
     
    
    std::vector< std::vector<int> > C2F; 
   
    
    
    for(int nelt = 0; nelt < nb_elt0; nelt++){
    
      int count = 0;
    
      Elt3D elt0 = mesh0[nelt];
      
      auto faces = FacesOf(elt0);
      
      std::vector<int> c2f(std::pow(8, (nmax - ntest)), 0);
      
      
      
      for(int n = 0; n < nb_elt; n++){
        
        Elt3D elt = mesh[n];
        
        R3 C1 = Ctr(elt);

        bool inside = 1;      
        for(int f = 0; f < 4; f++){
          R3 normal = NormalTo(faces[f]);
          inside *= ((normal, C1 - faces[f][0]) * (normal, elt0[f] - faces[f][0]) > 0. );
        }
        
        if(inside){
          c2f[count] = n;
          
          //std::cout << count << std::endl;
          count+=1;
        }
      
      }
      
      C2F.push_back(c2f);
      
      phi.Assign(nelt);
      
      R3x3 B = MatJac(elt0);
      
      for(int jj = 0; jj < std::pow(8, (nmax - ntest)); jj++){
      int idn = C2F[nelt][jj];
      
      Elt3D elt = mesh[idn];
      
      auto edges = EdgesOf(elt);
      
      for(int ll = 0; ll < nb_dof_loc; ll++){
      
      Elt1D edge_ll = edges[ll]; Order(edge_ll);
      
      R3 tangent = edge_ll[1] - edge_ll[0];         // ORIENTATION DEFINED THIS WAY --> ref. shapefct.hpp
      
      sol_interp[dof[idn][ll]] = 0.;
      
      R3   x = dof(idn)[ll] - elt0[0];
      
      R3   t = inv(B) * x;
      
      for(int ii = 0; ii < nb_dof_loc; ii++){
      sol_interp[dof[idn][ll]] += (sol_real[dof0[nelt][ii]] * phi(ii, t), tangent);
      }
      
      
      }
      
      
      }
      
    }
    
    
    // -----------------------------------------------------------------------//

    double norms    = 0.;
    double ref_norm = 0.;
    double error    = 0.;

    double norms0    = 0.;
    double ref_norm0 = 0.;
    double error0    = 0.;
    
    matrix<Cplx> diff(nb_dof, 1);
    for(int n = 0; n < nb_dof; n++){
      diff(n) = sol_interp[n] - sol_reference(n);
    }
    
    
    
    
    auto temp0 = mtimes(T, diff);
    
    auto temp1 = mtimes(T, sol_reference);


    auto temp00 = mtimes(M0, diff);
    
    auto temp01 = mtimes(M0, sol_reference);
    
    
    for(int n = 0; n < nb_dof; n++){
      ref_norm += std::abs(sol_reference(n) * temp1(n));
      // norm += std::abs(sol_interp[n] * temp1(n));
      error += std::abs(diff(n) * temp0(n));

      ref_norm0 += std::abs(sol_reference(n) * temp01(n));
      // norms0 += std::abs(sol_interp[n] * temp01(n));
      error0 += std::abs(diff(n) * temp00(n));
    }
    
    std::cout << std::sqrt(error0) / std::sqrt(ref_norm0) << "\t" <<  std::sqrt(error) / std::sqrt(ref_norm) << std::endl;
    // std::cout << norms0 << "\t" <<  norms << std::endl;
    std::cout << ref_norm0 << "\t" <<  ref_norm << std::endl;
  
  
  
  std::cout << "\n" << std::endl;
  
  
  }
  
  
  
  
  
  
  

};



int main(){

  std::string filename  = "mesh/oriented_cube";
  
  std::string extension = ".msh";
  
  const int nmax = 2;
  
  // Compute solutions
  for(int n = 0; n < nmax + 1; n++){
    Test<MA_KVL_3D_NEDxNED>::Launch(filename, extension, n);
  }
  
  
  // // Compute errors
  // for(int n = 0; n < nmax; n++){
  //   Test<MA_KVL_3D_NEDxNED>::Errors(filename, extension, n, nmax);
  // }
  
  
}


