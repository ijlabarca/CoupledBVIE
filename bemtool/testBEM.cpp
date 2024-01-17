// This file computes the block matrix for the Boundary Integral Equations for Maxwell Transmission Problems.
// The block matrix corresponds to
// - BIE : STF (or PMCHWT) Formulation.

// The unkwnowns are:
// alpha := (Rotated) Tangential trace of the Electric Field,                         \gammat u
// beta  := (scaled) (Rotated) Tangential trace of the Magnetic Field, i\omega \mu_1 \gammatau v


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


struct Test{
  static inline void Launch(std::string filename, std::string extension, int ntest){

    
    
    const int NUM = 8;
    
    double    tol = 1e-6;
    
    Geometry node(filename + std::to_string(ntest) + extension);

    // Geometry node("mesh/oriented_cube.msh");
    Mesh3D mesh; mesh.Load(node,1);
    
    Mesh2D mesh2D; mesh2D.Load(node,1); Orienting(mesh2D);
    int nb_elt = NbElt(mesh);
    int nb_elt2D = NbElt(mesh2D);
    
    std::cout << "nb_elt:\t" << nb_elt2D << std::endl;


    
//---------------------------------------------------------------------------------//
    
    BIOp<EFIE_RT0xRT0>        V0(mesh2D,mesh2D,kappa0);
    BIOp<MFIE_RT0xRT0>        K0(mesh2D,mesh2D,kappa0);
    BIOp<MFIE_RT0xRT0>        Kp0(mesh2D,mesh2D,kappa0);
    BIOp<EFIE_RT0xRT0>        W0(mesh2D,mesh2D,kappa0);
    
    BIOp<EFIE_RT0xRT0>        V1(mesh2D,mesh2D,kappa1);
    BIOp<MFIE_RT0xRT0>        K1(mesh2D,mesh2D,kappa1);
    BIOp<MFIE_RT0xRT0>        Kp1(mesh2D,mesh2D,kappa1);
    BIOp<EFIE_RT0xRT0>        W1(mesh2D,mesh2D,kappa1);

//---------------------------------------------------------------------------------//
    
    
    LocalVec<RT0_2D> Local0(mesh2D);
    
    
    Dof<RT0_2D> dofRT(mesh2D);
    
    
    int nb_dof2D = NbDof(dofRT); 
    
    int nb_dof_loc2D = RT0_2D::nb_dof_loc;
    
    std::cout << "nb_dof:\t" << nb_dof2D << std::endl;
    
    
    std::vector<double>  sol_real(2 * nb_dof2D, 0.);
    
      
    // -----------------------------------------------------------------------//
    // --------------- NEAR FIELD MATRIX (SPARSE) ----------------------------//
    // -----------------------------------------------------------------------//
    std::vector<N2> Close_Gamma;
    
    
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
    
    const int nb_close2D      = Close_Gamma.size();
    
    omp_set_num_threads(NUM);
    
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
             
         temp1 = (mu0 / mu1) * V0(j, k) + V1(j, k);
         temp2 = K0(j, k) + K1(j, k);
         temp3 = Kp0(j, k)+ Kp1(j, k); 
         temp4 = (mu1 / mu0) * kappa0*kappa0* W0(j, k) + kappa1 * kappa1 * W1(j, k);
          
          
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
    smatrix<Cplx> Kp_sparse(nb_dof2D,nb_dof2D,rowBEM, colBEM, valKpk);
    smatrix<Cplx> W_sparse(nb_dof2D,nb_dof2D, rowBEM, colBEM, valWk);
    
    // -----------------------------------------------------------------------//
    // --------------- END NEAR FIELD MATRIX ---------------------------------//
    // -----------------------------------------------------------------------//
    
    
    
    
    
    
    
    
    // -----------------------------------------------------------------------//
    // --------------------- H-MATRIX BEM -------------------------------//
    // -----------------------------------------------------------------------//
    matrix<double> Xb(static_cast<std::size_t>(nb_dof2D), 3);
    matrix<double> Yb(static_cast<std::size_t>(nb_dof2D), 3);
   
    matrix<double> Xbb(static_cast<std::size_t>(2*nb_dof2D), 3);
    matrix<double> Ybb(static_cast<std::size_t>(2*nb_dof2D), 3);
    

 
    for(int j = 0; j< nb_dof2D; j++){
    
    Xb(j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Xb(j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Xb(j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];
    
    
    Yb(j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Yb(j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Yb(j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];
    
 
    Xbb(j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Xbb(j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Xbb(j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];
    
        
    Ybb(j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Ybb(j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Ybb(j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];
    


    }
    

for(int j = 0; j< nb_dof2D; j++){
    
        
    Xbb(nb_dof2D + j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Xbb(nb_dof2D + j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Xbb(nb_dof2D + j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];

    Ybb(nb_dof2D + j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Ybb(nb_dof2D + j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Ybb(nb_dof2D + j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];
    

    }




    castor::tic();
    
    // -----------------------------------------------------------------------//
    auto fctV = [dofRT, &V_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> V(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
        BIOp<EFIE_RT0xRT0> K0(MeshOf(dofRT),MeshOf(dofRT),kappa0); 
        BIOp<EFIE_RT0xRT0> K1(MeshOf(dofRT),MeshOf(dofRT),kappa1); 
       
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) =  (mu0 / mu1) *K0(dofRT.ToElt(static_cast<int>(ix)),dofRT.ToElt(static_cast<int>(iy))) +  K1(dofRT.ToElt(static_cast<int>(ix)),dofRT.ToElt(static_cast<int>(iy)));
            }
        }
        
        V = full(eval(V_sparse(Ix, Iy)));
        M += V;
        return M;};
    
    // // -----------------------------------------------------------------------//
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
        BIOp<EFIE_RT0xRT0> K1(MeshOf(dofRT),MeshOf(dofRT),kappa1); 
       
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) =   (mu1/mu0) * kappa0 * kappa0 *K0(dofRT.ToElt(ix),dofRT.ToElt(iy)) 
                          +  kappa1 * kappa1*K1(dofRT.ToElt(ix),dofRT.ToElt(iy));
            }
        }
        
        V = full(eval(W_sparse(Ix, Iy)));
        M += V;

        return M;};
    // -----------------------------------------------------------------------//
     
    
    hmatrix<Cplx> Bh(Xb,Yb,tol,fctV);
    hmatrix<Cplx> Ah(Xb,Yb,tol,fctK); Ah *= (-1);
    hmatrix<Cplx> Dh(Xb,Yb,tol,fctKp);Dh *= (-1);
    hmatrix<Cplx> Ch(Xb,Yb,tol,fctW);
    

    
    hmatrix<Cplx> BEMh(Ch, Dh, Ah, Bh);


    disp(BEMh);
    castor::toc();


/*
    std::cout << "Difference BEMh and Concat: " << std::endl; 
    disp( norm(full(BEMh)-full(Concat),"inf")/norm(full(BEMh),"inf") );



    matrix<Cplx> Test = horzcat(transpose(Test0), transpose(Test1));
    matrix<Cplx> Testt = transpose(Test);

    std::cout << "Difference BEMh and Test: " << std::endl; 
    disp( norm(full(BEMh)-Test,"inf")/norm(full(BEMh),"inf") );


*/
    castor::tic();
    BEMh.hinv();
    castor::toc();

    disp(BEMh);
    // // -----------------------------------------------------------------------//
    // // --------------------- END H-MATRIX ------------------------------------//
    // // -----------------------------------------------------------------------//
    
    
    matrix<Cplx> b0(2*nb_dof2D, 1);
    for(int j=0; j<nb_elt2D; j++){
      
      C3 vect0 = Local0(Hx, Hy, Hz, j); vect0 = (mu1/mu0) * iu * kappa0 * vect0; 
      
      C3 vect1 = Local0(Ex, Ey, Ez, j);
      
      
      for(int kj = 0; kj < nb_dof_loc2D; kj++){
      
        b0(dofRT[j][kj])              += vect0[kj]; 
        b0(nb_dof2D + dofRT[j][kj])   += vect1[kj];  
      }
    }
    
    // hmatrix<Cplx> invBEMh = inv(BEMh);



    // std::cout << "Second inverse" << std::endl;
    // std::cout << "Works" << std::endl;

    castor::tic();
    // matrix<Cplx> sol = linsolve(BEMh, b0); // LU Factorization
    matrix<Cplx> sol = mtimes(BEMh, b0); // In-place inverse with BEMh.hinv(), then multiply
    // matrix<Cplx> sol = gmres(BEMh, b0, 1e-6, 1000); // GMRES iterative solver, no precond.
    castor::toc();
    
    
    std::ofstream file("solutions/STF-" + std::to_string(ntest) + ".txt");
    
    for(int n=0; n<nb_dof2D; n++){
      
      sol_real[n] = -sol(n).real();
    
      if (file.is_open()){
        file << sol_real[n] << '\n';
      }
      
    }
    
    file.close();
    
    WriteEltVectGmsh(dofRT, "outputs/valsSTF-" + std::to_string(ntest) + ".txt", sol_real);
    
    
   
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
    
    std::ifstream file0("solutions/STF-" + std::to_string(ntest) + ".txt");
    double value;
    
    for(int n = 0; n < nb_dof0; n++){
      file0 >> value;
      sol_real[n] = value;
    }
    
    file0.close();
    
    
    
    
    std::ifstream file("solutions/STF-" + std::to_string(nmax) + ".txt");

    
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
    
    
    std::vector<std::size_t> row, col; 
    
    std::vector<Cplx> val, valC; 
    
    for(int j=0; j<nb_elt; j++){
    
      C6x6 temp = Local(j);
      C6x6 tempC = LocalC(j);
      
    
      //std::cout << Local(j) << std::endl;
      for(int nx = 0; nx < nb_dof_loc; nx++){
        
      
        for(int ny = 0; ny < nb_dof_loc; ny++){
          row.push_back(dof[j][nx]);
          col.push_back(dof[j][ny]);
          
          val.push_back(temp(nx, ny));
          valC.push_back(tempC(nx, ny));
          
        }
      
        
      }
      
      
    }
    
    
    
    smatrix<Cplx> M0(nb_dof,nb_dof, row, col, val);
    smatrix<Cplx> M1(nb_dof,nb_dof, row, col, valC);
    
    smatrix<Cplx> T = M0;
    
    
    
    
     
    
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
      
      R3 tangent = edge_ll[0] - edge_ll[1];         // ORIENTATION DEFINED THIS WAY --> ref. shapefct.hpp
      
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

    double ref_norm = 0.;
    double error    = 0.;
    
    matrix<Cplx> diff(nb_dof, 1);
    for(int n = 0; n < nb_dof; n++){
      diff(n) = sol_interp[n] - sol_reference(n);
    }
    
    
    
    
    auto temp0 = mtimes(T, diff);
    
    auto temp1 = mtimes(T, sol_reference);
    
    
    for(int n = 0; n < nb_dof; n++){
      ref_norm += std::abs(sol_reference(n) * temp1(n));
      error += std::abs(diff(n) * temp0(n));
    }
    
    std::cout <<  std::sqrt(error) / std::sqrt(ref_norm) << std::endl;
  
  
  
  
  }
  
  
  
  
  
  
  

};



int main(){


  std::string filename  = "mesh/oriented_cube";
  
  std::string extension = ".msh";
  
  const int nmax = 2;
  
  // Compute solutions
  for(int n = 0; n < nmax + 1; n++){
    Test::Launch(filename, extension, n);
  }
  
  
  // // Compute errors
  // for(int n = 0; n < nmax; n++){
  //   Test::Errors(filename, extension, n, nmax);
  // }
  
  
}


