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
    
    double    tol = -1e-6;

    
    Geometry node(filename + std::to_string(ntest) + extension);

    // Geometry node("mesh/tetra0.msh");

    // Geometry node2("mesh/tetrab0.msh");
    Mesh3D mesh; mesh.Load(node,-1);
    
    Mesh2D mesh2D; mesh2D.Load(node,-1); Orienting(mesh2D);

    // Mesh2D mesh2D; mesh2D.Load(node2,1); Orienting(mesh2D);
    int nb_elt = NbElt(mesh);
    int nb_elt2D = NbElt(mesh2D);
    
    std::cout << "nb_elt:\t" << nb_elt << std::endl;

    bool isConstant = (testEps(mesh) && testMu(mesh));
    // bool isConstant = false;
    
    
//---------------------------------------------------------------------------------//
    
    BIOp<EFIE_RT0xRT0>        V0(mesh2D,mesh2D,kappa0);
    BIOp<MFIE_RT0xRT0>        K0(mesh2D,mesh2D,kappa0);
    

//---------------------------------------------------------------------------------//
    
    
    LocalMatrix<NED_3D, NED_3D>           Local(mesh);

    LocalMatrix<NED_3D, Curl_NED_3D>      LocalX(mesh);
    LocalMatrix<Curl_NED_3D, Curl_NED_3D> LocalC(mesh);
    

    LocalMatrix<RT0_2D, RT0_2D>           Local_Gamma(mesh2D);

    LocalVec<RT0_2D> Local0(mesh2D);

    LocalVec<NED_3D> Local00(mesh);

    LocalVec<Curl_NED_3D> Local11(mesh);
    
    
    Dof<NED_3D> dof(mesh);
    Dof<RT0_2D> dofRT(mesh2D);
    
    
    int nb_dof   = NbDof(dof);

    std::cout << "Number of edges: \t" << nb_dof << std::endl;

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
    



    
    
    
    
    // -----------------------------------------------------------------------//
    
    progress bar00("BEM Matrix Assembly\t",nb_close2D/NUM);
    
    std::vector<std::size_t>  rowBEM;
    std::vector<std::size_t>  colBEM;
    
    std::vector<Cplx>         valVk;
    std::vector<Cplx>         valKk;
    

    
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      
    std::vector<std::size_t>  rowBEM_priv;
    std::vector<std::size_t>  colBEM_priv;
    
    std::vector<Cplx>         valVk_priv;
    std::vector<Cplx>         valKk_priv;
    
      
      
      C3x3 temp1,temp2;

      #pragma omp for firstprivate(V0,K0,dofRT, Close_Gamma)
      for(int l=0; l<nb_close2D; l++){
      
         int j = Close_Gamma[l][0];
         int k = Close_Gamma[l][1];
             
         temp1 = V0(j, k);
         temp2 = K0(j, k);
          
          
         for(int nx = 0; nx < nb_dof_loc2D; nx++){
           for(int ny = 0; ny < nb_dof_loc2D; ny++){
            
             rowBEM_priv.push_back(dofRT[j][nx]);
             colBEM_priv.push_back(dofRT[k][ny]);
             
             valVk_priv.push_back(temp1(nx, ny));
             valKk_priv.push_back(temp2(nx, ny));
            
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
      }
      
      
      if(tid == 0){bar00.end();}
    
    }
    
    
    smatrix<Cplx> V_sparse(nb_dof2D,nb_dof2D, rowBEM, colBEM, valVk);
    smatrix<Cplx> K_sparse(nb_dof2D,nb_dof2D, rowBEM, colBEM, valKk);
    
    
    
    
    
    
    // -----------------------------------------------------------------------//
    // --------------- END NEAR FIELD MATRIX ---------------------------------//
    // -----------------------------------------------------------------------//
    
    matrix<double> Xb(nb_dof2D, 3);
    matrix<double> Yb(nb_dof2D, 3);
    
    for(int j = 0; j< nb_dof2D; j++){
    
    Xb(j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Xb(j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Xb(j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];
    // disp(eval(Xb(j, colon(0, 2))));
    
    Yb(j, 0) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][0];
    Yb(j, 1) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][1];
    Yb(j, 2) = dofRT(((dofRT.ToElt(j))[0])[0])[((dofRT.ToElt(j))[0])[1]][2];
    
    }

    // disp(Xb);

    // disp(" ");

    disp(" ");
    disp(" ");
    disp(" ");
    disp(" ");

    matrix<double> X(nb_dof, 3);
    matrix<double> Y(nb_dof, 3);
    
    for(int j = 0; j< nb_dof; j++){
    
    X(j, 0) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][0];
    X(j, 1) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][1];
    X(j, 2) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][2];
    
    // disp(eval(X(j, colon(0, 2))));
    Y(j, 0) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][0];
    Y(j, 1) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][1];
    Y(j, 2) = dof(((dof.ToElt(j))[0])[0])[((dof.ToElt(j))[0])[1]][2];
    
    }

    // disp(X);
    
    
    // -----------------------------------------------------------------------//
    // --------------- MASS MATRICES -----------------------------------------//
    // -----------------------------------------------------------------------//
    std::vector<std::size_t> row, col, rowJE, colJE; 
    
    std::vector<Cplx> val_eps, val, valC, valJE, valX, val_mu; 
    
    for(int j=0; j<nb_elt; j++){
    
      C6x6 temp = Local(j);
      C6x6 tempC = LocalC(j);

      C6x6 tempX = LocalX(j);
      
      C6x6 temp_eps = Local(epsilon, j);

      C6x6 temp_mu = LocalC(invmu, j);
    
      //std::cout << Local(j) << std::endl;
      for(int nx = 0; nx < nb_dof_loc; nx++){
        for(int ny = 0; ny < nb_dof_loc; ny++){
          row.push_back(dof[j][nx]);
          col.push_back(dof[j][ny]);
          
          val.push_back(temp(nx, ny));
          valC.push_back(tempC(nx, ny));
          val_eps.push_back(temp_eps(nx, ny));   
          val_mu.push_back(temp_mu(nx, ny));   
          valX.push_back(tempX(nx, ny));   
        }
      }
    }
    
    // disp("Check");
    
    smatrix<Cplx> M0(nb_dof,nb_dof, row, col, val);

    smatrix<Cplx> Mx(nb_dof,nb_dof, row, col, valX);
    // disp("Check");
    smatrix<Cplx> M1(nb_dof,nb_dof, row, col, valC);
    // disp("Check");

    smatrix<Cplx> Mcurl = M0 + M1;
    smatrix<Cplx> M0_eps(nb_dof,nb_dof, row, col, val_eps);

    smatrix<Cplx> M1_mu(nb_dof,nb_dof, row, col, val_mu);

    M1_mu -= M0_eps;



     std::vector<std::size_t> rowg, colg; 
    
    std::vector<Cplx> valg, valxg; 
    
    for(int j=0; j<nb_elt2D; j++){
    
      C3x3 tempg = Local_Gamma(j);

      C3x3 tempxg = Local_Gamma(j,"x");

      for(int nx = 0; nx < nb_dof_loc2D; nx++){
        for(int ny = 0; ny < nb_dof_loc2D; ny++){
          rowg.push_back(dofRT[j][nx]);
          colg.push_back(dofRT[j][ny]);
          
          valg.push_back(tempg(nx, ny)); 
          
          valxg.push_back(tempxg(nx, ny)); 

          
        }
      }
    }
    
    
    
    smatrix<Cplx> M0_Gamma(nb_dof2D,nb_dof2D, rowg, colg, valg);
    disp("Check");
    hmatrix<Cplx> M0h_Gamma(Xb, Yb, tol, M0_Gamma);



    smatrix<Cplx> JE(nb_dof2D,nb_dof2D, rowg, colg, valxg);
    matrix<Cplx> valuesJE = values(JE);
    // for(int n = 0; n<numel(valuesJE); n++){

    //   disp(valuesJE(n));

    // }

    // disp("Check");
    // disp(JE);

    disp(" ");

    // disp(M0_Gamma);


    hmatrix<Cplx> JEh(Xb, Yb, tol, JE);

    // -----------------------------------------------------------------------//
    // ------------------ Boundary-to-Volume DOFs-----------------------------//
    // -----------------------------------------------------------------------//
    std::vector<std::size_t> Pi, Pj;
    std::size_t countt = 0;

    for(int ii = 0; ii < nb_dof2D; ii++){
    
      for(int jj = 0; jj < nb_dof; jj++){

        if(norm(eval(Xb(ii, colon(0, 2))) - eval(X(jj, colon(0, 2))), "2") < 1e-8){

          Pi.push_back(ii);

          Pj.push_back(jj);

          countt++;
        }
      }
    }

    std::vector<Cplx> vec_P(nb_dof2D, 1.);


    // disp(X);
    // disp("\n\n");
    // disp(Xb);

    disp(Pi.size());
    disp(Pj.size());
    disp(vec_P.size());
    smatrix<Cplx> P_Gamma(nb_dof2D, nb_dof, Pi, Pj, vec_P);


    
    
    // -----------------------------------------------------------------------//
    // --------------- END MASS MATRICES -------------------------------------//
    // -----------------------------------------------------------------------//
    
    
    
    // -----------------------------------------------------------------------//
    // --------------------- H-MATRIX BEM -------------------------------//
    // -----------------------------------------------------------------------//
    

    matrix<double> Xbb = vertcat(Xb, Xb);
    matrix<double> Ybb = vertcat(Yb, Yb);
    
    castor::tic();
    
    // -----------------------------------------------------------------------//
    auto fctV = [dofRT, &V_sparse](matrix<std::size_t> Ix, matrix<std::size_t> Iy){
        matrix<Cplx> M(numel(Ix),numel(Iy));
        matrix<Cplx> V(numel(Ix),numel(Iy));
        std::size_t ix, iy;
        
        BIOp<EFIE_RT0xRT0> K0(MeshOf(dofRT),MeshOf(dofRT),kappa0); 
       
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = K0(dofRT.ToElt(ix),dofRT.ToElt(iy));
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
       
        for (std::size_t i=0; i<numel(Ix); ++i){
            ix = Ix(i);
            for (std::size_t j=0; j<numel(Iy); ++j){
                iy  = Iy(j);
                M(i,j) = K0(dofRT.ToElt(ix),dofRT.ToElt(iy));
            }
        }
        
        V = full(eval(K_sparse(Ix, Iy)));
        M += V;
        return M;};
        
        
    // -----------------------------------------------------------------------//
    
    
    // hmatrix<Cplx> Vop(Xb,Yb,(-1.) * tol,fctV);disp(Vop); //Vop *= (iu * kappa0);
    // hmatrix<Cplx> Kop(Xb,Yb,(-1.) * tol,fctK);disp(Kop);
    
    // hmatrix<Cplx> Vop(Xb,Yb,tol,fctV);disp(Vop); //Vop *= (iu * kappa0);
    // hmatrix<Cplx> Kop(Xb,Yb,tol,fctK);disp(Kop);

    //hmatrix<Cplx> BEMh(Ch, Dh, Ah, Bh);


    
    castor::toc();
    // -----------------------------------------------------------------------//
    // --------------------- END H-MATRIX ------------------------------------//
    // -----------------------------------------------------------------------//
    
    
    
    
    
    // -----------------------------------------------------------------------//
    
    std::cout << "\n\n" << std::endl;    
    


    hmatrix<Cplx> M0h(X,Y,tol,M0);
    hmatrix<Cplx> M1h(X,Y,tol,M1);


    // disp(M0);

    // disp("\n\n\n");

    // disp(M0_Gamma);

    // disp("\n\n\n");

    // disp(P_Gamma);

    // disp("\n\n\n");

    // disp(JE);

    // disp("\n\n\n");
    // double diffMatrix = norm(full(M0h_Gamma) - mtimes(mtimes(full(P_Gamma), full(M0h)), transpose(P_Gamma)), "2");

    // disp(diffMatrix);
    // disp("\n\n\n");
    // disp(M0_Gamma);

    // disp("\n\n\n");
    // disp(transpose(mtimes(mtimes(P_Gamma, M0), transpose(P_Gamma))));

    // disp("\n\n\n");
    // disp(M0);

    // hmatrix<Cplx> Mh = M0h + M1h; 

    //hmatrix<Cplx> MatFEM1(X,Y,tol,M1);

    //hmatrix<Cplx> MatFEM0(X,Y,tol,M0_eps); 
    
    hmatrix<Cplx> MatFEM(X,Y,1e-6,M1_mu); 

    // std::cout << "Computing inverse ..." << std::endl;

    castor::tic();


    // hmatrix<Cplx> M0hm1 = inv(M0h);

    // hmatrix<Cplx> FEMm1 = inv(MatFEM);
    
    // hmatrix<Cplx> Vopm1 = inv(Vop);

    std::function<matrix<Cplx>(matrix<Cplx>const&)> Schurfct, Am1fct, sparsePrecond, sparseM0;

    // Am1fct = [&Vopm1, &nb_dof, &nb_dof2D](matrix<Cplx> const B) {
        
    //     matrix<Cplx> X = mtimes(Vopm1, B);

    //     return X;

    // };

    // sparsePrecond = [&FEMm1](matrix<Cplx> const B) {
    // // sparsePrecond = [](matrix<Cplx> const B) {
        
    //     matrix<Cplx> X = mtimes(FEMm1, B);

    //     // matrix<Cplx> X = B;
    //     return X;

    // };


    // sparseM0 = [&M0hm1](matrix<Cplx> const B) {
    // // sparsePrecond = [](matrix<Cplx> const B) {
        
    //     matrix<Cplx> X = mtimes(M0hm1, B);

    //     // matrix<Cplx> X = B;
    //     return X;

    // };

    // Cplx half = 0.5;
    // hmatrix<Cplx> Ch  =  half*JEh;
    //               Ch +=  (-2*half) * Kop;

    // hmatrix<Cplx> Bh = -transpose(JEh);


    // Schurfct = [&Vop, &M1_mu, &sparsePrecond, &P_Gamma, &Bh, &Ch, &nb_dof, &nb_dof2D](matrix<Cplx> const B) {
        
    //     matrix<Cplx> Xa = mtimes(Vop, B);

    //     matrix<Cplx> Xb = mtimes(Bh, B);
    //                  Xb = mtimes(transpose(P_Gamma), Xb);

    //     if(norm(Xb, "2") > 1e-12){
    //     matrix<Cplx> Xbb= gmres(M1_mu, Xb, 1e-6, 2000, sparsePrecond, matrix<Cplx>(), 0);
    //                  //Xb = mtimes(FEMm1, Xb);


    //                  Xbb = mtimes(P_Gamma, Xbb);
    //                  Xbb = mtimes(Ch, Xbb);
    //                  Xa -= Xbb;
    //     }


    //     return Xa;

    // };




    // disp(inv(Matrix));
    // --------------------- END H-MATRIX ------------------------------------//
    // -----------------------------------------------------------------------//
    
    // matrix<Cplx> b0(nb_dof2D + nb_dof, 1);

    matrix<Cplx> Neu(nb_dof2D, 1);

    matrix<Cplx> Dir(nb_dof2D, 1);
    for(int j=0; j<nb_elt2D; j++){
      

      C3 vect0 = Local0(Hx, Hy, Hz, j); 
      C3 vect1 = Local0(Ex, Ey, Ez, j);
      
      
      for(int kj = 0; kj < nb_dof_loc2D; kj++){
      
        Dir(dofRT[j][kj])   += vect1[kj];  

        Neu(dofRT[j][kj])   += vect0[kj];  
      }
    }
    
    matrix<Cplx> UDir(nb_dof, 1);

    matrix<Cplx> UNeu(nb_dof, 1);

    // disp("Maybe here?");

    for(int j=0; j<nb_elt; j++){
      

      C6 vect1 = Local00(Ex, Ey, Ez, j);

      C6 vect0 = Local11(Hx, Hy, Hz, j); vect0 = iu * kappa0 * vect0;
      
      
      for(int kj = 0; kj < nb_dof_loc; kj++){
      
        UDir(dof[j][kj])   += vect1[kj];  

        UNeu(dof[j][kj])   += vect0[kj];  
 
      }
    }




    // matrix<Cplx> rhsDir_pre = mtimes(inv(M0h_Gamma), Dir);
    // matrix<Cplx> rhsDir  = mtimes(M0h_Gamma, rhsDir_pre); 
                //  rhsDir  = 0.5 * rhsDir;
                //  rhsDir += mtimes(Kop, rhsDir); 

    // matrix<Cplx> rhsNeu = mtimes(transpose(P_Gamma), Neu);
    //              rhsNeu = mtimes(MatFEM, rhsNeu);
    // matrix<Cplx> rhsNeuB = mtimes(P_Gamma, rhsNeu);
    //              rhsNeuB = mtimes(Ch, rhsNeuB);

    
    // rhsDir -= rhsNeuB; 
    
    // matrix<Cplx> solPre = gmres(Schurfct,Dir,1e-6,500,Am1fct,matrix<Cplx>(), 1);






    // solPre = mtimes(Bh, solPre);
    // solPre = mtimes(transpose(P_Gamma), solPre);
    // sol = mtimes(FEMm1, sol); 
    // matrix<Cplx> sol= gmres(M1_mu, solPre, 1e-6, 2000, sparsePrecond, matrix<Cplx>(), 1);

    // sol = real(sol);


    // castor::toc();
    // std::cout << "L2 Norm of solution: \t" << mtimes(transpose(sol), mtimes(M0, sol)) << std::endl;

    // std::cout << "Hcurl Norm of solution: \t" << mtimes(transpose(sol), mtimes(Mcurl, sol)) << std::endl;

    double mesh_size = 0.;

    for(int n = 0; n < nb_elt; n++){

      mesh_size = std::max(mesh_size, std::pow(std::abs(Vol(mesh[n])), 1./3.)); 

    }
    std::cout << "Mesh Size: \t" << mesh_size << std::endl;


    // disp("Maybe here?");
    matrix<Cplx> solref = mtimes(inv(full(M0)), UDir); solref = real(solref);

    // matrix<Cplx> curlRef = mtimes(inv(full(M0)), mtimes(Mx, solref));

    // matrix<Cplx> curlU = gmres(M0, mtimes(Mx, sol), 1e-6, 2000, sparseM0, matrix<Cplx>(), 1);curlU = real(curlU);

    // // matrix<Cplx> curlU = mtimes(inv(full(M0)), mtimes(Mx, sol)); curlU = real(curlU);

    // std::cout << "L2 Norm of curl(U): \t" << mtimes(transpose(curlU), mtimes(M0, curlU)) << std::endl;

    // std::cout << "Hcurl Norm of curl(U): \t" << mtimes(transpose(curlU), mtimes(Mcurl, curlU)) << std::endl;

    auto L2norm = mtimes(transpose(solref), mtimes(M0, solref));
    auto Hcurlnorm = L2norm + mtimes(transpose(solref), mtimes(M1, solref)); 
    std::cout << "L2 Norm of U: \t" << std::real(std::sqrt(L2norm(0))) << std::endl;

    std::cout << "Hcurl Norm of U: \t" << std::real(std::sqrt(Hcurlnorm(0))) << std::endl;


    // matrix<Cplx> diff = solref - sol;

    // std::cout << "L2 Error: \t" << mtimes(transpose(diff), mtimes(M0, diff)) << std::endl;

    // std::cout << "Hcurl Error: \t" << mtimes(transpose(diff), mtimes(Mcurl, diff)) << std::endl;
    

    // sol += rhsNeu;


    // matrix<Cplx> P = full(P_Gamma);

    // matrix<Cplx> Bfull(nb_dof2D, nb_dof2D);
    // Bh.hfull(Bfull);
    
    // matrix<Cplx> Bf = mtimes(transpose(P), Bfull);

    // matrix<Cplx> Cfull(nb_dof2D, nb_dof2D);
    // Ch.hfull(Cfull);

    // matrix<Cplx> Cf = mtimes(Cfull, P);

    // matrix<Cplx> Dfull(nb_dof2D, nb_dof2D);
    // Vop.hfull(Dfull);

    // Dfull *= (iu*kappa0);


    // matrix<Cplx> Matrix0 = horzcat((1./iu /kappa0) * full(MatFEM), Bf);

    // matrix<Cplx> Matrix0 = horzcat(full(MatFEM), Bf);
    // matrix<Cplx> Matrix1 = horzcat(Cf, Dfull);

    // matrix<Cplx> Matrix = vertcat(Matrix0, Matrix1);



    // matrix<Cplx> sol2 = mtimes(inv(Matrix), vertcat(zeros(nb_dof, 1), Dir));

    // matrix<Cplx> solDirect = eval(sol2(range(0,nb_dof), 0)); 

    // disp(norm(sol - solDirect, "2"));



    

    std::vector<double> solU(nb_dof, 0.);
    
    // std::ofstream fileU("Fichera/solutions2/FEM-BEM-U-" + std::to_string(ntest) + ".txt");
    
    // for(int n=0; n<nb_dof; n++){
      
    //   solU[n] = sol(n).real();
    
    //   if (fileU.is_open()){
    //     fileU << solU[n] << '\n';
    //   }
      
    // }
    
    // fileU.close();
    

  // // -------------------------------------------------------------------------------------------- //

  //   double norms    = 0.;
  //   double ref_norm = 0.;
  //   double error    = 0.;

  //   double norms0    = 0.;
  //   double ref_norm0 = 0.;
  //   double error0    = 0.;
    
  //   matrix<Cplx> diff(nb_dof, 1);
  //   for(int n = 0; n < nb_dof; n++){
  //     diff(n) = solU[n];
  //   }
    
    
    
    
  //   auto temp0 = mtimes(T, diff);
    
  //   auto temp00 = mtimes(M0, diff);
    
    
  //   for(int n = 0; n < nb_dof; n++){
  //     // ref_norm += std::abs(sol_reference(n) * temp1(n));
  //     // norm += std::abs(sol_interp[n] * temp1(n));
  //     error += std::abs(diff(n) * temp0(n));

  //     // ref_norm0 += std::abs(sol_reference(n) * temp01(n));
  //     // norms0 += std::abs(sol_interp[n] * temp01(n));
  //     error0 += std::abs(diff(n) * temp00(n));
  //   }
    
  //   // std::cout << std::sqrt(error0) / std::sqrt(ref_norm0) << "\t" <<  std::sqrt(error) / std::sqrt(ref_norm) << std::endl;
  //   std::cout << error0 << "\t" <<  error << std::endl;
  //   // std::cout << ref_norm0 << "\t" <<  ref_norm << std::endl;
  
  
  // -------------------------------------------------------------------------------------------- //






    std::vector<double> solUref(nb_dof, 0.);
    
    std::ofstream fileUref("solutions0/FEM-BEM-Uinc-" + std::to_string(ntest) + ".txt");
    
    for(int n=0; n<nb_dof; n++){
      
      solUref[n] = solref(n).real();
     
      if (fileUref.is_open()){
        fileUref << solUref[n] << '\n';
      }
      
    }
    
    fileUref.close();
    

    std::vector<double> sol_curlU(nb_dof, 0.);
    
       
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
          }
        }
        
      }
    }


    // WriteEltVectGmsh(dof, mesh2D, boundary, "Fichera/output2/valsFEM-BEM-U-" + std::to_string(ntest) + ".txt", solU);
    // WriteEltVectGmsh(dof, mesh2D, boundary, "Fichera/output2/valsFEM-BEM-curlU-" + std::to_string(ntest) + ".txt", sol_curlU);
    // WriteEltVectGmsh(dof, mesh2D, boundary, "output0/valsFEM-BEM-Uinc-" + std::to_string(ntest) + ".txt", solUref);
    // WriteEltVectGmsh(dof, mesh2D, boundary, "output0/valsFEM-BEM-curlUref-" + std::to_string(ntest) + ".txt", solCurlref);
   
    
   
    std::cout << "---------------------------------------" << std::endl;
  }
  
  
  
  
  
  
  
  static void Errors(std::string filename, std::string extension, const int ntest, const int nmax){
  
    
    Geometry node0(filename + std::to_string(ntest) + extension);
    Mesh3D mesh0; mesh0.Load(node0,-1);
    
    Geometry node(filename + std::to_string(nmax) + extension);
    Mesh3D mesh; mesh.Load(node,-1);
    
    
    
    int nb_elt  = NbElt(mesh);
    int nb_elt0 = NbElt(mesh0);
    

    Dof<NED_3D> dof(mesh);
    Dof<NED_3D> dof0(mesh0); NED_3D phi(mesh0);

    Dof<P1_3D> dof00(mesh0);
    
    
    int nb_dof = NbDof(dof);
    int nb_dof0 = NbDof(dof0);

    int nb_dof00 = NbDof(dof00);
    int nb_dof_loc = NED_3D::nb_dof_loc;
    


    std::cout << "Number of elements: \t" << nb_elt0 << std::endl;
    std::cout << "Number of edges: \t" << nb_dof0 << std::endl;
    std::cout << "Number of nodes: \t" << nb_dof00 << std::endl;

    std::vector<double>  sol_real(nb_dof0, 0.);
    std::vector<double>  curl_sol_real(nb_dof0, 0.);
    
    std::vector<double>  sol_interp(nb_dof, 0.);
    matrix<Cplx>         sol_reference(nb_dof, 1);


    std::vector<double>  curl_sol_interp(nb_dof, 0.);
    matrix<Cplx>         curl_sol_reference(nb_dof, 1);
    
    // ---------------------------------------------------------------//
    // ----------- LOAD SOLUTIONS INTO VECTORS -----------------------//
    // ---------------------------------------------------------------//



    double value;

    
    // std::ifstream file0("Fichera/solutions2/FEM-BEM-U-" + std::to_string(ntest) + ".txt");

    // std::ifstream file0("Tetra/solutions2/FEM-BEM-U-" + std::to_string(ntest) + ".txt");
    std::ifstream file0("Tetra/solutions2/STF-VIE-U-" + std::to_string(ntest) + ".txt");

    // std::ifstream file0("Fichera/solutions2/STF-VIE-U-" + std::to_string(ntest) + ".txt");
    
    omp_set_num_threads(1);
    
    // progress barSol("Load solution\t",nb_dof0);
    for(int n = 0; n < nb_dof0; n++){
      file0 >> value;
      sol_real[n] = value;

    }
    
    file0.close();

    // std::ifstream file0c("Fichera/solutions2/FEM-BEM-curlU-" + std::to_string(ntest) + ".txt");

    // std::ifstream file0c("Tetra/solutions2/FEM-BEM-curlU-" + std::to_string(ntest) + ".txt");
    std::ifstream file0c("Tetra/solutions2/STF-VIE-curlU-" + std::to_string(ntest) + ".txt");
    
    // std::ifstream file0c("Fichera/solutions2/STF-VIE-curlU-" + std::to_string(ntest) + ".txt");
    
    for(int n = 0; n < nb_dof0; n++){
      file0c >> value;
      curl_sol_real[n] = value;

    }
    
    file0c.close();
    
    
    
    // std::ifstream file("Fichera/solutions2/FEM-BEM-U-" + std::to_string(nmax) + ".txt");
    // std::ifstream file("Fichera/solutions2/STF-VIE-U-" + std::to_string(nmax) + ".txt");

    std::ifstream file("Tetra/solutions2/STF-VIE-U-" + std::to_string(nmax) + ".txt");
    // std::ifstream file("solutions8/STF-VIE-U-" + std::to_string(nmax) + ".txt");

    
    for(int n = 0; n < nb_dof; n++){
      file >> value;
      sol_reference(n) = value;

    }
    
    file.close();
    



    // std::ifstream filec("Fichera/solutions2/FEM-BEM-curlU-" + std::to_string(nmax) + ".txt");
    // std::ifstream filec("Fichera/solutions2/STF-VIE-curlU-" + std::to_string(nmax) + ".txt");

    // std::ifstream filec("Tetra/solutions2/STF-VIE-curlU-" + std::to_string(nmax) + ".txt");

    std::ifstream filec("Tetra/solutions2/STF-VIE-curlU-" + std::to_string(nmax) + ".txt");
    for(int n = 0; n < nb_dof; n++){
      filec >> value;
      curl_sol_reference(n) = value;

    }
    
    filec.close();
    
    
    // ---------------------------------------------------------------//
    // ----------- FINITE ELEMENT MATRIX ASSEMBLY --------------------//
    // ---------------------------------------------------------------//
    
    LocalMatrix<NED_3D, NED_3D>           Local(mesh);
    LocalMatrix<Curl_NED_3D, Curl_NED_3D> LocalC(mesh);
    
    // ---------------------------------------------------------------//
    // ---------------------------------------------------------------//
    
    
    std::vector<std::size_t> row0, col0, row1, col1; 
    
    std::vector<Cplx> val, valC; 
    
    omp_set_num_threads(1);
    
    // progress barMass("Mass Matrix\t",nb_elt);
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
      // barMass++;
      
    }
    
    
    
    smatrix<Cplx> M0(nb_dof,nb_dof, row0, col0, val);
    smatrix<Cplx> M1(nb_dof,nb_dof, row1, col1, valC);
    
    smatrix<Cplx> T = M0 + M1;
    
    // barMass.end();
    
    
     
    
    std::vector< std::vector<int> > C2F; 
    
    
    for(int nelt = 0; nelt < nb_elt0; nelt++){
      
      // std::cout << "\r";
      // std::cout << "H-Matrix Horzcat Assembly" << ": \t";
      // std::cout << l + 1 << " / " << ptr.size();
      // std::cout.flush();
    
      int count = 0;
    
      Elt3D elt0 = mesh0[nelt];
      
      auto faces = FacesOf(elt0);
      
      // std::vector<int> c2f(std::pow(8, (nmax - ntest)), 0);
      std::vector<int> c2f;

      // std::vector<int> c2f(8, 0);
      
      
      
      for(int n = 0; n < nb_elt; n++){
        
        Elt3D elt = mesh[n];
        
        R3 C1 = Ctr(elt);

        bool inside = 1;      
        for(int f = 0; f < 4; f++){
          R3 normal = NormalTo(faces[f]);
          inside *= ((normal, C1 - faces[f][0]) * (normal, elt0[f] - faces[f][0]) > 0. );
        }
        
        if(inside){
          // c2f[count] = n;

          c2f.push_back(n);
          
          //std::cout << count << std::endl;
          count+=1;
        }
      
      }
      
      C2F.push_back(c2f);
      
      phi.Assign(nelt);
      
      R3x3 B = MatJac(elt0);
      
      // for(int jj = 0; jj < std::pow(8, (nmax - ntest)); jj++){
      for(int jj = 0; jj < C2F[nelt].size(); jj++){

      // for(int jj = 0; jj < 8; jj++){
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



    std::vector< std::vector<int> > C2Fb; 
   
    
    for(int nelt = 0; nelt < nb_elt0; nelt++){
      
      // std::cout << "\r";
      // std::cout << "H-Matrix Horzcat Assembly" << ": \t";
      // std::cout << l + 1 << " / " << ptr.size();
      // std::cout.flush();
    
      int count = 0;
    
      Elt3D elt0 = mesh0[nelt];
      
      auto faces = FacesOf(elt0);
      
      // std::vector<int> c2f(std::pow(8, (nmax - ntest)), 0);
      std::vector<int> c2f;

      // std::vector<int> c2f(8, 0);
      
      
      
      for(int n = 0; n < nb_elt; n++){
        
        Elt3D elt = mesh[n];
        
        R3 C1 = Ctr(elt);

        bool inside = 1;      
        for(int f = 0; f < 4; f++){
          R3 normal = NormalTo(faces[f]);
          inside *= ((normal, C1 - faces[f][0]) * (normal, elt0[f] - faces[f][0]) > 0. );
        }
        
        if(inside){
          // c2f[count] = n;
          c2f.push_back(n);
          
          //std::cout << count << std::endl;
          count+=1;
        }
      
      }
      
      C2Fb.push_back(c2f);
      
      phi.Assign(nelt);
      
      R3x3 B = MatJac(elt0);
      
      // for(int jj = 0; jj < std::pow(8, (nmax - ntest)); jj++){
      for(int jj = 0; jj < C2Fb[nelt].size(); jj++){

      // for(int jj = 0; jj < 8; jj++){
      int idn = C2Fb[nelt][jj];
      
      Elt3D elt = mesh[idn];
      
      auto edges = EdgesOf(elt);
      
      for(int ll = 0; ll < nb_dof_loc; ll++){
      
      Elt1D edge_ll = edges[ll]; Order(edge_ll);
      
      R3 tangent = edge_ll[1] - edge_ll[0];         // ORIENTATION DEFINED THIS WAY --> ref. shapefct.hpp
      
      curl_sol_interp[dof[idn][ll]] = 0.;
      
      R3   x = dof(idn)[ll] - elt0[0];
      
      R3   t = inv(B) * x;
      
      for(int ii = 0; ii < nb_dof_loc; ii++){
      curl_sol_interp[dof[idn][ll]] += (curl_sol_real[dof0[nelt][ii]] * phi(ii, t), tangent);
      }
      
      
      }
      
      
      }
      
    }
    
    // -----------------------------------------------------------------------//

    // double ref_norm = 0.;
    // double error    = 0.;

    
    // matrix<Cplx> diff(nb_dof, 1);
    // for(int n = 0; n < nb_dof; n++){
    //   diff(n) = sol_interp[n] - sol_reference(n);
    // }

    
    matrix<Cplx> diff(nb_dof, 1);
    matrix<Cplx> curl_diff(nb_dof, 1);

    for(int n = 0; n < nb_dof; n++){
      diff(n) = sol_interp[n] - sol_reference(n);

      curl_diff(n) = curl_sol_interp[n] - curl_sol_reference(n);
    }
    
    
    // std::cout << std::sqrt(error0) / std::sqrt(ref_norm0) << "\t" <<  std::sqrt(error) / std::sqrt(ref_norm) << std::endl;
   auto error = mtimes(transpose(diff), mtimes(M0, diff));
   auto ref_norm = mtimes(transpose(sol_reference), mtimes(M0, sol_reference));

   auto error2 = mtimes(transpose(curl_diff), mtimes(M0, curl_diff));
   auto ref_norm2 = mtimes(transpose(curl_sol_reference), mtimes(M0, curl_sol_reference));



  //  std::cout << "L2 error: \t" << error << std::endl;
  //  std::cout << "L2 norm: \t" << ref_norm << std::endl;
  //  std::cout << "L2 relative error: \t" << error / ref_norm << std::endl;

    // std::cout << error0 << "\t" <<  error << std::endl;

    // std::cout << ref_norm0 << "\t" <<  ref_norm << std::endl;
  
 
    double mesh_size = 0.;

    for(int n = 0; n < nb_elt0; n++){

      mesh_size = std::max(mesh_size, std::pow(std::abs(Vol(mesh0[n])), 1./3.)); 

    }
    // std::cout << "Mesh Size: \t" << mesh_size << std::endl;




std::cout <<  std::real(std::sqrt(error(0))) << "\t" << std::real(std::sqrt(ref_norm(0))) << "\t" << std::real(std::sqrt(error(0)+error2(0))) << "\t" << std::real(std::sqrt(ref_norm(0)+ref_norm2(0))) << "\t"<< mesh_size << std::endl;



  // std::cout << "\n" << std::endl;
  
  }
  
  
  
  
  
  

};



int main(){

  std::string filename  = "mesh/oriented_cube";  
  // std::string filename  = "mesh/sphere";
  std::string extension = ".msh";
  
  const int nmax = 2;
  
  // Compute solutions 
  for(int n = 0; n < nmax+1; n++){
    Test::Launch(filename, extension, n);
  }
  
  
  // Compute errors
  // for(int n = 0; n < nmax+1; n++){
  //   Test::Errors(filename, extension, n, nmax);
  // }
  
  
}


