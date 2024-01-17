/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : hmatrix.hpp                                   |
 |    #    |   VERSION    : 1.0.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 14.02.2021                                    |
 | ( === ) |   SYNOPSIS   : Hierarchical matrix operations & algebra with |
 |  `---'  |                matlab-like standard functions and procedures |
 +========================================================================+
 */

#pragma once
#define CASTOR_HMATRIX_HPP
#include <castor/matrix.hpp>
#include <castor/linalg.hpp>
#include <castor/smatrix.hpp>
#include <memory>
#include <fstream>

namespace castor
{

//==========================================================================//
//                           BOUNDING BOX CLASS                             //
//==========================================================================//
// [bbox]
/// Bounding box
///
/// Bounding box data from n-d nodes X.
class bbox
{
public:
    // CONSTRUCTOR
    bbox(){};
    template<typename T>
    bbox(matrix<T>const& X);
    
    // FUNCTIONS
    bool isfar(bbox const& Yb, double eta=3.) const;
    
    // ACCESSORS
    matrix<double>const& ctr()  const {return m_ctr;};
    matrix<double>const& edg()  const {return m_edg;};
    matrix<double>const& maxi() const {return m_max;};
    matrix<double>const& mini() const {return m_min;};
    
private:
    matrix<double> m_ctr;    // CENTER
    matrix<double> m_edg;    // EDGE LENGTH
    matrix<double> m_max;    // MAXIMUM
    matrix<double> m_min;    // MINIMUM
};

//==========================================================================
//
/// Constructor of bounding box
template<typename T>
bbox::bbox(matrix<T>const& X)
{
    m_min = min(X,1);
    m_max = max(X,1);
    m_ctr = (m_min+m_max)/2.;
    m_edg = m_max-m_min;
}

//==========================================================================
//
/// Evaluate distance with another bounding box using admissibility
/// criterium fixed by eta (default is 3).
bool bbox::isfar(bbox const& Yb, double eta) const
{
    matrix<double> zero(1,3);
    matrix<double> d1 = maximum(zero,m_min-Yb.maxi());
    matrix<double> d2 = maximum(zero,Yb.mini()-m_max);
    double diamX  = std::sqrt(sum(pow(m_edg,2)));
    double diamY  = std::sqrt(sum(pow(Yb.edg(),2)));
    double distXY = std::sqrt(sum(pow(d1,2)+pow(d2,2)));
    return (std::min(diamX,diamY) <= eta*distXY);
}

//==========================================================================
// [operator<<]
/// Flux operator to display hmatrix size and statistics using standard
/// C++ convention of disp(Ah) operator.
std::ostream& operator<<(std::ostream& flux, bbox const& Xb)
{
    flux << "Bounding box data :" << std::endl;
    flux << "  Minimum : " << Xb.mini() << std::endl;
    flux << "  Center  : " << Xb.ctr()  <<std::endl;
    flux << "  Maximum : " << Xb.maxi() <<std::endl;
    flux << "  Edge    : " << Xb.edg();
    return flux;
}


//==========================================================================//
//                         BINARY TREE CLASS                                //
//==========================================================================//
// [bintree]
/// Binary tree
///
/// bintree(X,N) is a recursive tree with a binary structure and bounding boxes.
/// Each branch is divided into two new ones, with a median distribution of
/// the M tri-dimensional nodes X, stored in a M-by-3 matrix. Leaves are reached
/// when there are less than N nodes in a branch (default is N=100).
template<typename T>
class bintree
{
public:
    // CONSTRUCTOR
    bintree(){};
    bintree(matrix<T>const& X, std::size_t n=300);
    
    // FUNCTIONS
    matrix<T> leaf()  const;
    
    // ACCESSORS
    bool const&               isleaf()   const {return m_isleaf;};
    bbox const&               box()      const {return m_box;};
    matrix<T>const&           crd()      const {return m_crd;};
    matrix<std::size_t>const& ind(int i) const {return m_ind[i];};
    bintree<T>const&          sub(int i) const {return m_sub[i];};
    
private:
    bool                    m_isleaf; // TRUE IF LEAF
    bbox                    m_box;    // BOUNDING BOX DATA
    matrix<T>               m_crd;    // POINTS COORDINATES
    matrix<std::size_t>     m_ind[2]; // CHILDREN INDICES
    std::vector<bintree<T>> m_sub;    // CHILDREN (RECURSION)
};

//==========================================================================
//
/// Constructor of binary tree. Nodes coordinates are stored in a M-by-3 matrix
/// and maximum leaf size is fixed by n (default is 100).
template<typename T>
bintree<T>::bintree(matrix<T>const& X, std::size_t n)
{
    // Initialize
    std::size_t N = size(X,1);
    m_isleaf = (N<n);
    m_box    = bbox(X);
    m_crd    = X;
    
    // Children recursion
    if (!m_isleaf)
    {
        std::size_t         dim  = argmax(m_box.edg());
        matrix<std::size_t> iloc = argsort(eval(X(range(0,N),dim)));
        m_ind[0] = eval(iloc(range(0,N/2)));
        m_ind[1] = eval(iloc(range(N/2,N)));
        m_sub.push_back(bintree(eval(X(m_ind[0],col(X))),n));
        m_sub.push_back(bintree(eval(X(m_ind[1],col(X))),n));
    }
}

//==========================================================================
//
/// .leaf() gives a vector of random values, same for each leaf. Usefull to
/// visualize tree leaves.
template<typename T>
matrix<T> bintree<T>::leaf() const
{
    matrix<T> v(1,size(m_crd,1));
    if (m_sub.size()==2)
    {
        v(m_ind[0]) = m_sub[0].leaf();
        v(m_ind[1]) = m_sub[1].leaf();
    }
    else
    {
        v = rand(1)*ones(size(v));
    }
    return v;
}


//==========================================================================//
//                           H-MATRIX CLASS                                 //
//==========================================================================//
// [hmatrix]
/// Hierarchical matrices.
///
/// Array with element of type T, constructed using X and Y nodes of type S,
/// stored in hierarchical format using rank default compression.
///
/// Values are stored in matrix<T>, associated to hierarchical indices in
/// matrix<std::size_t> of the same dimensions of the values. Recursion
/// is constructed using std::vector<hmatrix<T>> object.
///
/// hmatrix<T> can be constructed using dense matrix, sparse matrix or
/// lambda function defining block matrix values.
template<typename T>
class hmatrix
{
public:
    // CONSTRUCTOR
    hmatrix(){};
    hmatrix(std::size_t m, std::size_t n, double tol, T v=0);
    hmatrix(hmatrix<T>& A, hmatrix<T>& B, hmatrix<T>& C, hmatrix<T>& D); // NEW!
    hmatrix(matrix<T>const& A, matrix<T>const& B, double tol);
    template<typename S>
    hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol, matrix<T>const& M);
    template<typename S>
    hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol, smatrix<T>const& Ms);
    template<typename S>
    hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol,
            std::function<matrix<T>(matrix<std::size_t>,matrix<std::size_t>)>const& fct);
    template<typename S>
    hmatrix(matrix<S>const& X, matrix<S>const& Y, hmatrix<T>const& Ml, hmatrix<T>const& Mr);
    template<typename S>
    hmatrix(bintree<S>const& Xb, bintree<S>const& Yb, double tol);
    template<typename S>
    hmatrix(bintree<S>const& Xb, bintree<S>const& Yb, double tol, smatrix<T>const& Ms);
    template<typename S>
    hmatrix(bintree<S>const& Xb, bintree<S>const& Yb, smatrix<T>const& Ml, hmatrix<T>const& Mh, smatrix<T>const& Mr);
    
    // FUNCTIONS
    void full2lowrank();
    void fusion();
    void hfull(matrix<T>& M, matrix<std::size_t> I={}, matrix<std::size_t> J={}) const;
    void hfull2(matrix<T>& M, matrix<std::size_t> idx, matrix<std::size_t> jdx,
              matrix<std::size_t> I={}, matrix<std::size_t> J={}) const;
    void hinv();
    void hllowsolve(hmatrix<T>const& Lh);
    void hllowsolve(matrix<T>& B, matrix<std::size_t> I={}) const;
    void hlmtimes(matrix<T>const& B, matrix<T>& C, matrix<std::size_t> I={}, matrix<std::size_t> J={}) const;
    void hlu(hmatrix<T>& U);
    void hlupsolve(matrix<T>& B, matrix<std::size_t> I={}) const;
    void hrmtimes(matrix<T>const& A, matrix<T>& C, matrix<std::size_t> I={}, matrix<std::size_t> J={}) const;
    void hrupsolve(hmatrix<T>const& Uh);
    void hrupsolve(matrix<T>& B, matrix<std::size_t> J={}) const;
    void htranspose();
    void leafptr(std::vector<hmatrix<T>*>& ptr,
                     std::vector<matrix<std::size_t>>& idx, std::vector<matrix<std::size_t>>& jdx,
                     matrix<std::size_t> I={}, matrix<std::size_t> J={});
    void recompress(double tol);
    void clear();
    void spydata(smatrix<logical>& Sf, smatrix<logical>& Sc, matrix<std::size_t> I={}, matrix<std::size_t> J={}) const;
    void stat(matrix<std::size_t>& info) const;
    void tgeabhm(T alpha, matrix<T>const& A, matrix<T>const& B, T beta, matrix<std::size_t> I={}, matrix<std::size_t> J={});
    void tgehmhm(T alpha, hmatrix<T>const& Ah, hmatrix<T>const& Bh, T beta);
    
    // OPERATORS
    hmatrix<T>& operator+=(T b);
    hmatrix<T>& operator*=(T b);
    hmatrix<T>& operator+=(hmatrix<T>const& Bh);
    matrix<T>   operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J) const;
    
    // ACCESSORS
    hmatrix<T>const&          getChildren(int i) const {return m_chd[i];};
    matrix<std::size_t>const& getColumn(int i)   const {return m_col[i];};
    matrix<T>const&           getData(int i)     const {return m_dat[i];};
    matrix<std::size_t>const& getRow(int i)      const {return m_row[i];};
    std::size_t               getSize(int i)     const {return m_siz(i);};
    double                    getTolerance()     const {return m_tol;};
    int                       getType()          const {return m_typ;};
    std::vector<hmatrix<T>>&  setChildren()      {return m_chd;};
    hmatrix<T>&               setChildren(int i) {return m_chd[i];};
    matrix<std::size_t>&      setColumn(int i)   {return m_col[i];};
    matrix<T>&                setData(int i)     {return m_dat[i];};
    matrix<std::size_t>&      setRow(int i)      {return m_row[i];};
    matrix<std::size_t>&      setSize()          {return m_siz;};
    double&                   setTolerance()     {return m_tol;};
    int&                      setType()          {return m_typ;};

    //std::vector<std::shared_ptr<hmatrix<T>>> ptr;
    
private:
    int                      m_typ=-1;    // LEAF TYPE (-1:UNDEF; 0:RECURSIVE; 1:COMPRESSED; 2:FULL)
    matrix<std::size_t>      m_siz={0,0}; // LEAF SIZE
    double                   m_tol;       // ACCURACY OF COMPRESSION
    matrix<T>                m_dat[2];    // MATRIX DATA (M or A*B)
    matrix<std::size_t>      m_row[4];    // CHILDREN ROWS INDICES (M11,M12,M21,M22)
    matrix<std::size_t>      m_col[4];    // CHILDREN COLUMNS INDICES
    std::vector<hmatrix<T>>  m_chd;       // CHILDREN (RECURSION)
};

//==========================================================================//
//                              CONSTRUCTORS                                //
//==========================================================================//
//==========================================================================
// [.hmatrix]
/// Single value constructor, fill leaf with tensor product. No tree is needed.
template<typename T>
hmatrix<T>::hmatrix(std::size_t m, std::size_t n, double tol, T v)
{
    m_typ    = 1;
    m_siz    = {m,n};
    m_tol    = tol;
    m_dat[0] = matrix<T>(m,1,v);
    m_dat[1] = matrix<T>(1,n,1.);
}



//==========================================================================
// [.hmatrix]
/// Block H-Matrix Constructor. No tree is needed.
template<typename T>
hmatrix<T>::hmatrix(hmatrix<T>& A, hmatrix<T>& B, hmatrix<T>& C, hmatrix<T>& D)
{
    m_typ    = 0;
    m_siz    = {A.getSize(0) + D.getSize(0), A.getSize(1) + D.getSize(1)};
    m_tol    = A.getTolerance();
    m_chd.resize(4,hmatrix<T>());
    
    
    std::vector<hmatrix<T>*> ptrB;
    std::vector<matrix<std::size_t>> idxB;
    std::vector<matrix<std::size_t>> jdxB;
    B.leafptr(ptrB, idxB, jdxB);
    
    
    for(std::size_t n = 0; n < ptrB.size(); n++){
      jdxB[n] = jdxB[n] + A.getSize(1);
    }
    
    
    std::vector<hmatrix<T>*> ptrC;
    std::vector<matrix<std::size_t>> idxC;
    std::vector<matrix<std::size_t>> jdxC;
    C.leafptr(ptrC, idxC, jdxC);
    
    
    for(std::size_t n = 0; n < ptrC.size(); n++){
      idxC[n] = idxC[n] + A.getSize(0);
    }
    
    
    std::vector<hmatrix<T>*> ptrD;
    std::vector<matrix<std::size_t>> idxD;
    std::vector<matrix<std::size_t>> jdxD;
    D.leafptr(ptrD, idxD, jdxD);
    
    
    for(std::size_t n = 0; n < ptrD.size(); n++){
      idxD[n] = idxD[n] + A.getSize(0);
      jdxD[n] = jdxD[n] + A.getSize(1);
    }
    
    
    
    
    
    
    m_chd[0] = A;
    m_chd[1] = B;
    m_chd[2] = C;
    m_chd[3] = D;
    
    m_row[0] = range(0, A.getSize(0));                                    m_col[0] = range(0, A.getSize(1));
    m_row[1] = range(0, A.getSize(0));                                    m_col[1] = range(0, B.getSize(1)) + A.getSize(1);
    m_row[2] = range(A.getSize(0), A.getSize(0) + C.getSize(0));          m_col[2] = range(0, A.getSize(1));
    m_row[3] = range(A.getSize(0), A.getSize(0) + C.getSize(0));          m_col[3] = range(0, B.getSize(1)) + A.getSize(1);
    
   
    // std::cout << "End of Block constructor" << std::endl;

}

//==========================================================================
// [.hmatrix]
/// Low-rank constructor, fill leaf with tensor product. No tree is needed.
template<typename T>
hmatrix<T>::hmatrix(matrix<T>const& A, matrix<T>const& B, double tol)
{
    m_typ    = 1;
    m_siz    = {size(A,1),size(B,2)};
    m_tol    = tol;
    m_dat[0] = A;
    m_dat[1] = B;
}

//==========================================================================
// [.hmatrix]
/// Dense matrix contructor using nodes, fill close interactions with dense
/// matrix and far with compressed. Parallelized leaves constructor.
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol, matrix<T>const& M)
{
    auto fct = [&M](matrix<std::size_t> Ix, matrix<std::size_t> Iy)
    {
        return eval(M(Ix,Iy));
    };
    (*this) = hmatrix(X,Y,tol,fct);
}

//==========================================================================
// [.hmatrix]
/// Sparse matrix contructor using nodes, fill non-zeros leaves with full
/// matrix and empty leaves with a zero tensor product. Recursive constructor.
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol, smatrix<T>const& Ms)
{
    (*this) = hmatrix<T>(bintree<S>(X),bintree<S>(Y),tol,Ms);
}


// //==========================================================================
// // [.hmatrix]
// /// Horzcat constructor
// /// 
// template<typename T>
// template<typename S>
// hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, hmatrix<T>& Ml, hmatrix<T>& Mr, 
//                     matrix<std::size_t> I, matrix<std::size_t> J){
 
//     double tol = max(Ml.getTolerance(), Mr.getTolerance());

//     bintree<S> Xb = bintree<S>(X);
//     bintree<S> Yb = bintree<S>(Yb);

//     (*this) = hmatrix<T>(Xb, Yb, tol);



//     matrix<T> Al, Bl, Ar, Br;

//     if(Xb.box().isfar(Yb.box())){

//         Mleft  = Ml(Ix, Iyl);

//         Mright = Mr(Ix, Iyr);

//         M = horzcat(Mleft, Mright);

//         bool flag;
//         // Full
//         if (nnz(Mleft)>numel(Mleft)/4)
//         {
//             // std::cout << "Number of elements: " << m_dat[0].size(1) << ", " << m_dat[0].size(2) << std::endl; 
//             std::tie(A,B,flag) = aca(m_dat[0],m_tol/10);
//         }
//         // "Sparse"
//         else
//         {
//             std::tie(A,B) = qrsvd(m_dat[0],eye<T>(m_siz(1)),m_tol/10);
//             (size(A,2)<min(m_siz)/2)?(flag = true):(flag = false);
//         }
//         std::tie(A, B, flag) = lowrank()



//     }







// }


//==========================================================================
// [.hmatrix]
/// Lambda function contructor using nodes, fill close interactions with dense
/// matrix and far with compressed. Parallelized leaves constructor.
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, double tol,
                   std::function<matrix<T>(matrix<std::size_t>,matrix<std::size_t>)>const& fct)
{
   // Build tree and block interactions
//    std::shared_ptr<hmatrix<T>> shrdptr = std::make_shared<hmatrix<T>>(bintree<S>(X),bintree<S>(Y),tol);
//    (*this) = *shrdptr;
   (*this) = hmatrix<T>(bintree<S>(X),bintree<S>(Y),tol);

   // Get leaves references
   std::vector<hmatrix<T>*> ptr;
   std::vector<matrix<std::size_t>> idx;
   std::vector<matrix<std::size_t>> jdx;
   leafptr(ptr,idx,jdx);
   bool flag;

   // Build leaf data with partial pivoting
   #pragma omp parallel
   #pragma omp single
   for (std::size_t l=0; l<ptr.size(); ++l)
   {
      std::cout << "\r";
      std::cout << "H-Matrix Assembly" << ": \t";
      std::cout << l + 1 << " / " << ptr.size();
      std::cout.flush();
   
       if (ptr[l]->m_siz(0)!=numel(idx[l]) || ptr[l]->m_siz(1)!=numel(jdx[l]))
       {
           error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
       }
       if (ptr[l]->m_typ==1 && tol>0)
       {
           #pragma omp task
           {
           std::tie(ptr[l]->m_dat[0],ptr[l]->m_dat[1],flag) = aca(idx[l],jdx[l],fct,tol);
           if (!flag)
           {
               std::cout << "ACA compression failed: " << l << std::endl;
               //error(__FILE__, __LINE__, __FUNCTION__,"ACA compression failed.");
           }
           }
       }
       else if (ptr[l]->m_typ==2 || tol<=0)
       {
           #pragma omp task
           ptr[l]->m_dat[0] = fct(idx[l],jdx[l]);
           ptr[l]->m_typ=2;

       }
       else
       {
           error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
       }
       
   }
   
   std::cout << "\n" << std::endl;

//    std::cout << "We are ok"<< std::endl;
}


//==========================================================================
// [.hmatrix]
/// Horzcat Version 2
///
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(matrix<S>const& X, matrix<S>const& Y, hmatrix<T>const& Ml, hmatrix<T>const& Mr)
{
    double tol = std::max(Ml.getTolerance(), Mr.getTolerance());

   (*this) = hmatrix<T>(bintree<S>(X),bintree<S>(Y),tol);

   std::vector<hmatrix<T>*> ptr;
   std::vector<matrix<std::size_t>> idx;
   std::vector<matrix<std::size_t>> jdx;

   leafptr(ptr,idx,jdx);
   bool flag;

  
   for (std::size_t l=0; l<ptr.size(); ++l)
   {
      std::cout << "\r";
      std::cout << "H-Matrix Horzcat Assembly" << ": \t";
      std::cout << l + 1 << " / " << ptr.size();
      std::cout.flush();
   
       if (ptr[l]->m_siz(0)!=numel(idx[l]) || ptr[l]->m_siz(1)!=numel(jdx[l]))
       {
           error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
       }


       matrix<std::size_t> jdxl, Jlh;
       matrix<std::size_t> jdxr, Jrh;
      
       matrix<std::size_t> Jl = range(0, Ml.getSize(1));
       matrix<std::size_t> Jr = range(0, Mr.getSize(1)) + Ml.getSize(1);
       std::tie(Jlh, jdxl) = argintersect(Jl, jdx[l]);
       std::tie(Jrh, jdxr) = argintersect(Jr, jdx[l]);


       Jr -= Ml.getSize(1) * ones(size(Jr));
       matrix<T> A, B;

       matrix<T> Fl(numel(idx[l]), numel(Jlh));


       if(numel(Jlh) > 0){Fl = Ml(idx[l], eval(Jl(Jlh)));}

       matrix<T> Fr(numel(idx[l]), numel(Jrh));
       if(numel(Jrh)> 0){Fr = Mr(idx[l], eval(Jr(Jrh)));}

       matrix<T> F  = zeros(length(idx[l]), length(jdx[l]));

       std::size_t jl = numel(Jlh);
       std::size_t jr = numel(Jrh);


       for (std::size_t i=0; i<numel(idx[l]); ++i){
        
       if(jl > 0){ 
	   for (std::size_t j=0; j<jl; ++j){
                F(i,jdxl(j)) = Fl(i, j); }
       }


       if(jr > 0){
            for (std::size_t j=0; j<jr; ++j){
                F(i,jdxr(j)) = Fr(i, j); }
       }
       
 
       }
        
       Fl.clear();
       Fr.clear();
   



       if (ptr[l]->m_typ==1 && tol>0)
       {
           {


            if (nnz(F)>numel(F)/4)
            {
                std::tie(A,B,flag) = aca(F, m_tol/10);
//		std::cout << "ACA Compression " << flag << std::endl;

            }
            else
            {
//		std::cout << "QRSVD Compression" << std::endl;
                std::tie(A,B) = qrsvd(F,eye<T>(size(F, 2)),m_tol/10);
                (size(A,2)<std::min(size(F, 1), size(F, 2))/2)?(flag = true):(flag = false);
            }
            if (flag)
            {
                matrix<T> tmp = F - mtimes(A,B);
                double errLinf = norm(tmp,"inf")/norm(F,"inf");
                double errL2   = norm(tmp,"2")/norm(F,"2");

//		std::cout << "Errors: " << errLinf << "\t" << errL2 << std::endl;
//		std::cout << "Dimensions: " << numel(idx[l]) << " " << numel(jdx[l]) << " " << numel(jl) << std::endl;
                if (errLinf<m_tol && errL2<m_tol)
                {
          	    ptr[l]->m_dat[0] = A;
                    ptr[l]->m_dat[1] = B;
                    ptr[l]->m_typ    = 1;
                }
                else{
                   
                    ptr[l]->m_dat[0] = F;
                    ptr[l]->m_typ    = 2;
                }
            }
            else if (!flag)
            {
               

                    ptr[l]->m_dat[0] = F;
		    ptr[l]->m_typ = 2;
            }



           }
       }
       else if (ptr[l]->m_typ==2 || tol<=0)
       {
           {
           ptr[l]->m_dat[0] = F;
           ptr[l]->m_typ    = 2;

           }
       }
       else
       {
           error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
       }


       
   }
   
   std::cout << "\n" << std::endl;

}

//==========================================================================
// [.hmatrix]
/// Empty contructor using tree, define far leaves as compressed and close
/// leaves as dense, fill with nothing.
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(bintree<S>const& Xb, bintree<S>const& Yb, double tol)
{
    // Global Indices
    m_siz = {size(Xb.crd(),1),size(Yb.crd(),1)};
    m_tol = tol;
    
    // Far leaves
    if (Xb.box().isfar(Yb.box()))
    {
        m_typ = 1;
    }
    // Full leaves
    else if (Xb.isleaf() || Yb.isleaf())
    {
        m_typ = 2;
    }
    // Recursive leaves
    else
    {
        m_chd.resize(4,hmatrix<T>());
        for (int h=0; h<4; ++h)
        {
            m_row[h] = Xb.ind(h/2);
            m_col[h] = Yb.ind(h%2);
            m_chd[h] = hmatrix<T>(Xb.sub(h/2),Yb.sub(h%2),tol);
        }
        m_typ = 0;
    }
}

//==========================================================================
// [.hmatrix]
/// Sparse matrix contructor using tree, fill non-zeros leaves with full
/// matrix and empty leaves with low-rank representation of a zeros block.
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(bintree<S>const& Xb, bintree<S>const& Yb, double tol, smatrix<T>const& Ms)
{
    // Dimensions compatibility
    if (size(Xb.crd(),1)!=size(Ms,1) || size(Yb.crd(),1)!=size(Ms,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    
    // Global data
    m_siz = size(Ms);
    m_tol = tol;
    
    // Empty leaf
    if (nnz(Ms)==0)
    {
        m_typ    = 1;
        m_dat[0] = zeros<T>(m_siz(0),1);
        m_dat[1] = zeros<T>(1,m_siz(1));
    }
    // Full leaf
    else if (Xb.isleaf() || Yb.isleaf())
    {
        m_typ    = 2;
        m_dat[0] = full(Ms);
    }
    // Recursion
    else
    {
        m_typ = 0;
        m_chd.resize(4,hmatrix<T>());
        for (int h=0; h<4; ++h)
        {
            m_row[h] = Xb.ind(h/2);
            m_col[h] = Yb.ind(h%2);
            m_chd[h] = hmatrix<T>(Xb.sub(h/2),Yb.sub(h%2),tol,eval(Ms(m_row[h],m_col[h])));
        }
    }
}

//==========================================================================
// [.hmatrix]
/// Projector constructor : Ml * Mh * Mr
template<typename T>
template<typename S>
hmatrix<T>::hmatrix(bintree<S>const& Xb, bintree<S>const& Yb, smatrix<T>const& Ml, hmatrix<T>const& Mh, smatrix<T>const& Mr)
{
    // Dimensions compatibility
    if (size(Xb.crd(),1)!=size(Ml,1) || size(Yb.crd(),1)!=size(Mr,2) ||
        size(Ml,2)!=size(Mh,1) || size(Mh,2)!=size(Mr,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    
    // Global Indices
    m_siz = {size(Xb.crd(),1),size(Yb.crd(),1)};
    m_tol = Mh.getTolerance();

    // Empty leaf
    if (nnz(Ml)==0 || nnz(Mr)==0)
    {
        m_typ       = 1;
        m_dat[0]    = zeros<T>(m_siz(0),1);
        m_dat[1]    = zeros<T>(1,m_siz(1));
    }
    // Far leaf
    else if (Mh.getType()==1)
    {
        m_typ    = 1;
        m_dat[0] = mtimes(Ml,Mh.getData(0));
        m_dat[1] = mtimes(Mh.getData(1),Mr);
    }
    // Full leaves
    else if (Mh.getType()==2 && (Xb.isleaf() || Yb.isleaf()))
    {
        m_typ    = 2;
        m_dat[0] = mtimes(mtimes(Ml,Mh.getData(0)),Mr);
    }
    // Full leaf from Mh
    //   M = (A1) (B1) (C1 C2)
    //       (A2)
    else if (Mh.getType()==2)
    {
        // Indices of sub matrices
        matrix<std::size_t> idx, jdx, tmp;
        std::tie(tmp,jdx) = ind2sub(size(Ml),index(Ml));
        std::tie(idx,tmp) = ind2sub(size(Mr),index(Mr));
        idx = unique(idx);
        jdx = unique(jdx);

        // Sub matrices
        smatrix<T> A = eval(Ml(row(Ml),jdx));
        smatrix<T> B = sparse(eval(Mh.getData(0)(jdx,idx)));
        smatrix<T> C = eval(Mr(idx,col(Mr)));

        // Product
        smatrix<T> Ms = mtimes(mtimes(A,B),C);
        (*this) = hmatrix<double>(Xb,Yb,m_tol,Ms);
    }
    // Full leaf from tree
    //   M = (A1 A2) (B1 B2) (C1)
    //               (B3 B4) (C2)
    else if (Xb.isleaf() || Yb.isleaf())
    {
        // Indices of sub matrices
        matrix<std::size_t> idx, jdx, tmp;
        std::tie(tmp,jdx) = ind2sub(size(Ml),index(Ml));
        std::tie(idx,tmp) = ind2sub(size(Mr),index(Mr));
        idx = unique(idx);
        jdx = unique(jdx);

        // Sub matrices
        matrix<T> Bf = zeros(length(jdx),length(idx));
        Mh.hfull2(Bf,jdx,idx);
        smatrix<T> As = eval(Ml(row(Ml),jdx));
        smatrix<T> Cs = eval(Mr(idx,col(Mr)));

        // Product
        m_typ = 2;
        m_dat[0] = mtimes(mtimes(As,Bf),Cs);
    }
    // Recursive leaf
    else if (Mh.getType()==0)
    {
        // Prepare blocks
        m_typ = 0;
        m_chd.resize(4,hmatrix<T>());
        for (int h=0; h<4; ++h)
        {
            m_row[h] = Xb.ind(h/2);
            m_col[h] = Yb.ind(h%2);
            m_chd[h] = hmatrix<T>(length(m_row[h]),length(m_col[h]),0);
        }
        
        // Recursion
        matrix<std::size_t> g = {0,2,1,3};
        for (int h=0; h<4; ++h)
        {
           for (int i=0; i<4; ++i)
            {
                m_chd[h] += hmatrix(Xb.sub(h/2), Yb.sub(h%2),
                                    eval(Ml(m_row[h],Mh.getRow(g(i)))),
                                    Mh.getChildren(g(i)),
                                    eval(Mr(Mh.getColumn(g(i)),m_col[h])));
            }
        }
        fusion();
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}


//==========================================================================//
//                               FUNCTIONS                                  //
//==========================================================================//
//==========================================================================
// [hmatrix.full2lowrank]
/// Convert full leaf to low-rank using ACA compression.
template<typename T>
void hmatrix<T>::full2lowrank()
{
    if (m_typ==2 && m_tol>0)
    {
        matrix<T> A, B;
        bool flag;
        // Full
        if (nnz(m_dat[0])>numel(m_dat[0])/4)
        {
            // std::cout << "Number of elements: " << m_dat[0].size(1) << ", " << m_dat[0].size(2) << std::endl; 
            std::tie(A,B,flag) = aca(m_dat[0],m_tol/10);
        }
        // "Sparse"
        else
        {
            std::tie(A,B) = qrsvd(m_dat[0],eye<T>(m_siz(1)),m_tol/10);
            (size(A,2)<min(m_siz)/2)?(flag = true):(flag = false);
        }
        // Check and update
        if (flag)
        {
            matrix<T> tmp = m_dat[0] - mtimes(A,B);
            double errLinf = norm(tmp,"inf")/norm(m_dat[0],"inf");
            double errL2   = norm(tmp,"2")/norm(m_dat[0],"2");
            if (errLinf<m_tol && errL2<m_tol)
            {
                m_dat[0] = A;
                m_dat[1] = B;
                m_typ    = 1;
            }
        }
        else if (!flag)
        {
            // std::cout << "ACA compression failed: " << std::endl;
        }
    }
}

//==========================================================================
// [hmatrix.fusion]
/// Low-rank fusion is possible if and only if hierarchical leaf has four
/// low-rank blocks.
template<typename T>
void hmatrix<T>::fusion()
{
    if (m_typ==0)
    {
        // Leaves
        matrix<int> leaf(1,4);
        for (int h=0; h<4; ++h)
        {
            leaf(h) = m_chd[h].m_typ;
        }

        // Low-rank fusion
        if (nnz(leaf==1)==4)
        {
            std::size_t l=0, n=0;
            matrix<std::size_t> r(1,4);
            for (int h=0; h<4; ++h)
            {
                r(h) = size(m_chd[h].m_dat[0],2);
                n += numel(m_chd[h].m_dat[0]) + numel(m_chd[h].m_dat[1]);
            }
            matrix<T> A(m_siz(0),sum(r)), B(sum(r),m_siz(1));
            for (int h=0; h<4; ++h)
            {
                A(m_row[h],range(l,l+r(h))) = m_chd[h].m_dat[0];
                B(range(l,l+r(h)),m_col[h]) = m_chd[h].m_dat[1];
                l += r(h);
            }
            std::tie(A,B) = qrsvd(A,B,m_tol);
            if (numel(A)+numel(B)<n)
            {
                m_typ    = 1;
                m_dat[0] = A;
                m_dat[1] = B;
                for (int h=0; h<4; ++h)
                {
                    m_row[h].resize(0,0);
                    m_col[h].resize(0,0);
                }
                m_chd.clear();
            }
        }
    }
}

//==========================================================================
// [hmatrix.hfull]
/// In-place conversion to dense matrix :
/// Ah -> M.
template<typename T>
void hmatrix<T>::hfull(matrix<T>& M, matrix<std::size_t> I, matrix<std::size_t> J) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].hfull(M,eval(I(m_row[h])),eval(J(m_col[h])));
        }
    }
    else if (m_typ==1) {M(I,J) = mtimes(m_dat[0],m_dat[1]);}
    else if (m_typ==2) {M(I,J) = m_dat[0];}
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.full2]
/// In-place conversion to dense sub-matrix :
/// Ah(I,J) -> M.
template<typename T>
void hmatrix<T>::hfull2(matrix<T>& M, matrix<std::size_t> idx, matrix<std::size_t> jdx,
                      matrix<std::size_t> I, matrix<std::size_t> J) const
{
    if (isempty(I)) {I = range(0,size(M,1));}
    if (isempty(J)) {J = range(0,size(M,2));}
    if (numel(I)!=numel(idx) || numel(J)!=numel(jdx))
    {
	std::cout << numel(I) << " " << numel(idx) << " " << numel(J) << " " << numel(jdx) << std::endl;
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            matrix<std::size_t> idxh, jdxh, Ih, Jh, Itmp, Jtmp;
            std::tie(idxh,Itmp) = argintersect(m_row[h],idx);
            std::tie(jdxh,Jtmp) = argintersect(m_col[h],jdx);
            if (!isempty(idxh) && !isempty(jdxh))
            {
                Ih = eval(I(Itmp));
                Jh = eval(J(Jtmp));
                m_chd[h].hfull2(M,idxh,jdxh,Ih,Jh);
            }
        }
    }
    else if (m_typ==1)
    {
        matrix<T> A = eval(m_dat[0](idx,col(m_dat[0])));
        matrix<T> B = eval(m_dat[1](row(m_dat[1]),jdx));
        M(I,J) = mtimes(A,B);
    }
    else if (m_typ==2)
    {
        M(I,J) = eval(m_dat[0](idx,jdx));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.hinv]
/// In-place hierarchical matrix inversion :
/// Ah -> Ah^(-1).
template<typename T>
void hmatrix<T>::hinv()
{
    if (m_siz(0)!=m_siz(1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    if (m_typ==0)
    {
        // Am1 -> M11
        // std::cout << "we are here" << std::endl;

        m_chd[0].hinv();
        
        // Usefull data
        hmatrix<T> Am1B = mtimes(m_chd[0],m_chd[1]);
        hmatrix<T> CAm1 = mtimes(m_chd[2],m_chd[0]);
        
        // std::cout << "Definitions with matrix product"<< std::endl; 
        // S = D - C * Am1 * B -> M22
        m_chd[3].tgehmhm((T)-1,m_chd[2],Am1B,(T)1);
        
        // Sm1 -> M22
        m_chd[3].hinv();
        
        // - Am1 * B * Sm1 -> M12
        m_chd[1].tgehmhm((T)-1,Am1B,m_chd[3],(T)0);
        
        // Am1 + Am1 * B * Sm1 * C * Am1 -> M11
        m_chd[0].tgehmhm((T)-1,m_chd[1],CAm1,(T)1);
        
        // - Sm1 * C * Am1 -> M21
        m_chd[2].tgehmhm((T)-1,m_chd[3],CAm1,(T)0);
        

        // std::cout << "This is fine"<< std::endl; 
        // Leaf fusion
        fusion();
    }
    else if (m_typ==2)
    {
        // std::cout << "Inverse of full matrix" << std::endl;
        m_dat[0] = inv(m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.hllowsolve]
/// In-place hmatrix-hmatrix lower-triangular left resolution :
/// Ah -> Lh^(-1) * Ah.
template<typename T>
void hmatrix<T>::hllowsolve(hmatrix<T>const& Lh)
{
    if (m_siz(0)!=Lh.getSize(0))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    // hmatrix \ hmatrix -> hmatrix
    if (Lh.getType()==0 && m_typ==0)
    {
        // X11 -> L11 \ B11
        m_chd[0].hllowsolve(Lh.getChildren(0));
        
        // X12 -> L11 \ B12
        m_chd[1].hllowsolve(Lh.getChildren(0));
        
        // X21 -> L22 \ (B21 - L21*X11)
        m_chd[2].tgehmhm((T)-1,Lh.getChildren(2),m_chd[0],(T)1);
        m_chd[2].hllowsolve(Lh.getChildren(3));
        
        // X22 -> L22 \ (B22 - L21 * X12)
        m_chd[3].tgehmhm((T)-1,Lh.getChildren(2),m_chd[1],(T)1);
        m_chd[3].hllowsolve(Lh.getChildren(3));
        
        // Fusion
        fusion();
    }
    // hmatrix \ compr or full -> compr or full
    else if (Lh.getType()==0 && m_typ>0)
    {
        Lh.hllowsolve(m_dat[0]);
    }
    // full \ compr or full -> compr or full
    else if (Lh.getType()==2 && m_typ>0)
    {
        m_dat[0] = linsolve(Lh.getData(0),m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.hllowsolve]
/// In-place hmatrix-matrix lower-triangular left resolution :
/// B -> Lh^(-1) * B.
template<typename T>
void hmatrix<T>::hllowsolve(matrix<T>& B, matrix<std::size_t> I) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(I))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    matrix<std::size_t> J = col(B);
    if (m_typ==0)
    {
        // X1 -> L11 \ B1
        matrix<std::size_t> I1 = eval(I(m_row[0]));
        m_chd[0].hllowsolve(B,I1);
        
        // X2 -> L22 \ (B2 - L21*X1)
        matrix<std::size_t> I2 = eval(I(m_row[2]));
        B(I2,J) = eval(B(I2,J)) - mtimes(m_chd[2],eval(B(I1,J)));
        m_chd[3].hllowsolve(B,I2);
    }
    else if (m_typ==2)
    {
        B(I,J) = linsolve(m_dat[0],eval(B(I,J)));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.hlupsolve]
/// In-place hmatrix-matrix upper-triangular left resolution :
/// B -> Uh^(-1) * B.
template<typename T>
void hmatrix<T>::hlupsolve(matrix<T>& B, matrix<std::size_t> I) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(I))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    matrix<std::size_t> J = col(B);
    if (m_typ==0)
    {
        // X2 -> U22 \ B2
        matrix<std::size_t> I2 = eval(I(m_row[2]));
        m_chd[3].hlupsolve(B,I2);
        
        // X1 -> U11 \ (B1 - U12*X2)
        matrix<std::size_t> I1 = eval(I(m_row[0]));
        B(I1,J) = eval(B(I1,J)) - mtimes(m_chd[1],eval(B(I2,J)));
        m_chd[0].hlupsolve(B,I1);
    }
    else if (m_typ==2)
    {
        B(I,J) = linsolve(m_dat[0],eval(B(I,J)));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.hlmtimes]
/// In-place hmatrix-matrix product :
/// C -> C + Ah * B.
template<typename T>
void hmatrix<T>::hlmtimes(matrix<T>const& B, matrix<T>& C, matrix<std::size_t> I, matrix<std::size_t> J) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].hlmtimes(B,C,eval(I(m_row[h])),eval(J(m_col[h])));
        }
    }
    else if (m_typ==1)
    {
        matrix<std::size_t> c = col(B);
        C(I,c) = eval(C(I,c)) + mtimes(m_dat[0],mtimes(m_dat[1],eval(B(J,c))));
    }
    else if (m_typ==2)
    {
        matrix<std::size_t> c = col(B);
        C(I,c) = eval(C(I,c)) + mtimes(m_dat[0],eval(B(J,c)));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.hlu]
/// In-place hmatrix-hmatrix lower-upper factorization :
/// Ah -> mtimes(Ah,Uh).
template<typename T>
void hmatrix<T>::hlu(hmatrix<T>& Uh)
{
    if (m_siz(0)!=m_siz(1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    Uh.setTolerance() = m_tol;
    Uh.setSize()      = m_siz;
    if (m_typ==0)
    {
        // Initialize Uh
        Uh.setType() = 0;
        Uh.setChildren().resize(4,hmatrix<T>());
        for (int h=0; h<4; ++h)
        {
            Uh.setRow(h)    = m_row[h];
            Uh.setColumn(h) = m_col[h];
        }
        
        // [L11,U11] -> M11
        
        // std::cout << "Maybe here" << std::endl;
        m_chd[0].hlu(Uh.setChildren(0));
        
        // U12 -> L11 \ M12   &&   L12 -> 0
        
        // std::cout << "Or here" << std::endl;
        Uh.setChildren(1) = std::move(m_chd[1]);
        
        
        // std::cout << "Or here maybe" << std::endl;
        Uh.setChildren(1).hllowsolve(m_chd[0]);
        m_chd[1] = hmatrix<T>(numel(m_row[1]),numel(m_col[1]),m_tol);
        
        // L21 -> M21 / U11   &&   U21 -> 0
        
        // std::cout << "hrupsolve" << std::endl;
        m_chd[2].hrupsolve(Uh.getChildren(0));
        Uh.setChildren(2) = hmatrix<T>(numel(m_row[2]),numel(m_col[2]),m_tol);
        
        // M22 -> M22 - L21*U12
        
        // std::cout << "tgehmhm" << std::endl;
        m_chd[3].tgehmhm((T)-1,m_chd[2],Uh.getChildren(1),(T)1);
        
        // [L22,U22] -> M22
        
        // std::cout << "Is it in the H-matrix LU?" << std::endl;
        m_chd[3].hlu(Uh.setChildren(3));
        
        // Fusion
        fusion();
        Uh.fusion();
    }
    else if (m_typ==2)
    {
        Uh.setType() = 2;
        // std::cout << "Is it in the full matrix LU?" << std::endl;
        std::tie(m_dat[0],Uh.setData(0)) = lu(m_dat[0]);
        
        // std::cout << "Not here" << std::endl;
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.hrmtimes]
/// In-place matrix-hmatrix product :
/// C -> C + A * Bh.
template<typename T>
void hmatrix<T>::hrmtimes(matrix<T>const& A, matrix<T>& C, matrix<std::size_t> I, matrix<std::size_t> J) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].hrmtimes(A,C,eval(I(m_row[h])),eval(J(m_col[h])));
        }
    }
    else if (m_typ==1)
    {
        matrix<std::size_t> r = row(A);
        C(r,J) = eval(C(r,J)) + mtimes(mtimes(eval(A(r,I)),m_dat[0]),m_dat[1]);
    }
    else if (m_typ==2)
    {
        matrix<std::size_t> r = row(A);
        C(r,J) = eval(C(r,J)) + mtimes(eval(A(r,I)),m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.hrupsolve]
/// In-place hmatrix-hmatrix upper-triangular right resolution :
/// Bh -> Bh * Uh^(-1).
template<typename T>
void hmatrix<T>::hrupsolve(hmatrix<T>const& Uh)
{
    if (m_siz(1)!=Uh.getSize(1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    // hmatrix \ hmatrix -> hmatrix
    if (Uh.getType()==0 && m_typ==0)
    {
        // X11 -> B11 / U11
        m_chd[0].hrupsolve(Uh.getChildren(0));
        
        // X21 -> B21 / U11
        m_chd[2].hrupsolve(Uh.getChildren(0));
        
        // X12 -> (B12 - X11*U12) / U22
        m_chd[1].tgehmhm((T)-1,m_chd[0],Uh.getChildren(1),(T)1);
        m_chd[1].hrupsolve(Uh.getChildren(3));
        
        // X22 -> (B22 - X21*U12) / U22
        m_chd[3].tgehmhm((T)-1,m_chd[2],Uh.getChildren(1),(T)1);
        m_chd[3].hrupsolve(Uh.getChildren(3));
        
        // Fusion
        fusion();
    }
    // compr / hmatrix -> compr
    else if (Uh.getType()==0 && m_typ==1)
    {
        Uh.hrupsolve(m_dat[1]);
    }
    // full / hmatrix -> full
    else if (Uh.getType()==0 && m_typ==2)
    {
        Uh.hrupsolve(m_dat[0]);
    }
    // compr / full -> compr
    else if (Uh.getType()==2 && m_typ==1)
    {
        m_dat[1] = linsolve(transpose(Uh.getData(0)),transpose(m_dat[1]));
        m_dat[1] = transpose(m_dat[1]);
    }
    else if (Uh.getType()==2 && m_typ==2)
    {
        m_dat[0] = linsolve(transpose(Uh.getData(0)),transpose(m_dat[0]));
        m_dat[0] = transpose(m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.hrupsolve]
/// In-place matrix-hmatrix upper-triangular right resolution :
/// B -> B * Uh^(-1).
template<typename T>
void hmatrix<T>::hrupsolve(matrix<T>& B, matrix<std::size_t> J) const
{
    if (isempty(J)) {J = range(0,m_siz(0));}
    if (m_siz(0)!=numel(J) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"H-matrix is not square.");
    }
    matrix<std::size_t> I = row(B);
    if (m_typ==0)
    {
        // X1 -> B1 / U11
        matrix<std::size_t> J1 = eval(J(m_col[0]));
        m_chd[0].hrupsolve(B,J1);
        
        // X2 -> (B2 - X1*U12) / U22
        matrix<std::size_t> J2 = eval(J(m_col[1]));
        B(I,J2) = eval(B(I,J2)) - mtimes(eval(B(I,J1)),m_chd[1]);
        m_chd[3].hrupsolve(B,J2);
    }
    else if (m_typ==2)
    {
        B(I,J) = transpose(linsolve(transpose(m_dat[0]),transpose(eval(B(I,J)))));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.htranspose]
/// In-place hmatrix transposition :
/// Bh -> Bh^t.
template<typename T>
void hmatrix<T>::htranspose()
{
    m_siz = {m_siz(1),m_siz(0)};
    if (m_typ==0)
    {
        auto Mh  = std::move(m_chd[1]);
        m_chd[1] = std::move(m_chd[2]);
        m_chd[2] = std::move(Mh);
        for (int h=0; h<4; ++h)
        {
            m_chd[h].htranspose();
        }
        matrix<std::size_t> row[4] = {m_row[0],m_row[1],m_row[2],m_row[3]};
        m_row[0] = m_col[0];  m_row[1] = m_col[2];
        m_row[2] = m_col[1];  m_row[3] = m_col[3];
        m_col[0] = row[0];  m_col[1] = row[2];
        m_col[2] = row[1];  m_col[3] = row[3];
    }
    else  if (m_typ==1)
    {
        matrix<T> At = transpose(m_dat[0]);
        m_dat[0] = transpose(m_dat[1]);
        m_dat[1] = At;
    }
    else  if (m_typ==2)
    {
        m_dat[0] = transpose(m_dat[0]);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}


//==========================================================================
// [hmatrix.leafptr]
/// Extract leaves pointers to compute externaly leaves data.
template<typename T>
void hmatrix<T>::leafptr(std::vector<hmatrix<T>*>& ptr,
                         std::vector<matrix<std::size_t>>& idx, std::vector<matrix<std::size_t>>& jdx,
                         matrix<std::size_t> I, matrix<std::size_t> J)
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].leafptr(ptr,idx,jdx,eval(I(m_row[h])),eval(J(m_col[h])));
        }
    }
    else
    {
        //std::shared_ptr<hmatrix<T>> sharedThis(this);
        ptr.push_back(this);
        idx.push_back(I);
        jdx.push_back(J);
    }
}

//==========================================================================
// [hmatrix.recompress]
/// Recompression of dense and low-rank leaves with fusion.
template<typename T>
void hmatrix<T>::recompress(double tol)
{
    m_tol = std::max(tol,m_tol);
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].recompress(tol);
        }
        fusion();
    }
    else if (m_typ==1)
    {
        std::tie(m_dat[0],m_dat[1]) = qrsvd(m_dat[0],m_dat[1],m_tol);
    }
    else if (m_typ==2)
    {
        full2lowrank();
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}



//==========================================================================
// [hmatrix.clear]
/// Clear matrix data.
template<typename T>
void hmatrix<T>::clear()
{
    
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].clear();
        }
    }
    else if (m_typ==1)
    {
        m_dat[0].clear();
        m_dat[1].clear();
    }
    else if (m_typ==2)
    {
        m_dat[0].clear();
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.spydata]
/// Extract sparse matrix Sf representing full leaves and Sc compressed ones.
template<typename T>
void hmatrix<T>::spydata(smatrix<logical>& Sf, smatrix<logical>& Sc, matrix<std::size_t> I, matrix<std::size_t> J) const
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        matrix<std::size_t> I1=range(0,numel(I)/2), I2=range(numel(I)/2,numel(I));
        matrix<std::size_t> J1=range(0,numel(J)/2), J2=range(numel(J)/2,numel(J));
        m_chd[0].spydata(Sf,Sc,eval(I(I1)),eval(J(J1)));
        m_chd[1].spydata(Sf,Sc,eval(I(I1)),eval(J(J2)));
        m_chd[2].spydata(Sf,Sc,eval(I(I2)),eval(J(J1)));
        m_chd[3].spydata(Sf,Sc,eval(I(I2)),eval(J(J2)));
    }
    else
    {
        smatrix<logical> S(m_siz(0),m_siz(1));
        S(row(S),{0,m_siz(1)-1}) = spones<logical>(m_siz(0),2);
        S({0,m_siz(0)-1},col(S)) = spones<logical>(2,m_siz(1));
        if (m_typ==1)
        {
            Sc(I,J) = S;
        }
        else if (m_typ==2)
        {
            Sf(I,J) = S;
        }
        else
        {
            error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
        }
    }
}

//==========================================================================
// [hmatrix.stat]
/// Extract hierarchical storage statistics, with :
/// Leaf : #hierarchical #compressed #full
/// Rate : #values/numel * 100
/// tol  : accuracy
template<typename T>
void hmatrix<T>::stat(matrix<std::size_t>& info) const
{
    if (m_typ==-1) {return;}
    info(m_typ) += 1;
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].stat(info);
        }
    }
    else if (m_typ==1) {info(3) += numel(m_dat[0])+numel(m_dat[1]);}
    else if (m_typ==2) {info(3) += numel(m_dat[0]);}
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.tgeabhm]
/// In-place hmatrix product with low-rank matrix representation :
/// Ch -> alpha*(A*B) + beta*Ch.
template<typename T>
void hmatrix<T>::tgeabhm(T alpha, matrix<T>const& A, matrix<T>const& B, T beta, matrix<std::size_t> I, matrix<std::size_t> J)
{
    if (isempty(I)) {I = range(0,m_siz(0));}
    if (isempty(J)) {J = range(0,m_siz(1));}
    if (m_siz(0)!=numel(I) || m_siz(1)!=numel(J))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h].tgeabhm(alpha,A,B,beta,eval(I(m_row[h])),eval(J(m_col[h])));
        }
        fusion();
    }
    else if (m_typ==1)
    {
        matrix<std::size_t> c = col(A);
        m_dat[0] *= beta;
        m_dat[0] = horzcat(m_dat[0],alpha*eval(A(I,c)));
        m_dat[1] = vertcat(m_dat[1],eval(B(c,J)));
        std::tie(m_dat[0],m_dat[1]) = qrsvd(m_dat[0],m_dat[1],m_tol);
    }
    else if (m_typ==2)
    {
        matrix<std::size_t> c = col(A);
        m_dat[0] *= beta;
        m_dat[0] += mtimes(alpha*eval(A(I,c)),eval(B(c,J)));
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}

//==========================================================================
// [hmatrix.tgehmhm]
/// In-place hmatrix product :
/// Ch -> alpha*(Ah*Bh) + beta*Ch.
template<typename T>
void hmatrix<T>::tgehmhm(T alpha, hmatrix<T>const& Ah, hmatrix<T>const& Bh, T beta)
{
    if (m_siz(0)!=Ah.getSize(0) || m_siz(1)!=Bh.getSize(1) || Ah.getSize(1)!=Bh.getSize(0))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    // -------------------------------------------------------------- Ch is hmatrix
    // alpha * hmatrix * hmatrix + beta * hmatrix -> hmatrix
    if (Ah.getType()==0 && Bh.getType()==0 && m_typ==0)
    {
        //        std::cout << " H1 ";
        // alpha * (A11*B11 + A12*B21) + beta C11 -> C11
        m_chd[0].tgehmhm(alpha,Ah.getChildren(0),Bh.getChildren(0),beta);
        m_chd[0].tgehmhm(alpha,Ah.getChildren(1),Bh.getChildren(2),(T)1);
        
        // alpha * (A11*B12 + A12*B22) + beta C12 -> C12
        m_chd[1].tgehmhm(alpha,Ah.getChildren(0),Bh.getChildren(1),beta);
        m_chd[1].tgehmhm(alpha,Ah.getChildren(1),Bh.getChildren(3),(T)1);
        
        // alpha * (A21*B11 + A22*B21) + beta C21 -> C21
        m_chd[2].tgehmhm(alpha,Ah.getChildren(2),Bh.getChildren(0),beta);
        m_chd[2].tgehmhm(alpha,Ah.getChildren(3),Bh.getChildren(2),(T)1);
        
        // alpha * (A21*B12 + A22*B22) + beta C22 -> C223
        m_chd[3].tgehmhm(alpha,Ah.getChildren(2),Bh.getChildren(1),beta);
        m_chd[3].tgehmhm(alpha,Ah.getChildren(3),Bh.getChildren(3),(T)1);
        
        // Fusion
        fusion();
    }
   
   
   
    // alpha * compr * hmatrix + beta * hmatrix -> hmatrix
    else if (Ah.getType()==1 && Bh.getType()==0 && m_typ==0)
    {
        //        std::cout << " H2 ";
        (*this).tgeabhm(alpha,Ah.getData(0),mtimes(Ah.getData(1),Bh),beta);
    }
    // alpha * hmatrix * compr + beta * hmatrix -> hmatrix
    else if (Ah.getType()==0 && Bh.getType()==1 && m_typ==0)
    {
        //        std::cout << " H3 ";
        (*this).tgeabhm(alpha,mtimes(Ah,Bh.getData(0)),Bh.getData(1),beta);
    }
    // alpha * compr * hmatrix + beta * hmatrix -> hmatrix  ------------- NEW!
    else if (Ah.getType()==2 && Bh.getType()==0 && m_typ==0)
    {
	std::cout << "I am new (2 0 0) here!" << std::endl;
        //        std::cout << " H2 NEW FULL";
        (*this).tgeabhm(alpha,Ah.getData(0),full(Bh),beta);
    }
    // alpha * hmatrix * compr + beta * hmatrix -> hmatrix -------------- NEW!
    else if (Ah.getType()==0 && Bh.getType()==2 && m_typ==0)
    {
	std::cout << "I am new (0 2 0) here!" << std::endl;
        //        std::cout << " H3 NEW FULL";
        (*this).tgeabhm(alpha,full(Ah),Bh.getData(0),beta);
    }
    // alpha * hmatrix * compr + beta * hmatrix -> hmatrix -------------- NEW!
    else if (Ah.getType()==2 && Bh.getType()==2 && m_typ==0)
    {
	std::cout << "I am new (2 2 0) here!" << std::endl;
        //        std::cout << " H3 NEW FULL";
        (*this).tgeabhm(alpha,Ah.getData(0),Bh.getData(0),beta);
    }
    // alpha * full * compr + beta * hmatrix -> hmatrix   
    else if (Ah.getType()==2 && Bh.getType()==1 && m_typ==0)
    {
        std::cout << "I am new (2 1 0) here!" << std::endl;
        //        std::cout << " H3 ";
        (*this).tgeabhm(alpha,mtimes(Ah.getData(0),Bh.getData(0)),Bh.getData(1),beta);
    }
    // alpha * compr * full + beta * hmatrix -> hmatrix   
    else if (Ah.getType()==1 && Bh.getType()==2 && m_typ==0)
    {
        std::cout << "I am new (1 2 0) here!" << std::endl;
        //        std::cout << " H3 ";
        (*this).tgeabhm(alpha,Ah.getData(0),mtimes(Ah.getData(1),Bh.getData(0)),beta);
    }
    else if (Ah.getType()==1 && Bh.getType()==1 && m_typ==0)
    {
        //        std::cout << " H4 ";
        (*this).tgeabhm(alpha,mtimes(Ah.getData(0),mtimes(Ah.getData(1),Bh.getData(0))),Bh.getData(1),beta);
    }
    // -------------------------------------------------------------- Ch is compr
    // alpha * hmatrix * hmatrix + beta * compr -> h-matrix
    else if (Ah.getType()==0 && Bh.getType()==0 && m_typ==1)
    {
        //        std::cout << " C1 ";
        m_typ = 0;
        m_chd.resize(4,hmatrix<T>());
        for (int h=0; h<4; ++h)
        {
            m_row[h] = Ah.getRow(h);
            m_col[h] = Bh.getColumn(h);
            m_chd[h] = hmatrix(numel(m_row[h]),numel(m_col[h]),m_tol);
            m_chd[h].setData(0) = eval(m_dat[0](m_row[h],col(m_dat[0])));
            m_chd[h].setData(1) = eval(m_dat[1](row(m_dat[1]),m_col[h]));
        }
        m_dat[0].clear();
        m_dat[1].clear();
        (*this).tgehmhm(alpha,Ah,Bh,beta);
        //        (*this).tgeabhm(alpha,full(Ah),full(Bh),beta);
    }
    // alpha * compr * hmatrix + beta * compr -> compr
    else if (Ah.getType()==1 && Bh.getType()==0 && m_typ==1)
    {
        //        std::cout << " C2 ";
        (*this).tgeabhm(alpha,Ah.getData(0),mtimes(Ah.getData(1),Bh),beta);
    }
    // alpha * hmatrix * compr + beta * compr -> compr
    else if (Ah.getType()==0 && Bh.getType()==1 && m_typ==1)
    {
        //        std::cout << " C3 ";
        (*this).tgeabhm(alpha,mtimes(Ah,Bh.getData(0)),Bh.getData(1),beta);
    }
    // alpha * compr * compr + beta * compr -> compr
    else if (Ah.getType()==1 && Bh.getType()==1 && m_typ==1)
    {
        //        std::cout << " C4 ";
        (*this).tgeabhm(alpha,mtimes(Ah.getData(0),mtimes(Ah.getData(1),Bh.getData(0))),Bh.getData(1),beta);
    }
    // alpha * full * compr + beta * compr -> compr
    else if (Ah.getType()==2 && Bh.getType()==1 && m_typ==1)
    {
        //        std::cout << " C5 ";
        (*this).tgeabhm(alpha,mtimes(Ah.getData(0),Bh.getData(0)),Bh.getData(1),beta);
    }
    // alpha * compr * full + beta * compr -> compr
    else if (Ah.getType()==1 && Bh.getType()==2 && m_typ==1)
    {
        //        std::cout << " C6 ";
        (*this).tgeabhm(alpha,Ah.getData(0),mtimes(Ah.getData(1),Bh.getData(0)),beta);
    }
    // alpha * full * full + beta * compr -> full or compr
    else if (Ah.getType()==2 && Bh.getType()==2 && m_typ==1)
    {
        //        std::cout << " C7 ";
        m_typ = 2;
        m_dat[0] = mtimes(m_dat[0],m_dat[1]);
        tgemm(alpha,Ah.getData(0),Bh.getData(0),beta,m_dat[0]);
        m_dat[1].clear();
        full2lowrank();
    } 
    // alpha * full * full + beta * compr -> full or compr   -------------- NEM! 0 2 1
    else if (Ah.getType()==0 && Bh.getType()==2 && m_typ==1)
    {
	std::cout << "I am new (0 2 1) here!" << std::endl;
        //        std::cout << " C7 ";
        m_typ = 2;
        m_dat[0] = mtimes(m_dat[0],m_dat[1]);
        tgemm(alpha,full(Ah),Bh.getData(0),beta,m_dat[0]);
        m_dat[1].clear();
        full2lowrank();
    }
    // alpha * full * full + beta * compr -> full or compr   -------------- NEM! 2 0 1
    else if (Ah.getType()==2 && Bh.getType()==0 && m_typ==1)
    {
         std::cout << "I am new (2 0 1) here!" << std::endl;
         //        std::cout << " C7 ";
         m_typ = 2;
         m_dat[0] = mtimes(m_dat[0],m_dat[1]);
         tgemm(alpha,Ah.getData(0),full(Bh),beta,m_dat[0]);
         m_dat[1].clear();
         full2lowrank();
     }
    // -------------------------------------------------------------- Ch is full
    // alpha * compr * compr + beta * full -> full
    else if (Ah.getType()==1 && Bh.getType()==1 && m_typ==2)
    {
        //        std::cout << " F1 ";
        matrix<T> tmp = mtimes(Ah.getData(0),mtimes(Ah.getData(1),Bh.getData(0)));
        tgemm(alpha,tmp,Bh.getData(1),beta,m_dat[0]);
    }
    // alpha * full * compr + beta * full -> full
    else if (Ah.getType()==2 && Bh.getType()==1 && m_typ==2)
    {
        //        std::cout << " F2 ";
        matrix<T> tmp = mtimes(Ah.getData(0),Bh.getData(0));
        tgemm(alpha,tmp,Bh.getData(1),beta,m_dat[0]);
    }
    // alpha * compr * full + beta * full -> full
    else if (Ah.getType()==1 && Bh.getType()==2 && m_typ==2)
    {
        //        std::cout << " F3 ";
        matrix<T> tmp = mtimes(Ah.getData(1),Bh.getData(0));
        tgemm(alpha,Ah.getData(0),tmp,beta,m_dat[0]);
    }
    // alpha * full * full + beta * full -> full
    else if (Ah.getType()==2 && Bh.getType()==2 && m_typ==2)
    {
        //        std::cout << " F4 ";
        tgemm(alpha,Ah.getData(0),Bh.getData(0),beta,m_dat[0]);
    }
    else if(m_typ == 2)
    {
	std::cout << "I am new (X X 2) here!" << std::endl;
	tgemm(alpha,full(Ah), full(Bh), beta, m_dat[0]);

    }

    // Unavailable cases
    else
    {
        std::cout << Ah.getType() << " " << Bh.getType() << " " << m_typ << std::endl;
        std::cout << Ah.getSize(0) << " " << Ah.getSize(1) << " " << Bh.getSize(0) << " " << Bh.getSize(1) << " "  << m_siz(0) << " " << m_siz(1) << std::endl;
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
}


//==========================================================================//
//                               OPERATORS                                  //
//==========================================================================//
//==========================================================================
// [operator())]
/// Access to sub matrix with dense conversion.
template<typename T>
matrix<T> hmatrix<T>::operator()(matrix<std::size_t>const& I, matrix<std::size_t>const& J) const
{
    matrix<T> A(length(I),length(J));
    hfull2(A,I,J);
    return A;
}

//==========================================================================
// [operator+=]
/// In-place addition of hmatrix or scalar :
/// Ah -> Ah + b,
/// Ah -> Ah + Bh.
template<typename T>
hmatrix<T>& hmatrix<T>::operator+=(T b)
{
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h] += b;
        }
    }
    else if (m_typ==1)
    {
        m_dat[0] = horzcat(m_dat[0],b*ones<T>(size(m_dat[0],1),1));
        m_dat[1] = vertcat(m_dat[1],ones<T>(1,size(m_dat[1],2)));
    }
    else  if (m_typ==2)
    {
        m_dat[0] += b;
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
    return (*this);
}
template<typename T>
hmatrix<T>& hmatrix<T>::operator+=(hmatrix<T>const& Bh)
{
    if (m_siz(0)!=Bh.getSize(0) || m_siz(1)!=Bh.getSize(1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Dimensions must agree.");
    }
    // hmatrix + hmatrix -> hmatrix
    if (m_typ==0 && Bh.getType()==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h] += Bh.getChildren(h);
        }
        fusion();
    }
    // hmatrix + compr -> hmatrix
    else if (m_typ==0 && Bh.getType()==1)
    {
        tgeabhm((T)1,Bh.getData(0),Bh.getData(1),(T)1);
    }
    // compr + --- -> ---
    else if (m_typ==1)
    {
        matrix<T> A = m_dat[0];
        matrix<T> B = m_dat[1];
        (*this) = Bh;
        tgeabhm((T)1,A,B,(T)1);
    }
    // full + compr -> full or compr
    else if (m_typ==2 && Bh.getType()==1)
    {
        tgeabhm((T)1,Bh.getData(0),Bh.getData(1),(T)1);
    }
    // full + full -> full
    else if (m_typ==2 && Bh.getType()==2)
    {
        m_dat[0] += Bh.getData(0);
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
    return (*this);
}

//==========================================================================
// [operator*=]
/// In-place hmatrix scalar multiplication :
/// Ah -> Ah*b.
template<typename T>
hmatrix<T>& hmatrix<T>::operator*=(T b)
{
    if (m_typ==0)
    {
        for (int h=0; h<4; ++h)
        {
            m_chd[h] *= b;
        }
    }
    else  if (m_typ==1)
    {
        m_dat[0] *= b;
    }
    else  if (m_typ==2)
    {
        m_dat[0] *= b;
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Unavailable case.");
    }
    return (*this);
}

//==========================================================================
// [operator+]
/// External hmatrix addition :
/// Ch = Ah + b,
/// Ch = Ah + Bh.
template<typename T>
inline hmatrix<T> operator+(hmatrix<T>const& Ah, T b)
{
    hmatrix<T> Ch = Ah;
    Ch += b;
    return Ch;
}
template<typename T>
inline hmatrix<T> operator+(T a, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Bh;
    Ch += a;
    return Ch;
}
template<typename T>
inline hmatrix<T> operator+(hmatrix<T>const& Ah, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Ah;
    Ch += Bh;
    return Ch;
}

//==========================================================================
// [operator-]
/// External hmatrix unitary minus or substraction :
/// Ch = -Ah,
/// Ch = Ah - b,
/// Ch = Ah - Bh.
template<typename T>
inline hmatrix<T> operator-(hmatrix<T>const& Ah)
{
    hmatrix<T> Bh = Ah;
    Bh *= ((T)-1);
    return Bh;
}
template<typename T>
inline hmatrix<T> operator-(hmatrix<T>const& Ah, T b)
{
    hmatrix<T> Ch = Ah;
    Ch += (-b);
    return Ch;
}
template<typename T>
inline hmatrix<T> operator-(T a, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Bh;
    Ch *= ((T)-1);
    Ch += a;
    return Ch;
}
template<typename T>
inline hmatrix<T> operator-(hmatrix<T>const& Ah, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Bh;
    Ch *= -1;
    Ch += Ah;
    return Ch;
}

//==========================================================================
// [operator*]
/// External hmatrix terms by terms scalar multiplication :
/// Ch = Ah * b.
template<typename T>
inline hmatrix<T> operator*(hmatrix<T>const& Ah, T b)
{
    hmatrix<T> Ch = Ah;
    Ch *= b;
    return Ch;
}
template<typename T>
inline hmatrix<T> operator*(T a, hmatrix<T>const& Bh)
{
    hmatrix<T> Ch = Bh;
    Ch *= a;
    return Ch;
}

//==========================================================================
// [operator/]
/// External hmatrix terms by terms scalar division :
/// Ch = Ah / b.
template<typename T>
inline hmatrix<T> operator/(hmatrix<T>const& Ah, T b)
{
    hmatrix<T> Ch = Ah;
    Ch *= ((T)1)/b;
    return Ch;
}

//==========================================================================
// [operator<<]
/// Flux operator to display hmatrix size and statistics using standard
/// C++ convention of disp(Ah) operator.
template<typename T>
std::ostream& operator<<(std::ostream& flux, hmatrix<T> const& Ah)
{
    matrix<std::size_t> info(1,4);
    Ah.stat(info);
    double memBytes = info(3)*sizeof(T);
    flux << "H-matrix " << size(Ah,1) << "x" << size(Ah,2);
    flux << " of type '" << typeid(T).name() << "' (";
    if (memBytes<1e3) {flux << memBytes << " B)";}
    else if (memBytes<1e6) {flux << memBytes/1e3 << " KB)";}
    else {flux << memBytes/1e6 << " MB)";}
    flux << std::endl;
    flux << "  Leaf : " << info(0) << "h " << info(1) << "c " << info(2) << "f " << std::endl;
    flux << "  Rate : " << 100.*info(3)/(size(Ah,1)*size(Ah,2)) << " %" << std::endl;
    flux << "  tol  : " << Ah.getTolerance();
    return flux;
}


//==========================================================================//
//                       MATLAB-LIKE FUNCTION                               //
//==========================================================================//
//==========================================================================
// [disp]
/// Display H-Matrix informations
template<typename T>
void disp(hmatrix<T>const& Ah, int info=2, std::ostream& flux=std::cout)
{
    if (info<=0)
    {
        flux << Ah;
    }
    else if (info==1)
    {
        flux << Ah << std::endl;
    }
    else if (info>=2)
    {
        flux << Ah << std::endl;
    }
}

//==========================================================================
// [full]
/// Hierarchical to dense matrix conversion.
template<typename T>
inline matrix<T> full(hmatrix<T>const& Ah)
{
    matrix<T> A(size(Ah,1),size(Ah,2));
    Ah.hfull(A);
    return A;
}

//==========================================================================
// [gmres]
/// Iterative solver using hierarchical matrix-vector product, eventually
/// using hierarchical preconditionners.
template<typename T>
matrix<T> gmres(hmatrix<T>const& Ah, matrix<T>const& B,
                double tol = 1e-6, std::size_t maxit = 10,
                std::function<matrix<T>(matrix<T>const&)>const& Am1 = std::function<matrix<T>(matrix<T>const&)>(),
                matrix<T>const& X0 = matrix<T>())
{
    if (size(Ah,1)!=size(B,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    std::function<matrix<T>(matrix<T>const&)> Afct;
    Afct = [&Ah](matrix<T>const& X) {return mtimes(Ah,X);};
    return gmres(Afct,B,tol,maxit,Am1,X0);
}
template<typename T>
matrix<T> gmres(hmatrix<T>const& Ah, matrix<T>const& B, double tol, std::size_t maxit,
                hmatrix<T>const& Ahm1, matrix<T>const& X0 = matrix<T>())
{
    if (size(Ah,1)!=size(B,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    std::function<matrix<T>(matrix<T>const&)> Afct, Am1fct;
    Afct   = [&Ah](matrix<T>const& X) {return mtimes(Ah,X);};
    Am1fct = [&Ahm1](matrix<T>const& X) {return mtimes(Ahm1,X);};
    return gmres(Afct,B,tol,maxit,Am1fct,X0);
}
template<typename T>
matrix<T> gmres(hmatrix<T>const& Ah, matrix<T>const& B, double tol, std::size_t maxit,
                hmatrix<T>const& Lh, hmatrix<T>const& Uh, matrix<T>const& X0 = matrix<T>())
{
    if (size(Ah,1)!=size(B,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    std::function<matrix<T>(matrix<T>const&)> Afct, Am1fct;
    Afct   = [&Ah](matrix<T>const& X) {return mtimes(Ah,X);};
    Am1fct = [&Lh,&Uh](matrix<T>const& B)
    {
        matrix<T> X = B;
        Lh.hllowsolve(X);
        Uh.hlupsolve(X);
        return X;
    };
    return gmres(Afct,B,tol,maxit,Am1fct,X0);
}

//==========================================================================
// [hones]
/// Ones hierarchical matrix using low-rank representation.
template<typename T=double>
hmatrix<T> hones(std::size_t m, long n=-1)
{
    if (n==-1) {n=m;}
    return hmatrix<T>(m,n,0,(T)1);
}
template<typename T=double>
hmatrix<T> hones(matrix<std::size_t>const& S)
{
    return hones<T>(S(0),S(1));
}

//==========================================================================
// [hzeros]
/// Zeros hierarchical matrix using low-rank representation.
template<typename T=double>
hmatrix<T> hzeros(std::size_t m, long n=-1)
{
    if (n==-1) {n=m;}
    return hmatrix<T>(m,n,0,(T)0);
}
template<typename T=double>
hmatrix<T> hzeros(matrix<std::size_t>const& S)
{
    return hzeros<T>(S(0),S(1));
}

//==========================================================================
// [inv]
/// Inverse of hierarchical matrix.
template<typename T>
inline hmatrix<T> inv(hmatrix<T>const& Ah)
{
    if (size(Ah,1)!=size(Ah,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    hmatrix<T> Ahm1 = Ah;
    Ahm1.hinv();
    return Ahm1;
}

//==========================================================================
// [linsolve]
/// Solve linear system Ah*X=B where Ah is a hierarchical matrix, using
/// hierarchical LU factorization.
template<typename T>
matrix<T> linsolve(hmatrix<T>const& Ah, matrix<T>const& B)
{
    if (size(Ah,1)!=size(B,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    hmatrix<T> Lh, Uh;
    // std::cout << "First we compute the LU factorization " << std::endl;
    std::tie(Lh,Uh) = lu(Ah);
    
    
    // std::cout << "We copy the RHS into a matrix " << std::endl;
    matrix<T> X = B;
    
    
    // std::cout << "We solve for L " << std::endl;
    Lh.hllowsolve(X);
    
    // std::cout << "We solve for U " << std::endl;
    Uh.hlupsolve(X);
    return X;
}

//==========================================================================
// [lu]
/// Hierarchical lower/upper factorization.
template<typename T>
inline auto lu(hmatrix<T>const& Ah, double tol=0)
{
    if (size(Ah,1)!=size(Ah,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    hmatrix<T> Lh=Ah, Uh;
    if (tol>0) {Lh.recompress(tol);}
    Lh.hlu(Uh);
    return std::make_tuple(Lh,Uh);
}

//==========================================================================
// [mtimes]
/// Hierarchical-matrix multiply with dense or hierarchical matrix.
template<typename T>
inline matrix<T> mtimes(hmatrix<T>const& Ah, matrix<T>const& B)
{
    if (size(Ah,2)!=size(B,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    matrix<T> C(size(Ah,1),size(B,2));
    Ah.hlmtimes(B,C);
    return C;
}
template<typename T>
inline matrix<T> mtimes(matrix<T>const& A, hmatrix<T>const& Bh)
{
    if (size(A,2)!=size(Bh,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    matrix<T> C(size(A,1),size(Bh,2));
    Bh.hrmtimes(A,C);
    return C;
}
template<typename T>
inline hmatrix<T> mtimes(hmatrix<T>const& Ah, hmatrix<T>const& Bh)
{
    if (size(Ah,2)!=size(Bh,1))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    hmatrix<T> Ch(size(Ah,1),size(Bh,2),std::max(Ah.getTolerance(),Bh.getTolerance()));
    tgemm((T)1,Ah,Bh,(T)1,Ch);
    return Ch;
}

//==========================================================================
// [recompress]
/// Leaves recompression according to a new tolerance.
template<typename T>
inline hmatrix<T> recompress(hmatrix<T>const& Ah, double tol)
{
    hmatrix<T> Bh = Ah;
    Bh.recompress(tol);
    return Bh;
}

//==========================================================================
// [size]
/// Size of H-matrix.
template<typename T>
inline matrix<std::size_t> size(hmatrix<T>const& Ah)
{
    return {size(Ah,1),size(Ah,2)};
}
template<typename T>
inline std::size_t size(hmatrix<T>const& Ah, int dim)
{
    return Ah.getSize(dim-1);
}

//==========================================================================
// [spydata]
/// H-matrix leaf structure in sparse matrix format, both for full and
/// compressed leaves.
template<typename T>
inline auto spydata(hmatrix<T>const& Ah)
{
    smatrix<logical> Sf(size(Ah,1),size(Ah,2)), Sc(size(Ah,1),size(Ah,2));
    Ah.spydata(Sf,Sc);
    return std::make_tuple(Sf,Sc);
}

//==========================================================================
// [tgeabm]
/// In-place low-rank matrix product with hierarchical matrix:
/// C = alpha*A*B + beta*Ch,
template<typename T>
inline void tgeabm(T alpha, matrix<T>const& A, matrix<T>const& B, T beta, hmatrix<T>& Ch)
{
    if (size(A,2)!=size(B,1) || size(A,1)!=size(Ch,1) || size(B,2)!=size(Ch,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    Ch.tgeabhm(alpha,A,B,beta);
}

//==========================================================================
// [tgemm]
/// In-place hierarchical matrix product:
/// C = alpha*Ah*Bh + beta*Ch,
template<typename T>
inline void tgemm(T alpha, hmatrix<T>const& Ah, hmatrix<T>const& Bh, T beta, hmatrix<T>& Ch)
{
    if (size(Ah,2)!=size(Bh,1) || size(Ah,1)!=size(Ch,1) || size(Bh,2)!=size(Ch,2))
    {
        error(__FILE__, __LINE__, __FUNCTION__,"Matrix dimensions must agree.");
    }
    Ch.tgehmhm(alpha,Ah,Bh,beta);
}

//==========================================================================
// [transpose]
/// Non-conjugate transpose of hierarchical matrix.
template<typename T>
inline hmatrix<T> transpose(hmatrix<T>const& Ah)
{
    hmatrix<T> Aht = Ah;
    Aht.htranspose();
    return Aht;
}

}

