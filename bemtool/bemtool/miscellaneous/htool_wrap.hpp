#ifndef BEMTOOL_MISC_HTOOLWRAP_HPP
#define BEMTOOL_MISC_HTOOLWRAP_HPP

#include "miscellaneous/htool/htool.hpp"
#include <bemtool/fem/dof.hpp>
#include <bemtool/operator/operator.hpp>
#include <bemtool/operator/block_op.hpp>
#include <bemtool/potential/potential.hpp>

namespace bemtool{

template <typename KernelType, typename Discretization>
class BIO_Generator : public htool::VirtualGenerator<Cplx>{
  Dof<Discretization> dof;
  SubBIOp<BIOp<KernelType>> subV;
  // std::vector<int> boundary;

public:
    BIO_Generator(const Dof<Discretization>& dof0, const double& kappa):VirtualGenerator(NbDof(dof0),NbDof(dof0)), dof(dof0),subV(dof,dof,kappa) {}
    // {boundary=is_boundary_nodes(dof);}

  void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {
    SubBIOp<BIOp<KernelType>> subV_local = subV;
    subV_local.compute_block(M,N,rows,cols,ptr);
  }

};


template <typename KernelType, typename Discretization>
class VIO_Generator : public htool::VirtualGenerator<Cplx>{
  Dof<Discretization> dof;
  SubVIOp<VIOp<KernelType>> V;

public:
    VIO_Generator(const Dof<Discretization>& dof0, const double& kappa):VirtualGenerator(NbDof(dof0),NbDof(dof0)), dof(dof0),V(dof,dof,kappa) {}

  void copy_submatrix(int M, int N, const int *const rows, const int *const cols, Cplx *ptr) const {
    SubVIOp<VIOp<KernelType>> V_local = V;
    V_local.compute_block(M,N,rows,cols,ptr);
  }

};


}// namespace

#endif
