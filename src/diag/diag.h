/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef EXACTDIAG_H
#define EXACTDIAG_H

#include <iostream>
#include <memory>
#include <map>
#include "../scheduler/worker.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "basis.h"
#include "datafile.h"
#include <boost/filesystem.hpp>

namespace diag {

class ExactDiag : public scheduler::Worker
{
public:
  using idx_t = Basis::idx_t;
  ExactDiag(const input::Parameters& parms); 
  ~ExactDiag() {}
  int start(const input::Parameters& parms) override { return 0; }
  int run(const input::Parameters& parms) override;
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
  //const var::parm_vector& vp(void) { return varparms; }
private:
  lattice::LatticeGraph graph_;
  model::Hamiltonian model_;
  Basis basis_;
  idx_t dim_;
  intMatrix orb_id_;
  RealMatrix Hmat_;
  mutable Eigen::SelfAdjointEigenSolver<RealMatrix> solver_;

  // properties
  std::vector<std::string> pnames_;
  std::vector<double> pvals_;
  std::vector<double> energy_;
  std::vector<double> Sx_;
  std::vector<double> Sy_;
  std::vector<double> Sz_;
  std::vector<double> Ssq_;

  std::map<std::string,double> props_;
  std::vector<std::shared_ptr<file::DataFile>> dfiles_;
  bool heading_printed_{false};
  bool replace_mode_{true};
  mutable std::ostringstream info_str_;

  void make_info_str(const input::Parameters& inputs);
  int construct_hmatrix(void);
  int diagonalize(void);
  int compute_properties(const RealVector&  eigvec);
  void print_out(void);
  const int& orbital(const int& s, const int& m) const { return orb_id_(s,m); }
};

//double wrapper(const std::vector<double>& x, std::vector<double>& grad, void *my_data);

} // end namespace vmc

#endif
