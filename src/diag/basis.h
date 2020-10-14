/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef BASIS_H
#define BASIS_H

#include <Eigen/Sparse>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <utility>
#include "../lattice/graph.h"
#include "../model/hamiltonian_term.h"
//#include "qn.h"
//#include "quantum_operator.h"

namespace diag {

class BasisState 
{
public:
  using idx_t = unsigned long long;
  using op_output = std::pair<int,BasisState>;
  enum {max_bits=std::numeric_limits<idx_t>::digits};
  using bitstring = std::bitset<max_bits>;
  BasisState(); 
  BasisState(const bitstring& uket, const bitstring& dket); 
  ~BasisState() {}
  static int set_counts(const int& norb, const int& ntotal, const bool& n_consvd);
  static int set_counts(const int& norb, const int& nup, const int& ndn, const bool& n_consvd);
  int reset(void); 
  int init_state(void); 
  int shift_next(void);
  int apply(const model::op& site_op, const int& i);
  int op_n_up(const int& i) const;
  int op_n_dn(const int& i) const;
  int op_cdag_up(const int& i);
  int op_cdag_dn(const int& i);
  int op_c_up(const int& i);
  int op_c_dn(const int& i);
  int op_hop_up(const int& fr, const int& to);
  int op_hop_dn(const int& fr, const int& to);
  std::pair<idx_t,idx_t> id(void) const; 
  static std::pair<idx_t,idx_t> max_id(void); 
  friend std::ostream& operator<<(std::ostream& os, const BasisState& bs);
private:
  bitstring uket_;
  bitstring dket_;
  static int Norb;
  static int Ntotal;
  static int Nup;
  static int Ndn;
  static bool N_Consvd;
  static bool Sz_Consvd;
  int shift_next(bitstring& bs);
};

class Basis : private std::vector<BasisState>
{
public:
  using idx_t = BasisState::idx_t;
  Basis(const input::Parameters& parms, const lattice::LatticeGraph& graph);
  ~Basis() {}
  const idx_t& dim(void) const { return dim_; }
  const BasisState& state(const idx_t& i) const { return operator[](i); }
  const idx_t& idx(const BasisState& state) const 
    { return idx_map_[ivalue(state.id())]; };
  bool is_null(const BasisState& state) const 
    { 
      if (ivalue(state.id())>=idx_map_.size()) return true;
      else return  idx_map_[ivalue(state.id())]==null_idx_; 
    }
  bool is_null(const idx_t& idx) const { return idx==null_idx_; }
private:
  int num_sites_;
  int num_orbitals_;
  //int num_upspins_;
  //int num_dnspins_;
  idx_t dim_;
  idx_t maxid_up_;
  idx_t maxid_dn_;
  idx_t null_idx_;
  std::vector<idx_t> idx_map_;
  int combination(const int& n, const int& k);
  inline idx_t ivalue(const std::pair<idx_t,idx_t>& id) const 
    { return id.first*maxid_dn_+id.second; }
};



} // end namespace diag

#endif
