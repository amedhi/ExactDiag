/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include <tuple>
#include <algorithm>
#include <stdexcept>
#include "basis.h"

namespace diag {

int BasisState::Norb = 0;
int BasisState::Ntotal = 0;
int BasisState::Nup = 0;
int BasisState::Ndn = 0;
bool BasisState::N_Consvd = true;
bool BasisState::Sz_Consvd = true;

BasisState::BasisState()
{
  uket_.reset();
  dket_.reset();
}

BasisState::BasisState(const bitstring& uket, const bitstring& dket)
  : uket_{uket}, dket_{dket}
{
}

int BasisState::set_counts(const int& norb, const int& ntotal, const bool& n_consvd)
{
  if (norb < 0) {
    throw std::range_error("BasisState::set_counts: invalid argument");
  }
  if (ntotal<0 || ntotal>2*norb) {
    throw std::range_error("BasisState::set_counts: out-of-range particle number");
  } 
  Norb = norb;
  Ntotal = ntotal;
  if (Ntotal > Norb) {
    Ndn = Norb;
    Nup = Ntotal-Ndn;
  }
  else {
    Ndn = Ntotal;
    Nup = 0;
  }
  N_Consvd = n_consvd;
  Sz_Consvd = false;
  return 0;
}

int BasisState::set_counts(const int& norb, const int& nup, const int& ndn, const bool& n_consvd)
{
  if (norb<0 || nup<0 || ndn<0) {
    throw std::range_error("BasisState::set_counts: invalid argument");
  }
  if (nup>norb || ndn>norb) {
    throw std::range_error("BasisState::set_counts: out-of-range particle number");
  }
  if ((nup+ndn)>2*norb) {
    throw std::range_error("BasisState::set_counts: out-of-range particle number");
  } 
  Norb = norb;
  Nup = nup;
  Ndn = ndn;
  Ntotal = nup + ndn;
  N_Consvd = n_consvd;
  Sz_Consvd = true;
  return 0;
}


std::ostream& operator<<(std::ostream& os, const BasisState& bs)
{
  os << "|";
  for (int i=bs.Norb-1; i>=0; --i) {
    os << bs.uket_[i];
  }
  os << " ";
  for (int i=bs.Norb-1; i>=0; --i) {
    os << bs.dket_[i];
  }
  os << ">";
  return os;
}

int BasisState::reset(void) 
{
  uket_.reset();
  dket_.reset();
  return 0;
}

int BasisState::init_state(void) 
{
  uket_.reset();
  dket_.reset();
  for (int i=0; i<Nup; ++i) uket_.set(i);
  for (int i=0; i<Ndn; ++i) dket_.set(i);
  return 0;
}

int BasisState::shift_next(void) 
{
  int m = shift_next(dket_);
  if (m > 0) return m;
  int n = shift_next(uket_);
  if (n > 0) return n;
  // higher Sz sector
  int n_up = uket_.count();
  int n_dn = dket_.count();
  if (n_up < Ntotal && n_up < Norb && n_dn>0) {
    BasisState::Nup += 1;
    BasisState::Ndn -= 1;
    init_state();
    return 1;
  }
  else {
    set_counts(Norb, Ntotal, true);
    init_state();
    return 0;
  }
}

int BasisState::shift_next(bitstring& ket) 
{
  int n = ket.count();
  if (n == Norb) return 0;
  if (n > 0) {
    bitstring bs;
    for (int i=Norb-n; i<Norb; ++i) bs.set(i);
    idx_t imax = bs.to_ullong();
    idx_t istate = ket.to_ullong();
    while (true) {
      if (istate >= imax) {
        ket.reset();
        for (int i=0; i<n; ++i) ket.set(i);
        return -1;
      }
      ket = bitstring(++istate);
      int m = 0;
      for (int i=Norb-1; i>=0; --i) m += ket[i];
      if (m == n) return 1;
    }
  }
  return 0;
}


std::pair<BasisState::idx_t,BasisState::idx_t> BasisState::id(void) const
{
  return std::make_pair(uket_.to_ullong(),dket_.to_ullong());
}

std::pair<BasisState::idx_t,BasisState::idx_t> BasisState::max_id(void) 
{
  bitstring up, dn;
  for (int i=Norb-Nup; i<Norb; ++i) up.set(i);
  for (int i=Norb-Ndn; i<Norb; ++i) dn.set(i);
  idx_t id = std::max(up.to_ullong(),dn.to_ullong());
  return std::make_pair(id,id);
}

int BasisState::apply(const model::op& site_op, const int& i)
{
  switch (site_op) {
    case model::op::n_up: return op_n_up(i);
    case model::op::n_dn: return op_n_dn(i);
    case model::op::n: return (op_n_up(i)+op_n_dn(i));
    case model::op::n_ud: return (op_n_up(i)*op_n_dn(i));
    case model::op::cdag_up: return op_cdag_up(i);
    case model::op::cdag_dn: return op_cdag_dn(i);
    case model::op::c_up: return op_c_up(i);
    case model::op::c_dn: return op_c_dn(i);
    //case model::op::hop_dn: return op_hop_up(i);
    //case model::op::hop_dn: return op_hop_dn(i);
    default: throw std::range_error("BasisState::apply: unknown site operator");
  }
}

int BasisState::op_n_up(const int& i) const
{
  return uket_.test(i);
}

int BasisState::op_n_dn(const int& i) const
{
  return dket_.test(i);
}

int BasisState::op_cdag_up(const int& i)
{
  if (uket_.test(i)) {
    return 0;
  }
  else {
    int sign = 1;
    for (int b=i+1; b<Norb; ++b) if (uket_.test(b)) sign = -sign;
    uket_.set(i);
    return sign;
  }
}

int BasisState::op_cdag_dn(const int& i) 
{
  if (dket_.test(i)) {
    return 0;
  }
  else {
    int sign = 1;
    for (int b=0; b<Norb; ++b) if (uket_.test(b)) sign = -sign;
    for (int b=i+1; b<Norb; ++b) if (dket_.test(b)) sign = -sign;
    dket_.set(i);
    return sign;
  }
}

int BasisState::op_c_up(const int& i) 
{
  if (uket_.test(i)) {
    int sign = 1;
    for (int b=i+1; b<Norb; ++b) if (uket_.test(b)) sign = -sign;
    uket_.reset(i);
    return sign;
  }
  else {
    return 0;
  }
}

int BasisState::op_c_dn(const int& i) 
{
  if (dket_.test(i)) {
    dket_.reset(i);
    int sign = 1;
    for (int b=0; b<Norb; ++b) if (uket_.test(b)) sign = -sign;
    for (int b=i+1; b<Norb; ++b) if (dket_.test(b)) sign = -sign;
    return sign;
  }
  else {
    return 0;
  }
}

int BasisState::op_hop_up(const int& fr, const int& to)
{
  int n = op_c_up(fr); 
  if (n == 0) return n;
  else return n*op_cdag_up(to);
}

int BasisState::op_hop_dn(const int& fr, const int& to)
{
  int n = op_c_dn(fr); 
  if (n == 0) return n;
  else return n*op_cdag_dn(to);
}

//----------------------------------------------------------------------------
Basis::Basis(const input::Parameters& parms, const lattice::LatticeGraph& graph)
{
  num_sites_ = graph.num_sites();
  num_orbitals_ = graph.lattice().num_unitcells()
                   *graph.lattice().num_basis_orbitals();
  int particles_per_site = parms.set_value("particles_per_site",1);
  int ntotal = particles_per_site*num_sites_;
  std::cout << "(sites, orbitals, particles) = " << num_sites_ <<", "<<
    num_orbitals_<<", "<<ntotal<<"\n"; 
  BasisState::set_counts(num_orbitals_,ntotal,true);
  // assuming 'ntotal' conserved
  dim_ = 0;
  int nmax = std::min(ntotal,num_orbitals_);
  for (int n=0; n<=nmax; ++n) {
    int nup = n;
    int ndn = ntotal-n;
    int sdim = combination(num_orbitals_,nup)*combination(num_orbitals_,ndn);
    dim_ += sdim;
    //std::cout << "Basis::dimension = "<<nup<<" "<<ndn<<" "<< sdim << "\n";
  }
  null_idx_ = dim_;
  //std::cout << "Basis::dimension = " << dim_ << "\n";
  //getchar();

  // max id of DN states
  auto id = BasisState::max_id();
  maxid_up_ = id.first;
  maxid_dn_ = id.second;
  idx_t max_id = maxid_up_*maxid_dn_+maxid_dn_;
  //std::cout << "max_id " << max_id << "\n";
  //getchar();

  // generate and store the states
  resize(dim_);
  idx_map_.resize(max_id+1);
  for (auto& idx : idx_map_) idx = null_idx_;
  // initial state
  BasisState state;
  state.init_state();
  operator[](0) = state;
  idx_map_[ivalue(state.id())] = 0;
  for (idx_t i=1; i<dim_; ++i) {
    state.shift_next();
    operator[](i) = state;
    idx_map_[ivalue(state.id())] = i;
    //std::cout << "s["<<i<<"] = "<<state << ": ";
    //std::cout << "("<<state.id().first<<", "<<state.id().second<<")"<< 
    //" = " << ivalue(state.id()) << "\n\n";
  }
  /* 
  for (idx_t i=0; i<dim_; ++i) {
    std::cout << "s["<<i<<"] = "<< (*this)[i] << ": ";
    idx_t idx = idx_map_[ivalue((*this)[i].id())];
    std::cout << "idx = "<< idx << "\n";
    if (i != idx) {
      throw std::logic_error("BasisState: incorrect indexing");
    }
  }
  getchar();
  /*/
}


int Basis::combination(const int& n, const int& k) {
  int m = k;
  // use symmetry property C(n,k)=C(n, n-k)
  if ( m > n-m ) m = n-m;
  if (n < 0 || m<0) return 0;
  else if (n < m) return 0;
  else if ( m == n || m == 0 ) return 1;
  else return combination(n-1,m-1)+combination(n-1,m);
}

} // end namespace diag

