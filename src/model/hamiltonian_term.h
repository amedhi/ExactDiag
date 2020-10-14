/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef HAMILTONIAN_TERM_H
#define HAMILTONIAN_TERM_H

#include <string>
#include <vector>
#include <complex>
#include <array>
#include <map>
//#include <unordered_map>
#include <stdexcept>
#include <Eigen/Core>
#include "modelparams.h"
#include "strmatrix.h"
#include "../diag/matrix.h"
#include "../lattice/lattice.h"
//#include "./quantum_op.h"
//#include "./hilbertspace.h"

namespace model {

enum class op {
  n_up, n_dn, n, n_ud, cdag_up, c_up, cdag_dn, c_dn, hop_up, hop_dn, null
};

class op_string : public std::vector<op>
{
public:
  op_string() {}
  op_string(const op& op) 
    {
      clear();
      push_back(op);
    }
  op_string(const op& op1, const op& op2) 
    {
      clear();
      push_back(op1);
      push_back(op2);
    }
  op_string(const op& op1, const op& op2, const op& op3, const op& op4) 
    {
      clear();
      push_back(op1);
      push_back(op2);
      push_back(op3);
      push_back(op4);
    }
  ~op_string() {}
};

class TermOp 
{
public:
  TermOp() {}
  TermOp(const std::string& name, const std::pair<op_string,op_string>& op_pair);
  TermOp(const std::string& name, const std::pair<op_string,op_string>& op_pair,
    const intMatrix& action_mat);
  ~TermOp() {}
  void set_action_matrix(const intMatrix& action_mat) 
  { action_matrix_=action_mat; }
  const std::string& name(void) const { return name_; }
  const op_string& front_part(void) const { return front_part_; }
  const op_string& back_part(void) const { return back_part_; }
  const intMatrix& action_matrix(void) const { return action_matrix_; }
  bool act_on(const int& i, const int& j) const { return action_matrix_(i,j); }
protected:
  int site_dim_;
  intMatrix action_matrix_;
private:
  std::string name_;
  op_string back_part_;
  op_string front_part_;
  bool is_diagonal_{true};
};


class IntraOrb_UD : public TermOp
{
public:
  IntraOrb_UD(const int& num_orbitals);
  ~IntraOrb_UD() {}
};

class InterOrb_UD : public TermOp
{
public:
  InterOrb_UD(const int& num_orbitals);
  ~InterOrb_UD() {}
};

class InterOrb_UU : public TermOp
{
public:
  InterOrb_UU(const int& num_orbitals);
  ~InterOrb_UU() {}
};

class InterOrb_DD : public TermOp
{
public:
  InterOrb_DD(const int& num_orbitals);
  ~InterOrb_DD() {}
};

class SpinFlip : public TermOp
{
public:
  SpinFlip(const int& num_orbitals);
  ~SpinFlip() {}
};

class PairHop : public TermOp
{
public:
  PairHop(const int& num_orbitals);
  ~PairHop() {}
};

class hop_UU : public TermOp
{
public:
  hop_UU(const int& num_orbitals);
  ~hop_UU() {}
};

class hop_DD : public TermOp
{
public:
  hop_DD(const int& num_orbitals);
  ~hop_DD() {}
};

//---------------------------------------------------------------
class CouplingConstant : public std::unordered_map<int, strMatrix>
{
public:
  using super_type = std::unordered_map<int,strMatrix>;
  using iterator = super_type::iterator;
  using const_iterator = super_type::const_iterator;
  using value_type = super_type::value_type;

  CouplingConstant() {}
  CouplingConstant(const std::string& expr); 
  CouplingConstant(const strMatrix::row_t& expr_vec);
  ~CouplingConstant() {}
  CouplingConstant& operator=(const std::string& expr); 
  CouplingConstant& operator=(const std::vector<std::string>& expr_vec); 
  CouplingConstant& operator=(const strMatrix& expr_mat); 
  //std::pair<iterator, bool> insert(const value_type& val);
  void create(const unsigned& num_type);
  void add_type(const unsigned& type, const std::string& expr);
  void add_type(const unsigned& type, const std::vector<std::string>& expr_vec);
  void add_type(const unsigned& type, const strMatrix& expr_mat);
  //void create(const value_type& type0, const value_type& type1={0,"_null_"}, 
  //  const value_type& type2={0,"_null_"}, const value_type& type3={0,"_null_"}, 
  //  const value_type& type4={0,"_null_"}, const value_type& type5={0,"_null_"});
  //void add_type(const value_type& val);
  void clear(void); 
  void clear_map(void) { super_type::clear(); } 
  const bool& valid(void) const { return valid_; } 
  const std::string& expression(const unsigned& type) const;
  static const int global_type;
private:
  int num_types_{0}; 
  bool valid_{false};
};

class HamiltonianTerm 
{
public:
  using ccval_t = ComplexMatrix;
  //using BondSiteMap = std::map<unsigned, std::pair<unsigned, unsigned> >;
  HamiltonianTerm() {}
  HamiltonianTerm(const std::string& name, const TermOp& op, const CouplingConstant& cc, 
    const int& size);
  ~HamiltonianTerm() {}
  void init();
  void eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals);
  const TermOp& op_string(void) const { return op_; }
  bool is_defined_for(const unsigned& operand_type) const 
    { return is_defined_[operand_type]; }
  const ccval_t& coupling(const unsigned& operand_type) const
    { return cc_values_[operand_type]; }
  const strMatrix& coupling_expr(const unsigned& operand_type) const
    { return cc_.at(operand_type); }
  const std::string& name(void) const { return name_; }
private:
  std::string name_;
  TermOp op_;
  CouplingConstant cc_;
  int max_operand_types_;
  std::vector<bool> is_defined_; // operator defined for a given operand 'type'
  std::vector<ccval_t> cc_values_;
};


} // end namespace model

#endif
