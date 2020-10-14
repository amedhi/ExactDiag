/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef HILBERTSPACE_H
#define HILBERTSPACE_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace hs {

class QN 
{
public:
  using value_type = int;
  QN(const std::string& name, const value_type& min=0, const value_type& max=0, 
    const value_type& step=1, const bool& fermionic=true, const bool& half_int=false);
  ~QN() {};
  const value_type& min(void) const { return min_; }
  const value_type& max(void) const { return max_; }
  const value_type& step(void) const { return step_; }
  const int& id(void) const { return id_; }
  const std::string& name(void) const { return name_; }
  const std::size_t& num_states(void) const { return num_states_; }
  const bool& fermionic(void) const { return fermionic_; }
  const bool& half_int(void) const { return half_int_; }
private:
  static int qn_count;
  int id_;
  std::string name_;
  value_type min_;
  value_type max_;
  value_type step_;
  bool fermionic_;
  bool half_int_;
  //bool valid_;
  std::size_t num_states_;
};

class qn_I : public QN
{
public:
  qn_I() : QN("Identity",1,1,0,false,false) {}
};

class qn_Nup : public QN
{
public:
  qn_Nup() : QN("n_up",0,1,1,true,false) {}
};

class qn_Ndn : public QN
{
public:
  qn_Ndn() : QN("n_dn",0,1,1,true,false) {}
};

class qn_Sz : public QN
{
public:
  qn_Sz() : QN("Sz",-1,1,1,false,true) {}
};


class OP 
{
public:
  OP() {}
  OP(const std::string& name, const QN& qn, const int& change, 
     const double& matrix_elem=0);
  ~OP() {}
  const int& id(void) const { return id_; }
  const std::string& name(void) const { return name_; }
  const int& qn_id(void) const { return qn_id_; }
  const int& change(void) const { return qn_change_; }
  const double& matrix_elem(void) const { return matrix_elem_; }
private:
  static int op_count;
  int id_;
  std::string name_;
  int qn_id_;
  int qn_change_;
  double matrix_elem_;
};

class op_I : public OP
{
public:
  op_I() : OP("Identity",qn_I(),0) {}
};

class op_Nup : public OP
{
public:
  op_Nup() : OP("Nup",qn_Nup(),0) {}
};

class op_Ndn : public OP
{
public:
  op_Ndn() : OP("Ndn",qn_Ndn(),0) {}
};

class op_c_up : public OP
{
public:
  op_c_up() : OP("c_up",qn_Nup(),-1,1) {}
};

class op_cdag_up : public OP
{
public:
  op_cdag_up() : OP("cdag_up",qn_Nup(),1,1) {}
};

class op_c_dn : public OP
{
public:
  op_c_dn() : OP("c_dn",qn_Ndn(),-1,1) {}
};

class op_cdag_dn : public OP
{
public:
  op_cdag_dn() : OP("cdag_dn",qn_Ndn(),1,1) {}
};


class site_basis 
{
public:
  site_basis(const std::string& name) : name_{name} { clear(); } 
  ~site_basis() {}
  void clear(void) { qn_list_.clear(); op_list_.clear(); }
  void add_qn(const QN& qn) { qn_list_.push_back(qn); }
  void add_op(const OP& op) { op_list_.push_back(op); }
  const std::string& name(void) const { return name_; }
  const std::vector<QN>& qn_list(void) const { return qn_list_; }
  const std::vector<OP>& op_list(void) const { return op_list_; }
  int num_qn(void) const { return qn_list_.size(); }
  int num_op(void) const { return op_list_.size(); }
private:
  std::string name_;
  std::vector<QN> qn_list_;
  std::vector<OP> op_list_;
};


class op_string : public std::vector<OP>
{
public:
  op_string() {}
  op_string(const OP& op) 
    {
      clear();
      push_back(op);
    }
  op_string(const OP& op1, const OP& op2) 
    {
      clear();
      push_back(op1);
      push_back(op2);
    }
  op_string(const OP& op1, const OP& op2, const OP& op3, const OP& op4) 
    {
      clear();
      push_back(op1);
      push_back(op2);
      push_back(op3);
      push_back(op4);
    }
  ~op_string() {}
};

class fermion_basis : public site_basis
{
public:
  fermion_basis() : site_basis("fermion") 
  {
    add_qn(qn_Nup()); 
    add_qn(qn_Ndn()); 
    add_op(op_Nup()); 
    add_op(op_Ndn()); 
  }
  ~fermion_basis();
};


} // end namespace hs

#endif
