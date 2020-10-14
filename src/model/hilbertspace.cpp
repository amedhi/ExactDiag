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
#include "hilbertspace.h"

namespace hs {

int QN::qn_count = 0;
int OP::op_count = 0;

QN::QN(const std::string& name, const value_type& min, const value_type& max, 
  const value_type& step, const bool& fermionic, const bool& half_int) 
  : name_{name}, min_{min}, max_{max}, step_{step}, fermionic_{fermionic},
    half_int_{half_int}
{
  if (step<0 || min_>max_)  
    throw std::invalid_argument("QN::QN: invalid argument");
  if (step_==0 && min_ != max_) 
    throw std::invalid_argument("QN::QN: invalid argument");
  if (max_>min_ && (max_-min_)%step_ !=0)
    throw std::invalid_argument("QN::QN: invalid argument");
  if (min_ == max_) num_states_ = 1;
  else num_states_ = (max_-min_)/step_ + 1;
  if (half_int_) {
    step_ *= 2;
    num_states_ = (max_-min_)/step_ + 1;
  }
  id_ = QN::qn_count++;
}

OP::OP(const std::string& name, const QN& qn, const int& change, 
     const double& matrix_elem)
  : name_{name}, qn_id_{qn.id()}, qn_change_{change}, matrix_elem_{matrix_elem}
{
}



} // end namespace hs
