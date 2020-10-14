/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-11 13:02:35
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-17 12:02:29
*----------------------------------------------------------------------------*/
#include <cmath>
#include "model.h"
#include <boost/algorithm/string.hpp>

namespace model {

int Hamiltonian::define_model(const input::Parameters& inputs, 
  const lattice::Lattice& lattice)
{
  double defval;
  std::string name, path; //, matrixelem, op, qn, site, src, tgt, fact;
  CouplingConstant cc;

  // define the models 
  model_name = inputs.set_value("model", "HUBBARD");
  boost::to_upper(model_name);

  int num_bands, num_orb;
  strMatrix expr_mat;
  strMatrix::row_t expr_vec;

  if (model_name == "HUBBARD") {
    mid = model_id::HUBBARD;
    if (lattice.id()==lattice::lattice_id::SQUARE) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      add_bondterm(name="hop_UU", cc="-t", hop_UU(1));
      add_bondterm(name="hop_DD", cc="-t", hop_DD(1));
      // site terms
      cc.create(1);
      expr_mat.resize(1,1);
      expr_mat(0,0) = "U";
      cc.add_type(0, expr_mat);
      add_siteterm(name="IntraOrb_UD", cc, IntraOrb_UD(1));
    }
    else if (lattice.id()==lattice::lattice_id::CHAIN) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);

      // hopping term
      add_bondterm(name="hop_UU", cc="-t", hop_UU(1));
      add_bondterm(name="hop_DD", cc="-t", hop_DD(1));
      // site terms
      cc.create(1);
      expr_mat.resize(1,1);
      expr_mat(0,0) = "U";
      cc.add_type(0, expr_mat);
      add_siteterm(name="IntraOrb_UD", cc, IntraOrb_UD(1));
    }
    else {
      throw std::range_error("*error: modellibrary: model not defined for the lattice"); 
    }
  }

  else if (model_name == "HUBBARD_NBAND") {
    mid = model_id::HUBBARD_NBAND;
    if (lattice.id()==lattice::lattice_id::SQUARE_NBAND ||
        lattice.id()==lattice::lattice_id::CHAIN_NBAND ) {
      num_bands = inputs.set_value("num_bands", 1);
      num_orb = num_bands;
      add_parameter("num_bands", num_bands);
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      add_parameter(name="J", defval=0.0, inputs);
      //double J = get_parameter_value("J")*get_parameter_value("U");
      //change_parameter_value("J",J);

      // bond operators (Diagonal in either product basis or SOC basis)
      cc.create(1);
      expr_mat.resize(num_orb,num_orb);
      for (int i=0; i<num_orb; ++i) {
        for (int j=0; j<num_orb; ++j) {
          if (i==j) expr_mat(i,j) = "-t";
          else expr_mat(i,j) = "0";
        }
      }
      cc.add_type(0, expr_mat);
      add_bondterm(name="hop_UU", cc, hop_UU(num_orb));
      add_bondterm(name="hop_DD", cc, hop_DD(num_orb));

      // Interaction
      add_siteterm(name="IntraOrb_UD", cc=strMatrix("U",num_orb), IntraOrb_UD(num_orb));
      add_siteterm(name="InterOrb_UD", cc=strMatrix("U-2*J",num_orb), InterOrb_UD(num_orb));
      add_siteterm(name="InterOrb_UU", cc=strMatrix("U-3*J",num_orb), InterOrb_UU(num_orb));
      add_siteterm(name="InterOrb_DD", cc=strMatrix("U-3*J",num_orb), InterOrb_DD(num_orb));
      add_siteterm(name="SpinFlip", cc=strMatrix("-J",num_orb), SpinFlip(num_orb));
      add_siteterm(name="PairHop", cc=strMatrix("J",num_orb), PairHop(num_orb));
      // Define standard operators
    }
    else {
      throw std::range_error("*error: modellibrary: model not defined for the lattice"); 
    }
  }

  /*------------- undefined model--------------*/
  else {
    throw std::range_error("*error: modellibrary: undefined model");
  }

  // if the model has site disorder
  /*
  if (site_disorder) {
    add_disorder_term(name="disorder", op::ni_sigma());
  }*/
  
  return 0;
}

int Hamiltonian::construct(const input::Parameters& inputs, 
  const lattice::Lattice& lattice)
{
  init(lattice);
  define_model(inputs, lattice);
  finalize(lattice);
  return 0;
}


} // end namespace model
