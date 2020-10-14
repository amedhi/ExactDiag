/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 11:14:38
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include "diag.h"
//#include <nlopt.hpp>
//#include "../optimizer/LBFGS.h"

namespace diag {

ExactDiag::ExactDiag(const input::Parameters& inputs) 
  : graph_(inputs) 
  , model_(inputs, graph_.lattice())
  , basis_(inputs, graph_)
{
	dim_ = basis_.dim();
	Hmat_.resize(dim_,dim_);
	// orbital indices
	int num_sites =  graph_.lattice().num_sites();
	int max_orb = 0;
  for (auto s=graph_.sites_begin(); s!=graph_.sites_end(); ++s) {
  	int norb = graph_.site_dim(s);
  	if (max_orb < norb) max_orb = norb;
  }
  orb_id_.resize(num_sites,max_orb);
  orb_id_.setConstant(-1);
  int n = 0;
  for (auto s=graph_.sites_begin(); s!=graph_.sites_end(); ++s) {
  	int norb = graph_.site_dim(s);
  	for (int i=0; i<norb; ++i) {
  		orb_id_(graph_.site(s),i) = n++;
  	}
  }
  // parameter values
  if (model_.have_parameter("U")) {
    pnames_.push_back("U");
    pvals_.push_back(model_.get_parameter_value("U"));
  }
  if (model_.have_parameter("J")) {
    pnames_.push_back("J");
    pvals_.push_back(model_.get_parameter_value("J"));
  }

  // files
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  if (mode=="APPEND") replace_mode_ = false;
  else replace_mode_ = true;
  std::string prefix = inputs.set_value("prefix", "results");
  boost::algorithm::to_lower(prefix);
  if (prefix=="./"||prefix==""||prefix==".") prefix = "./";
  if (prefix != "./") {
  	prefix = "./"+prefix+"/";
  	boost::filesystem::path prefix_dir(prefix);
  	boost::filesystem::create_directories(prefix_dir);
  }
  make_info_str(inputs);
  dfiles_.push_back(std::make_shared<file::DataFile>());
  dfiles_.back()->init(prefix, "energy",info_str_.str(),replace_mode_);
  dfiles_.push_back(std::make_shared<file::DataFile>());
  dfiles_.back()->init(prefix, "Sx",info_str_.str(),replace_mode_);
  dfiles_.push_back(std::make_shared<file::DataFile>());
  dfiles_.back()->init(prefix, "Sy",info_str_.str(),replace_mode_);
  dfiles_.push_back(std::make_shared<file::DataFile>());
  dfiles_.back()->init(prefix, "Sz",info_str_.str(),replace_mode_);
  dfiles_.push_back(std::make_shared<file::DataFile>());
  dfiles_.back()->init(prefix, "Ssq",info_str_.str(),replace_mode_);
}

int ExactDiag::run(const input::Parameters& inputs) 
{
  //std::cout << "ExactDiag::run\n";
  model_.update_parameters(inputs);
  // parameter values
  double U = 0.0;
  double J = 0.0;
  pvals_.clear();
  if (model_.have_parameter("U")) {
  	U = model_.get_parameter_value("U");
  	pvals_.push_back(U);
  }
  if (model_.have_parameter("J")) {
  	J = model_.get_parameter_value("J");
  	J *= U;
  	if ((U-3*J) < 0.0) {
    	throw std::range_error("ExactDiag::run: out-of-range value of 'J'");
  	}
  	model_.change_parameter_value("J",J);
  	pvals_.push_back(J);
  }
  model_.update_terms();
	construct_hmatrix();
	diagonalize();
  return 0;
}

int ExactDiag::diagonalize(void) 
{
	solver_.compute(Hmat_);
	// pick up ground states & one first excited state
	int nstates = 0;
	for (int i=0; i<dim_; ++i) {
		nstates += 1;
		double E = solver_.eigenvalues()(i);
		if (i > 0 && std::abs(E-solver_.eigenvalues()(i-1))>1.0E-8) break;
	}

	// properties
	energy_.clear();
	Sx_.clear();
	Sy_.clear();
	Sz_.clear();
	Ssq_.clear();
	for (int i=0; i<nstates; ++i) {
		compute_properties(solver_.eigenvectors().col(i));
		energy_.push_back(solver_.eigenvalues()(i));
		Sx_.push_back(props_["Sx"]);
		Sy_.push_back(props_["Sy"]);
		Sz_.push_back(props_["Sz"]);
		Ssq_.push_back(props_["Ssq"]);
	}
	print_out();

	/*
  std::cout<<std::scientific<<std::uppercase<<std::setprecision(6);
	for (int i=0; i<dim_; ++i) {
		double E = solver_.eigenvalues()(i);
		std::cout << "E["<<i<<"] = "<<std::setw(12)<<E<< "\n"; 
		if (i > 0 && std::abs(E-solver_.eigenvalues()(i-1))>1.0E-6) break;
	}
	// eigen states
	int n = 0;
	std::cout << "Groundstate = \n";
	std::cout << "  "<<std::setw(12)<<solver_.eigenvectors()(0,n) 
		<<" "<<basis_.state(0)<<"\n";
	for (int i=1; i<dim_; ++i) {
		double c = solver_.eigenvectors()(i,n);
		if (std::abs(c)>1.0E-4) {
			std::cout << "+ "<<std::setw(14)<<c<<" "<<basis_.state(i)<<"\n";
		}
	}
	*/


  return 0;
}

int ExactDiag::compute_properties(const RealVector&  eigvec)
{
	using namespace model;
	// Sx
	double Sx, Sy, Sz;
	double Splus = 0.0;
	double Sminus = 0.0;
	Sz = 0.0;
  int mat_elem;
	BasisState::idx_t i, j;
  BasisState state;
	for (j=0; j<basis_.dim(); ++j) {
    for (auto s=graph_.sites_begin(); s!=graph_.sites_end(); ++s) {
      int site = graph_.site(s);
      int norb = graph_.site_dim(s);
    	for (int m=0; m<norb; ++m) {
    		// Splus 
    		state = basis_.state(j);
  			mat_elem = state.apply(op::c_dn,orbital(site,m));
  			mat_elem *= state.apply(op::cdag_up,orbital(site,m));
    		if (!basis_.is_null(state)) {
    			i = basis_.idx(state);
    			Splus += mat_elem*eigvec(i)*eigvec(j);
    		} 
    		// Sminus 
    		state = basis_.state(j);
  			mat_elem = state.apply(op::c_up,orbital(site,m));
  			mat_elem *= state.apply(op::cdag_dn,orbital(site,m));
    		if (!basis_.is_null(state)) {
    			i = basis_.idx(state);
    			Sminus += mat_elem*eigvec(i)*eigvec(j);
    		} 
    		// Sz 
    		state = basis_.state(j);
  			mat_elem = state.apply(op::n_up,orbital(site,m));
  			mat_elem -= state.apply(op::n_dn,orbital(site,m));
  			Sz += mat_elem*eigvec(j)*eigvec(j);
    	}
    }
  }
  Sx = 0.5*(Splus+Sminus);
  Sy = 0.5*(Splus-Sminus);
  Sz = 0.5*Sz;

  // Square operators
 	int mat_elem2;
  double SmSp = 0.0;
  double Sz_sq = 0.0;
	for (j=0; j<basis_.dim(); ++j) {
    for (auto s=graph_.sites_begin(); s!=graph_.sites_end(); ++s) {
      int site1 = graph_.site(s);
      int norb1 = graph_.site_dim(s);
    	for (int m=0; m<norb1; ++m) {
    		for (auto s2=graph_.sites_begin(); s2!=graph_.sites_end(); ++s2) {
      		int site2 = graph_.site(s2);
      		int norb2 = graph_.site_dim(s2);
    			for (int n=0; n<norb2; ++n) {
    				// SminusSplus
    				state = basis_.state(j);
  					mat_elem = state.apply(op::c_dn,orbital(site2,n));
  					mat_elem *= state.apply(op::cdag_up,orbital(site2,n));
  					mat_elem *= state.apply(op::c_up,orbital(site1,m));
  					mat_elem *= state.apply(op::cdag_dn,orbital(site1,m));
    				if (!basis_.is_null(state)) {
    					i = basis_.idx(state);
    					SmSp += mat_elem*eigvec(i)*eigvec(j);
    				}

    				// Sz^2
    				state = basis_.state(j);
  					mat_elem = state.apply(op::n_up,orbital(site2,n));
  					mat_elem -= state.apply(op::n_dn,orbital(site2,n));
  					mat_elem2 = state.apply(op::n_up,orbital(site1,m));
  					mat_elem2 -= state.apply(op::n_dn,orbital(site1,m));
  					Sz_sq += mat_elem*mat_elem2*eigvec(j)*eigvec(j);
    			}
    		}
    	}
    }
  }
  Sz_sq *= 0.25;
  double S_sq = SmSp + Sz_sq + Sz;
  props_["Sx"] = Sx;
  props_["Sy"] = Sy;
  props_["Sz"] = Sz;
  props_["Ssq"] = S_sq;
  //std::cout << "S^2, Sx, Sy, Sz = "<<S_sq<<", "<<Sx<<", "<<Sy<<", "<<Sz<< "\n";
	return 0;
}

int ExactDiag::construct_hmatrix(void) 
{
	Hmat_.setZero();
	BasisState::idx_t i, j;
	// site terms
	for (j=0; j<basis_.dim(); ++j) {
    for (auto s=graph_.sites_begin(); s!=graph_.sites_end(); ++s) {
      int site = graph_.site(s);
      int type = graph_.site_type(s);
      int norb = graph_.site_dim(s);
  		for (auto sterm=model_.siteterms_begin(); sterm!=model_.siteterms_end(); ++sterm) {
    		for (int m=0; m<norb; ++m) {
    			for (int n=0; n<norb; ++n) {
    				if (sterm->op_string().act_on(m,n)) {
    					BasisState state = basis_.state(j);
    					int mat_elem = 1;
  						for (const auto& op : sterm->op_string().front_part()) {
  							int p = state.apply(op,orbital(site,n));
  							mat_elem *= p;
  						}
  						for (const auto& op : sterm->op_string().back_part()) {
  							int p = state.apply(op,orbital(site,m));
  							mat_elem *= p;
  						}
    					i = basis_.idx(state);
    					if (!basis_.is_null(i)) {
    						ComplexMatrix coeff_mat = sterm->coupling(type);
    						Hmat_(i,j) += mat_elem * std::real(coeff_mat(m,n));				
    					} 
    				}
    			}
    		}
  		}
    }
	}
	// bond terms
	for (j=0; j<basis_.dim(); ++j) {
    for (auto b=graph_.bonds_begin(); b!=graph_.bonds_end(); ++b) {
    	int s = graph_.source(b);
    	int t = graph_.target(b);
      int sdim = graph_.site_dim(s);
      int tdim = graph_.site_dim(t);
      int type = graph_.bond_type(b);
  		for (auto bterm=model_.bondterms_begin(); bterm!=model_.bondterms_end(); ++bterm) {
    		for (int m=0; m<sdim; ++m) {
    			for (int n=0; n<tdim; ++n) {
    				if (bterm->op_string().act_on(m,n)) {
    					BasisState state = basis_.state(j);
    					int mat_elem = 1;
  						for (const auto& op : bterm->op_string().front_part()) {
  							int p = state.apply(op,orbital(t,n));
  							mat_elem *= p;
  						}
  						for (const auto& op : bterm->op_string().back_part()) {
  							int p = state.apply(op,orbital(s,m));
  							mat_elem *= p;
  						}
    					i = basis_.idx(state);
    					if (!basis_.is_null(i)) {
    						ComplexMatrix coeff_mat = bterm->coupling(type);
    						auto amplitude = static_cast<double>(mat_elem) * coeff_mat(m,n);
    						Hmat_(i,j) += std::real(amplitude);
    						Hmat_(j,i) += std::real(amplitude); // hc part
    					} 
    				}
    			}
    		}
  		}
    }
	}

	// print
	if (false) {
		for (i=0; i<dim_; ++i) {
			std::cout << "s["<<i<<"] = "<<basis_.state(i)<< ": ";
			std::cout << "H["<<i<<","<<i<<"] = "<<Hmat_(i,i) << "\n";
		}
	}	
	/*
	for (i=0; i<dim_; ++i) {
		for (j=0; j<dim_; ++j) {
			std::cout << "H["<<i<<","<<j<<"] = "<<Hmat_(i,j) << "\n";
		}
	}
	*/
	return 0;
}

void ExactDiag::print_out(void)
{

  for (auto& file : dfiles_) {
		file->open();
  	if (!heading_printed_ && replace_mode_) {
    	for (const auto& pname : pnames_) {
      	file->fs()<<std::left<<std::setw(15)<< pname;
    	}
      file->fs()<<"   ";
    	for (int n=0; n<energy_.size(); ++n) {
    		if (file->dname()=="energy")
      		file->fs()<<std::left<<std::setw(15)<<"E"+std::to_string(n);
    		if (file->dname()=="Sx")
      		file->fs()<<std::left<<std::setw(15)<<"Sx"+std::to_string(n);
    		if (file->dname()=="Sy")
      		file->fs()<<std::left<<std::setw(15)<<"Sy"+std::to_string(n);
    		if (file->dname()=="Sz")
      		file->fs()<<std::left<<std::setw(15)<<"Sz"+std::to_string(n);
    		if (file->dname()=="Ssq")
      		file->fs()<<std::left<<std::setw(15)<<"Ssq"+std::to_string(n);
      }
    	file->fs() << "\n";
    	file->fs()<<"#"<< std::string(72, '-') << "\n";
    }
  	file->fs()<<std::scientific<<std::uppercase<<std::setprecision(6);
  	for (const auto& pval : pvals_) {
    	file->fs()<<std::left<<std::setw(15)<<pval;
  	}
    if (file->dname()=="energy") {
  		for (const auto& en : energy_) {
  			file->fs()<<std::right<<std::setw(15)<<en;
  		}
    }
    if (file->dname()=="Sx") {
  		for (const auto& elem : Sx_) {
  			file->fs()<<std::right<<std::setw(15)<<elem;
  		}
    }
    if (file->dname()=="Sy") {
  		for (const auto& elem : Sy_) {
  			file->fs()<<std::right<<std::setw(15)<<elem;
  		}
    }
    if (file->dname()=="Sz") {
  		for (const auto& elem : Sz_) {
  			file->fs()<<std::right<<std::setw(15)<<elem;
  		}
    }
    if (file->dname()=="Ssq") {
  		for (const auto& elem : Ssq_) {
  			file->fs()<<std::right<<std::setw(15)<<elem;
  		}
    }
  	file->fs() << "\n";
		file->close();
  }

  heading_printed_ = true;
}

void ExactDiag::make_info_str(const input::Parameters& inputs)
{
  info_str_.clear();
  print_copyright(info_str_);
  info_str_ << "# "<< inputs.job_id() <<"\n"; 
  info_str_ << model_.info_str(); 
  //info_str_ << "#" << std::string(72, '-') << "\n";
}

void ExactDiag::print_copyright(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: Exact Diagonalization\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}


} // end namespace mc
