// Copyright (C) 2016 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh, Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#ifndef QuadTMINLP_HPP
#define QuadTMINLP_HPP

#include "BonTMINLP.hpp"
#include "IpTNLP.hpp"

using namespace  Ipopt;
using namespace Bonmin;
    
class QuadTMINLP : public TMINLP
{
	private:
		Index numVars_;			// Number of variables.
	
		Index numCons_; 		// Number of constraints.
		
		Index intconSize_;    // Number of integer constraints

		const Number *qMatrix_ = NULL;	//qMatrix_ is a pointer to matrix of size numVars X numVars_ 
						// with coefficents of quadratic terms in objective function.

		const Number *lMatrix_ = NULL;//lMatrix_ is a pointer to matrix of size 1*numVars_
						// with coefficents of linear terms in objective function.
		
		const Number *intcon_  = NULL; // The matrix containing the integer constraints	
	
		const Number *conMatrix_ = NULL;//conMatrix_ is a pointer to matrix of size numCons X numVars
						// with coefficients of terms in a each objective in each row.

		const Number *conUB_= NULL;	//conUB_ is a pointer to a matrix of size of 1*numCons_
						// with upper bounds of all constraints.

		const Number *conLB_ = NULL;	//conLB_ is a pointer to a matrix of size of 1*numConsn_ 
						// with lower bounds of all constraints.

		const Number *varUB_= NULL;	//varUB_ is a pointer to a matrix of size of 1*numVar_ 
						// with upper bounds of all variables.

		const Number *varLB_= NULL;	//varLB_ is a pointer to a matrix of size of 1*numVar_
						// with lower bounds of all variables.

		const Number *varGuess_= NULL;	//varGuess_ is a pointer to a matrix of size of 1*numVar_
						// with initial guess of all variables.
	
		Number *finalX_= NULL;		//finalX_ is a pointer to a matrix of size of 1*numVar_
						// with final value for the primal variables.

		Number *finalZl_= NULL;		//finalZl_ is a pointer to a matrix of size of 1*numVar_
						// with final values for the lower bound multipliers

		Number *finalZu_= NULL;		//finalZu_ is a pointer to a matrix of size of 1*numVar_
						// with final values for the upper bound multipliers

		Number *finalLambda_= NULL;	//finalLambda_ is a pointer to a matrix of size of 1*numConstr_
						// with final values for the upper bound multipliers

		Number finalObjVal_;		//finalObjVal_ is a scalar with the final value of the objective.

		int status_;			//Solver return status
 
public:
  // Constructor
  QuadTMINLP(Index nV, Index nC, Index intconSize,Number *qM, Number *lM, Number *intcon,Number *cM, Number *cLB, Number *cUB, Number *vLB, Number *vUB,Number *vG):
			numVars_(nV),numCons_(nC),intconSize_(intconSize),qMatrix_(qM),lMatrix_(lM),intcon_(intcon),conMatrix_(cM),conLB_(cLB),conUB_(cUB),varLB_(vLB),varUB_(vUB),varGuess_(vG), finalObjVal_(0){	}
  
  // virtual destructor.
  virtual ~QuadTMINLP(){}

  /* Copy constructor.*/   
  QuadTMINLP(const QuadTMINLP &other){}
  
  // Go to http://coin-or.org/Bonmin for the details of the below methods
  
  virtual bool get_variables_types(Index n, VariableType* var_types);
 
  virtual bool get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types);

  virtual bool get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types);

    
  virtual bool get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,
                            Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);
  
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);
  
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);
  
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);
  
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  virtual void QuadTMINLP::finalize_solution(TMINLP::SolverReturn status,
				Index n, const Number* x,Number obj_value);

  virtual const SosInfo * sosConstraints() const{return NULL;}
  virtual const BranchingInfo* branchingInfo() const{return NULL;}

    // Bonmin methods end here

		virtual const double * getX();		//Returns a pointer to a matrix of size of 1*numVar 
						// with final value for the primal variables.

		virtual double getObjVal();		//Returns the output of the final value of the objective.

		virtual int returnStatus();		//Returns the status count
};

#endif
