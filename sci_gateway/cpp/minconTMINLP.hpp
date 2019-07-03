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

#ifndef minconTMINLP_HPP
#define minconTMINLP_HPP

#include "BonTMINLP.hpp"
#include "IpTNLP.hpp"
#include "api_scilab.h"

using namespace  Ipopt;
using namespace Bonmin;

class minconTMINLP : public TMINLP
{
private:
    
    scilabEnv env_;                //Scilab Environment Variable
    
    Index numVars_;             //Number of variables
    
    Index numCons_;             //Number of constraints
    
    Index numLC_;               //Number of Linear constraints
    
    Index intconSize_;
    
    Number *x0_= NULL;          //lb_ is a pointer to a matrix of size of 1*numVars_ with lower bound of all variables.
    
    Number *lb_= NULL;          //lb_ is a pointer to a matrix of size of 1*numVars_ with lower bound of all variables.
    
    Number *ub_= NULL;          //ub_ is a pointer to a matrix of size of 1*numVars_ with upper bound of all variables.
    
    Number *conLb_= NULL;       //conLb_ is a pointer to a matrix of size of numCon_*1 with lower bound of all constraints.
    
    Number *conUb_= NULL;       //conUb_ is a pointer to a matrix of size of numCon_*1 with upper bound of all constraints.
    
    Number *finalX_= NULL;          //finalX_ is a pointer to a matrix of size of 1*numVars_ with final value for the primal variables.
    
    Number finalObjVal_;        //finalObjVal_ is a scalar with the final value of the objective.
    
    Number *intcon_ = NULL;
    
    int status_;                     //Solver return status
    minconTMINLP(const minconTMINLP&);
    minconTMINLP& operator=(const minconTMINLP&);
    
public:
    // Constructor
    minconTMINLP(scilabEnv env, Index nV, Number *x0, Number *lb, Number *ub, Index nLC, Index nCons, Number *conlb, Number *conub, Index intconSize, Number *intcon): env_(env),  numVars_(nV),x0_(x0),lb_(lb),ub_(ub),numLC_(nLC),numCons_(nCons),conLb_(conlb),conUb_(conub),intconSize_(intconSize),intcon_(intcon),finalX_(0),finalObjVal_(1e20){    }
    
    /** default destructor */
    virtual ~minconTMINLP();
    
    virtual bool get_variables_types(Index n, VariableType* var_types);
    
    virtual bool get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types);
    
    virtual bool get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types);
    
    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);
    
    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u);
    
    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda);
    
    /** Method to return the objective value */
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
    
    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
    
    /** Method to return the constraint residuals */
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
    
    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values);
    
    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    virtual bool eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values);
    
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(SolverReturn status,Index n, const Number* x, Number obj_value);
    
    virtual const SosInfo * sosConstraints() const{return NULL;}
    virtual const BranchingInfo* branchingInfo() const{return NULL;}
    
    virtual bool getHessFromScilab(Index n, wchar_t* name, const Number* x, double obj, const Number* lambda, double * resCh);
    
    const double * getX();        //Returns a pointer to a matrix of size of 1*numVars_
    //with final value for the primal variables.
    
    const double * getGrad();       //Returns a pointer to a matrix of size of 1*numVars_
    //with final value of gradient for the primal variables.
    
    const double * getHess();       //Returns a pointer to a matrix of size of numVars_*numVars_
    //with final value of hessian for the primal variables.
    
    double getObjVal();        //Returns the output of the final value of the objective.
    
    double iterCount();        //Returns the iteration count
    
    int returnStatus();        //Returns the status count
};

#endif
