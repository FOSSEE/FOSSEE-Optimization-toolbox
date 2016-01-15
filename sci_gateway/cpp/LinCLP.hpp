/*
 * Linear Solver Toolbox for Scilab using CLP library
 * Authors :
	Guru Pradeep Reddy
        Bhanu Priya Sayal

* Optimizing (minimizing) the linear objective function having any number of 
  variables and linear constraints(equality/inequality).
 *
*/

#ifndef __LinCLP_HPP__
#define __LinCLP_HPP__

#include"OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"

class LinCLP
{
	private:

		int numVars_;				//Number of variables

		int numCons_;          		//Number of inequality constraints

		int numEqCons_;        		//Number of equality constraints

		double objMatrix_[];   		//Objective function vector 
 
		double conMatrix_[];   		//Inequality constraint matrix
  
		double bMatrix_[];     		//Inequality constraint vector

		double AeqMatrix_[];   		//Equality constraint matrix

		double beqMatrix_[];   		//Equality constraint vector
 
		double conLB_[];       		//Lower bounds for all variables

		double conUB_[];       		//Upper bounds for all variables

        double options_[];          //options for setting maximum iterations and writing mps and lp files

        double flagMps_;            //if flag is set, then write into mps

        double flagLp_;             //if flag is set, then write into lp
      
        char mpsPath_[];            //File path to store mps file

        char lpPath_[];             //File path to store lp file



	public:
/*
 * Constructor 
*/
		LinCLP(int numVars_ , int numCons_ , int numEqCons , double objMatrix_[] , double conMatrix_[] , double bMatrix_[] , double AeqMatrix_[] , double beqMatrix_[] , double conLB_[] , double conUB_[], double options_[], double flagMps_, double flagLp_, char mpsPath_[],char lpPath_[]);
		
		const double* getX();   	//Returns a pointer to matrix of size 
									//1*numVars with final values for the objective variables

	    double getObjVal();     	//Returns the output of the final value of the objective

		int returnStatus();     	//Returns the status of the problem 
           
        double consViolation(); 	//Returns maximum of constraint violation 

	    double iterCount();     	//Returns the iteration count

        const double* getReducedCost();   //Returns a pointer to matrix of size 
									      //1*numVars with values for lower dual vector

	    double* getInEqlin();   	//Returns a pointer to matrix of size
									//1*numCons with values for inequality dual vector
  
        double* getEqlin();  		//Returns a pointer to matrix of size 
									//1*numEqCons with values for equality dual vector

};

#endif __LinCLP_HPP__

