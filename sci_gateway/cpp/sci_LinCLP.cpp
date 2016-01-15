/*
 * Linear Solver Toolbox for Scilab using CLP library
 * Authors :
	Guru Pradeep Reddy
    Bhanu Priya Sayal
*/

#include "LinCLP.hpp"
extern "C"{
#include "api_scilab.h"
#include "Scierror.h"
#include "localization.h"
#include "sciprint.h"
#include "sci_iofunc.hpp"

//creating a problem pointer using Base class of OsiSolverInterface and
//Instantiate the object using specific derived class of ClpSolver
OsiSolverInterface* si = new OsiClpSolverInterface();
//declaration of global variables
int numVars_ , numCons_ , numEqCon_ , numRows_ , numCols_;
double* con;                                          //inequality constraint matrix
double* b;                                            //inequality vector
double* Aeq;                                          //equality constraint matrix
double* beq;                                          //equality vector
const double* xValue;                                 //

//Clp Solver function definition
LinCLP::LinCLP(int numVars_ , int numCons_ , int numEqCons_ , double objMatrix_[] , double conMatrix_[] , double bMatrix_[] , double AeqMatrix_[] , double beqMatrix_[] , double conLB_[] , double conUB_[], double options_[], double flagMps_, double flagLp_, char mpsFile_[], char lpFile_[])
{
   numCols_  = numVars_;                              //assigning number of colums = number of variables
   numRows_  = numCons_;                              //assigning number of rows = number of inequality constraints
   numEqCon_ = numEqCons_;                            //number of equality constraints
   double *objective = new double[numCols_];          //the objective coefficients
   double *colLB_    = new double[numCols_];          //the lower bounds for variables
   double *colUB_    = new double[numCols_];          //the upper bounds for variables
   //allocating memory to the inequality constraint matrix
   con = new double[numCols_*numRows_]; 
   for(int i =0; i <numCols_*numRows_;i++)
   	con[i] = conMatrix_[i];
   //allocating memory to the inequality vector
   b = new double[numRows_];
   for(int i =0; i <numRows_;i++)
   	b[i] = bMatrix_[i];
   //allocating memory to the equality constraint matrix
   Aeq = new double[numCols_*numEqCons_];
   for(int i =0; i <numCols_*numEqCons_;i++)
   	Aeq[i] = AeqMatrix_[i];
   //allocating memory to the equality vector 
   beq = new double[numEqCons_];
   	for(int i =0; i <numEqCons_;i++)
   		beq[i] = beqMatrix_[i];

   //Defining the objective coefficients
   for(int i=0 ; i<numCols_ ; i++)
   	{
		objective[i] = objMatrix_[i];
   	}
   
   //Defining the variable lower/upper bounds
   for(int i=0 ; i<numCols_ ; i++)
   	{
   		colLB_[i] = conLB_[i];
   		colUB_[i] = conUB_[i];
   	}
  
   double *rowLB_ = new double[numRows_];             //the row lower bounds
   double *rowUB_ = new double[numRows_];             //the row upper bounds

   //Defining the constraint matrix
   CoinPackedMatrix *matrix =  new CoinPackedMatrix(false , 0 , 0);
   matrix->setDimensions(0 , numCols_);
   for(int i=0 ; i<numCons_ ; i++)
   	{
    	CoinPackedVector row;
	 	for(int j=0 ; j<numVars_ ; j++)
	 		{
   				row.insert(j, conMatrix_[i+j*numCons_]);
   	 		}

   	 	rowLB_[i] = -si->getInfinity();       	      //Setting the row lower bound 
   	 	rowUB_[i] = bMatrix_[i];               		  //Setting the row upper bound 
        matrix->appendRow(row);
   	}

   //Adding equality constraints
   for(int i=0 ; i<numEqCons_ ; i++)
   	{
    	CoinPackedVector row;
       	for(int j=0 ; j<numVars_ ; j++)
        	{
        		row.insert(j , AeqMatrix_[i+j*numEqCons_]);
			}

        rowLB_[numCons_+i] = beqMatrix_[i];           //Setting the row lower bound  
        rowUB_[numCons_+i] = beqMatrix_[i];           //Setting the row upper bound
        matrix->appendRow(row);
   	} 

   //setting options for maximum iterations
   si->setIntParam(OsiMaxNumIteration,options_[0]);

   //Load the problem to OSI
   si->loadProblem(*matrix , colLB_ , colUB_ , objective , rowLB_ , rowUB_);
    
   //writing mps if flag is set
   if(flagMps_ == 1)
   	si->writeMps(mpsFile_);

   //writing lp if flag is set
   if(flagLp_ == 1)
   	si->writeLp(lpFile_);

   //Solve the problem
   si->initialSolve();
  
}

   //Output the solution to Scilab
   //get solution for x
   const double* LinCLP::getX()
	{ 
        xValue = si->getColSolution();
		return xValue;
	}

   //get objective value
   double LinCLP::getObjVal()
	{
		const double objValue = si->getObjValue();
		return objValue;
	}
   
   //get exit status 
   int LinCLP::returnStatus()
	{
   		double status;
   		if(si->isProvenOptimal())
    			status=0;
   		else if(si->isProvenPrimalInfeasible())
        		status=1;
   		else if(si->isProvenDualInfeasible())
        		status=2;
   		else if(si->isIterationLimitReached())
        		status=3;
   		else if(si->isAbandoned())
        		status=4;
   		else if(si->isPrimalObjectiveLimitReached())
        		status=5;
   		else if(si->isDualObjectiveLimitReached())
        		status=6;

		return status;
	}

   //get constraint violation
   double LinCLP::consViolation()
	{ 
        double maxViolation = 0.0;
   		double sum = 0;
   		double temp;
		for(int i=0 ; i<numRows_ ; i++)
			{  
				sum=0;
				for(int j=0 ; j<numCols_ ; j++)
    					{
				      		sum = sum+con[i+j*numRows_]*xValue[j];
    					} 
    				temp = sum-b[i];
    				if(maxViolation<temp)
       					{
							maxViolation = temp;
                         }
   			}
   		for(int i=0 ; i<numEqCon_ ; i++)
   			{
				sum = 0;
				for(int j=0 ; j<numCols_ ; j++)
    					{
      						sum = sum+Aeq[i+j*numEqCon_]*xValue[j];
   	    				} 
    				if(sum-beq[i]<0)
    					temp = -(sum-beq[i]);
    				else
    					temp = sum-beq[i];
    				if(maxViolation<temp)
       					maxViolation = temp;
   			}
		return maxViolation;

	}

   //get number of iterations
   double LinCLP::iterCount()
	{
		const double iterations = si->getIterationCount(); 
		return iterations;
	}

   //get lower vector
   const double* LinCLP::getReducedCost()
	{
		const double* reducedCost = si->getReducedCost();
		return reducedCost;
	}

   //get inequality dual vector
   double* LinCLP::getInEqlin()
	{
		double* dual = si->getRowPrice();
        double* ineqlin;
        ineqlin = new double[numRows_];
        for(int j=0 ; j<numRows_ ; j++)
	       	{
            	ineqlin[j]=dual[j];
            }
        return ineqlin;
	}

   //get equality dual vector
    double* LinCLP::getEqlin()
	{
		const double* dual = si->getRowPrice();
        double* eqlin;
        eqlin = new double[numEqCons_];
        for(int i=numRows_ ; i<numRows_+numEqCons_ ; i++)
			{
        		eqlin[i-numRows_]=dual[i];
			}	
		return eqlin;            
	}   
}
   
