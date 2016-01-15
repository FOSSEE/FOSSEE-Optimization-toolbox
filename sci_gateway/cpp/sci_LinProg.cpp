/*
 * Linear Solver Toolbox for Scilab using CLP library
 * Authors :
	Guru Pradeep Reddy
    Bhanu Priya Sayal
*/

#include "sci_iofunc.hpp"
#include "LinCLP.hpp"

extern "C"{
#include <api_scilab.h>
#include <Scierror.h>
#include <localization.h>
#include <sciprint.h>

//Solver function
int sci_linearprog(char *fname) 
{
  //Objective function
  double* obj_;  
  //Inequality constraint matrix coefficients
  double* A_;  
  //Inequality vector 
  double* b_;  
  //Equality constraint matrix coefficients
  double* Aeq_;  
  //Equality vector
  double* beq_;  
  //Lower bounds for variables
  double* conLB_;  
  //Upper bounds for variables
  double* conUB_;
  //options for maximum iterations and writing mps
  double* options_;
  //mps flag value
  double flagMps_;
  //lp flag value
  double flagLp_;
  //mps file path
  char * mpsFile_;
  //lp file path
  char * lpFile_;
  //Error structure in Scilab  
  SciErr sciErr;
  //Number of rows and columns in objective function
  int objRows_, objCols_;
        
  CheckInputArgument(pvApiCtx , 12 , 12);             //Checking the input arguments
  CheckOutputArgument(pvApiCtx , 8, 8);               //Checking the output arguments

  //Getting the input arguments from Scilab

  //Objective function from Scilab
  if(getDoubleMatrixFromScilab(1 , &objRows_ , &objCols_ , &obj_))
  	{
		return 1;		
  	}  
  int numVars_ = objCols_;                            //Assigning number of variables = number of columns in objective function
  
  //Inequality constraint matrix from Scilab
  int ARows_;                                         //Number of rows in Inequality constraint matrix
  int ACols_;                                         //Number of columns in Inequality constraint matrix
  if(getDoubleMatrixFromScilab(2 , &ARows_ , &ACols_ , &A_))
	{
		return 1;
	}
  int numCons_ = ARows_;           //Assigning number of inequality constraints = number of rows in Inequality constraint matrix
  int bRows_;

  //Inequality vector from Scilab
  if(getDoubleMatrixFromScilab(3 , &bRows_ , &ARows_ , &b_))
	{
		return 1;
		
	}

  //Equality constraint matrix from Scilab 
  int AeqRows_;                                       //Number of rows in Equality constraint matrix 
  int AeqCols_;                                       //Number of columns in Equality constraint matrix
  int numEqCons_;                                     //Indicates number of equality constraints
  if(getDoubleMatrixFromScilab(4 , &AeqRows_ , &AeqCols_ , &Aeq_))
	{
		return 1;
	}
  numEqCons_ = AeqRows_;               //Assigning number of equality constraints = number of rows in Equality constraint matrix

  //Equality vector from Scilab
  int beqRows_;                                       //Number of rows in Equality vector
  int beqCols_;                                       //specifying number of columns in Equality vector as it can be null matrix
  if(getDoubleMatrixFromScilab(5 , &beqRows_ , &beqCols_ , &beq_))
	{
		return 1;	
	}

  //Lower bounds for variables from Scilab
  int lbCols_;                                        //Number of columns in lower bound matrix
  if(getFixedSizeDoubleMatrixFromScilab(6 , 1 , numVars_ , &conLB_))
	{
		return 1;	
	}

  //Upper bounds for variables from Scilab
  int ubCols_;                                        //Number of columns in upper bound matrix
  if(getFixedSizeDoubleMatrixFromScilab(7 , 1 , numVars_ , &conUB_))
	{
		return 1;	
	}

  //get options from scilab
  if(getFixedSizeDoubleMatrixInList(8 , 2 , 1 , 1 , &options_))
  	{
		return 1;      
	}

  //get mps flag value
  if(getDoubleFromScilab(9 , &flagMps_))
  	{
		return 1;	
	}

  if(flagMps_ == 1)
  	{
    	int *piAddressVarNine = NULL;                //pointer used to access argument of the function
    	char mpsPath_[100];
  		mpsFile_ = mpsPath_;

        //load address of 1st argument into piAddressVarOne
	    sciErr = getVarAddressFromPosition(pvApiCtx, 11, &piAddressVarNine);

	    //check whether there is an error or not.
	    if (sciErr.iErr)
        	{
        		printError(&sciErr, 0);
        		return 1;
	        }
	    if ( !isStringType(pvApiCtx,piAddressVarNine) )
    		{
				Scierror(999,"Wrong type for input argument 1: A file name is expected.\n");
				return 1;
			}

        //read the value in that pointer pointing to file name
	    int err=getAllocatedSingleString(pvApiCtx, piAddressVarNine, &mpsFile_);
	}
  else 
	{
		mpsFile_ = "";
	}

   //get lp flag value
  if(getDoubleFromScilab(10 , &flagLp_))
  	{
		return 1;	
	}

  if(flagLp_ == 1)
  	{
    	int *piAddressVarTen = NULL;                //pointer used to access argument of the function
    	char lpPath_[100];
  		lpFile_ = lpPath_;

        //load address of 1st argument into piAddressVarOne
	    sciErr = getVarAddressFromPosition(pvApiCtx, 12, &piAddressVarTen);

	    //check whether there is an error or not.
	    if (sciErr.iErr)
        	{
        		printError(&sciErr, 0);
        		return 1;
	        }
	    if ( !isStringType(pvApiCtx,piAddressVarTen) )
    		{
				Scierror(999,"Wrong type for input argument 1: A file name is expected.\n");
				return 1;
			}

        //read the value in that pointer pointing to file name
	    int err=getAllocatedSingleString(pvApiCtx, piAddressVarTen, &lpFile_);
	}
  else 
	{
		lpFile_ = "";
	}
   
   //Call to the Clp Solver
   LinCLP* Prob = new LinCLP(numVars_,numCons_,numEqCons_,obj_,A_,b_,Aeq_,beq_,conLB_,conUB_,options_,flagMps_,flagLp_,mpsFile_,lpFile_);

   //Output the solution to Scilab
   //get solution for x
   double* xValue = Prob->getX();
   
   //get objective value
   double objValue = Prob->getObjVal();

   //get Status value
   double status = Prob->returnStatus();

   //get maximum constraint violation
   double maxViolation = Prob->consViolation();

   //get number of iterations
   double iterations = Prob->iterCount();

   //get reduced cost
   double* reducedCost = Prob->getReducedCost();

   //get dual vector
   double* ineqlin = Prob->getInEqlin();
   
   //get dual vector
   double* eqlin = Prob->getEqlin();

   returnDoubleMatrixToScilab(1 , 1 , numVars_ , xValue);
   returnDoubleMatrixToScilab(2 , 1 , 1 , &objValue);
   returnDoubleMatrixToScilab(3 , 1 , 1 , &status);
   returnDoubleMatrixToScilab(4 , 1 , 1 , &maxViolation);
   returnDoubleMatrixToScilab(5 , 1 , 1 , &iterations);
   returnDoubleMatrixToScilab(6 , 1 , numVars_ , reducedCost);
   returnDoubleMatrixToScilab(7 , 1 , numCons_ , ineqlin);
   returnDoubleMatrixToScilab(8 , 1 , numEqCons_ , eqlin);

}
}


  	
   
