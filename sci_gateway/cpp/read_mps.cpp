/*
 * Linear Solver Toolbox for Scilab using CLP library
 * Authors :
	Guru Pradeep Reddy
        Bhanu Priya Sayal
*/

#include "sci_iofunc.hpp"
#include "OsiClpSolverInterface.hpp"

extern "C"{
#include <api_scilab.h>
#include <Scierror.h>
#include <localization.h>
#include <sciprint.h>

//Solver function
int sci_rmps(char *fname) 
{
    //creating a problem pointer using base class of OsiSolverInterface and
    //instantiate the object using derived class of ClpSolverInterface
    OsiSolverInterface* si = new OsiClpSolverInterface();

    // Error management variable
	SciErr sciErr;

	//data declarations
	int *piAddressVarOne = NULL;                 //pointer used to access argument of the function
	char file[100];                              //string to hold the name of .mps file
	char* ptr=file;                              //pointer to point to address of file name
    double* options_;                            //options to set maximum iterations 
	CheckInputArgument(pvApiCtx, 2,2 );          //Check we have exactly two arguments as input or not
	CheckOutputArgument(pvApiCtx, 6, 6);         //Check we have exactly six arguments on output side or not

    //Getting the input arguments from Scilab
    //Getting the MPS file path
	//load address of 1st argument into piAddressVarOne
	sciErr = getVarAddressFromPosition(pvApiCtx, 1, &piAddressVarOne);

	//check whether there is an error or not.
	if (sciErr.iErr)
    	{
        	printError(&sciErr, 0);
        	return 1;
		}
	if ( !isStringType(pvApiCtx,piAddressVarOne) )
		{
			Scierror(999,"Wrong type for input argument 1: A file name is expected.\n");
			return 1;
		}
    //read the value in that pointer pointing to file name
	int err=getAllocatedSingleString(pvApiCtx, piAddressVarOne, &ptr);
    
    //get options from Scilab
    if(getFixedSizeDoubleMatrixInList(2 , 2 , 1 , 1 , &options_))
  		{
			return 1;      
		}

    //Read the MPS file
    si->readMps(ptr);

    //setting options for maximum iterations
    si->setIntParam(OsiMaxNumIteration,options_[0]);

    //Solve the problem
    si->initialSolve();
  
    //Quering about the problem
    //get number of variables
    double numVars_;
    numVars_ = si->getNumCols();
  
    //get number of constraint equations
    double numCons_;
    numCons_ = si->getNumRows();
   
    //Output the solution to Scilab
    //get solution for x
    double* xValue = si->getColSolution();
   
    //get objective value
    double objValue = si->getObjValue();

    //get Status value
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

    //get number of iterations
    double iterations = si->getIterationCount();

    //get reduced cost 
    double* reducedCost = si->getReducedCost();
   
    //get dual vector
    double* dual = si->getRowPrice();
  
    returnDoubleMatrixToScilab(1 , 1 , numVars_ , xValue);
    returnDoubleMatrixToScilab(2 , 1 , 1 , &objValue);
    returnDoubleMatrixToScilab(3 , 1 , 1 , &status);
    returnDoubleMatrixToScilab(4 , 1 , 1 , &iterations);
    returnDoubleMatrixToScilab(5 , 1 , numVars_ , reducedCost);
    returnDoubleMatrixToScilab(6 , 1 , numCons_ , dual);

}
}


  	
   
