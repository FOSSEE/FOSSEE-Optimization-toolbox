// Symphony Toolbox for Scilab
// (Definition of) Functions for input and output from Scilab
// By Keyur Joshi

#include "api_scilab.h"
#include "Scierror.h"
#include "sciprint.h"
#include "BOOL.h"
#include <localization.h>
#include "call_scilab.h"
#include <string.h>


using namespace std;
int getFunctionFromScilab(int argNum, int **dest)
{	
	//data declarations
	SciErr sciErr;
	int iRet,*varAddress, iType;
	double inputDouble;
	const char errMsg[]="Wrong type for input argument #%d: A function is expected.\n";
	const int errNum=999;
	//get variable address
	sciErr = getVarAddressFromPosition(pvApiCtx, 1, dest);
	if(sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	//check that the variable is necessarily a function
	sciErr = getVarType(pvApiCtx, *dest, &iType);
  	if(sciErr.iErr || iType != sci_c_function)
	{
		Scierror(errNum,errMsg,argNum);
		return 1;
	}
	return 0;
	
}

int getDoubleFromScilab(int argNum, double *dest)
{
	//data declarations
	SciErr sciErr;
	int iRet,*varAddress;
	const char errMsg[]="Wrong type for input argument #%d: A double is expected.\n";
	const int errNum=999;
	//get variable address
	sciErr = getVarAddressFromPosition(pvApiCtx, argNum, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	//check that it is a non-complex double
	if ( !isDoubleType(pvApiCtx,varAddress) ||  isVarComplex(pvApiCtx,varAddress) )
	{
		Scierror(errNum,errMsg,argNum);
		return 1;
	}
	//retrieve and store
	iRet = getScalarDouble(pvApiCtx, varAddress, dest);
	if(iRet)
	{
		Scierror(errNum,errMsg,argNum);
		return 1;
	}
	return 0;
}

int getUIntFromScilab(int argNum, int *dest)
{
	SciErr sciErr;
	int iRet,*varAddress;
	double inputDouble;
	const char errMsg[]="Wrong type for input argument #%d: A nonnegative integer is expected.\n";
	const int errNum=999;
	//same steps as above
	sciErr = getVarAddressFromPosition(pvApiCtx, argNum, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	if ( !isDoubleType(pvApiCtx,varAddress) ||  isVarComplex(pvApiCtx,varAddress) )
	{
		Scierror(errNum,errMsg,argNum);
		return 1;
	}
	iRet = getScalarDouble(pvApiCtx, varAddress, &inputDouble);
	//check that an unsigned int is stored in the double by casting and recasting
	if(iRet || ((double)((unsigned int)inputDouble))!=inputDouble)
	{
		Scierror(errNum,errMsg,argNum);
		return 1;
	}
	*dest=(unsigned int)inputDouble;
	return 0;
}

int getIntFromScilab(int argNum, int *dest)
{
	SciErr sciErr;
	int iRet,*varAddress;
	double inputDouble;
	const char errMsg[]="Wrong type for input argument #%d: An integer is expected.\n";
	const int errNum=999;
	//same steps as above
	sciErr = getVarAddressFromPosition(pvApiCtx, argNum, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	if ( !isDoubleType(pvApiCtx,varAddress) ||  isVarComplex(pvApiCtx,varAddress) )
	{
		Scierror(errNum,errMsg,argNum);
		return 1;
	}
	iRet = getScalarDouble(pvApiCtx, varAddress, &inputDouble);
	//check that an int is stored in the double by casting and recasting
	if(iRet || ((double)((int)inputDouble))!=inputDouble)
	{
		Scierror(errNum,errMsg,argNum);
		return 1;
	}
	*dest=(int)inputDouble;
	return 0;
}

int getFixedSizeDoubleMatrixFromScilab(int argNum, int rows, int cols, double **dest)
{
	int *varAddress,inputMatrixRows,inputMatrixCols;
	SciErr sciErr;
	const char errMsg[]="Wrong type for input argument #%d: A matrix of double of size %d by %d is expected.\n";
	const int errNum=999;
	//same steps as above
	sciErr = getVarAddressFromPosition(pvApiCtx, argNum, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	if ( !isDoubleType(pvApiCtx,varAddress) ||  isVarComplex(pvApiCtx,varAddress) )
	{
		Scierror(errNum,errMsg,argNum,rows,cols);
		return 1;
	}
	sciErr = getMatrixOfDouble(pvApiCtx, varAddress, &inputMatrixRows, &inputMatrixCols,NULL);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	//check that the matrix has the correct number of rows and columns
	if(inputMatrixRows!=rows || inputMatrixCols!=cols)
	{
		Scierror(errNum,errMsg,argNum,rows,cols);
		return 1;
	}
	getMatrixOfDouble(pvApiCtx, varAddress, &inputMatrixRows, &inputMatrixCols, dest);
	return 0;
}

int getDoubleMatrixFromScilab(int argNum, int *rows, int *cols, double **dest)
{
	int *varAddress;
	SciErr sciErr;
	const char errMsg[]="Wrong type for input argument #%d: A matrix of double is expected.\n";
	const int errNum=999;
	//same steps as above
	sciErr = getVarAddressFromPosition(pvApiCtx, argNum, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	if ( !isDoubleType(pvApiCtx,varAddress) ||  isVarComplex(pvApiCtx,varAddress) )
	{
		Scierror(errNum,errMsg,argNum);
		return 1;
	}
	getMatrixOfDouble(pvApiCtx, varAddress, rows, cols, dest);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	return 0;
}

int getFixedSizeDoubleMatrixInList(int argNum, int itemPos, int rows, int cols, double **dest)
{
	int *varAddress,inputMatrixRows,inputMatrixCols;
	SciErr sciErr;
	const char errMsg[]="Wrong type for input argument #%d: A matrix of double of size %d by %d is expected.\n";
	const int errNum=999;
	//same steps as above
	sciErr = getVarAddressFromPosition(pvApiCtx, argNum, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}

	getMatrixOfDoubleInList(pvApiCtx, varAddress, itemPos, &rows, &cols, dest);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	return 0;
}

int getStringFromScilab(int argNum,char **dest)
{
	int *varAddress,inputMatrixRows,inputMatrixCols;
	SciErr sciErr;
	sciErr = getVarAddressFromPosition(pvApiCtx, argNum, &varAddress);

	//check whether there is an error or not.
	if (sciErr.iErr)
    	{
        	printError(&sciErr, 0);
        	return 1;
		}
	if ( !isStringType(pvApiCtx,varAddress) )
		{
			Scierror(999,"Wrong type for input argument 1: A file name is expected.\n");
			return 1;
		}
    //read the value in that pointer pointing to file name
	getAllocatedSingleString(pvApiCtx, varAddress, dest);
    
}

bool getFunctionFromScilab1(int n,char name[], double *x,int posFirstElementOnStackForSF,int nOfRhsOnSF,int nOfLhsOnSF, double **dest)
{	
	double check;
  	createMatrixOfDouble(pvApiCtx, posFirstElementOnStackForSF, 1, n, x);  
  	C2F(scistring)(&posFirstElementOnStackForSF,name,&nOfLhsOnSF,&nOfRhsOnSF,(unsigned long)strlen(name));
    
  	if(getDoubleFromScilab(posFirstElementOnStackForSF+1,&check))
  	{
		return true;
	}
	if (check==1)
	{
		return true;
	}	
	else
	{ 
		int x_rows, x_cols;
		if(getDoubleMatrixFromScilab(posFirstElementOnStackForSF, &x_rows, &x_cols, dest))
  		{
			sciprint("No results ");
			return true;
			
  		}
	}	
	return 0;
}

bool getHessFromScilab(int n,int numConstr_,char name[], double *x,double *obj,double *lambda,int posFirstElementOnStackForSF,int nOfRhsOnSF,int nOfLhsOnSF, double **dest)
{	
	double check;
  	createMatrixOfDouble(pvApiCtx, posFirstElementOnStackForSF, 1, n, x);
  	createMatrixOfDouble(pvApiCtx, posFirstElementOnStackForSF+1, 1, 1, obj);
	createMatrixOfDouble(pvApiCtx, posFirstElementOnStackForSF+2, 1, numConstr_, lambda);
  	C2F(scistring)(&posFirstElementOnStackForSF,name,&nOfLhsOnSF,&nOfRhsOnSF,(unsigned long)strlen(name));
                               
  	if(getDoubleFromScilab(posFirstElementOnStackForSF+1,&check))
  	{
		return true;
	}
	if (check==1)
	{
		return true;
	}	
	else
	{ 
		int x_rows, x_cols;
		if(getDoubleMatrixFromScilab(posFirstElementOnStackForSF, &x_rows, &x_cols, dest))
  		{
			sciprint("No results ");
			return 1;	
  		}
	}	
	return 0;
}

int getIntMatrixFromScilab(int argNum, int *rows, int *cols, int **dest)
{
	int *varAddress;
	SciErr sciErr;
	const char errMsg[]="Wrong type for input argument #%d: A matrix of integer is expected.\n";
	const int errNum=999;
	//same steps as above
	sciErr = getVarAddressFromPosition(pvApiCtx, argNum, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	// if ( !isIntegerType(pvApiCtx,varAddress) ||  isVarComplex(pvApiCtx,varAddress) )
	// {
	// 	Scierror(errNum,errMsg,argNum);
	// 	return 1;
	// }
	getMatrixOfInteger32(pvApiCtx, varAddress, rows, cols, dest);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	return 0;
}

int return0toScilab()
{
	int iRet;
	//create variable in scilab
	iRet = createScalarDouble(pvApiCtx, nbInputArgument(pvApiCtx)+1,0);
	if(iRet)
	{
		/* If error, no return variable */
		AssignOutputVariable(pvApiCtx, 1) = 0;
		return 1;
	}
	//make it the output variable
	AssignOutputVariable(pvApiCtx, 1) = nbInputArgument(pvApiCtx)+1;
	//return it to scilab
	//ReturnArguments(pvApiCtx);
	return 0;
}

int returnDoubleToScilab(double retVal)
{
	int iRet;
	//same steps as above
	iRet = createScalarDouble(pvApiCtx, nbInputArgument(pvApiCtx)+1,retVal);
	if(iRet)
	{
		/* If error, no return variable */
		AssignOutputVariable(pvApiCtx, 1) = 0;
		return 1;
	}
	AssignOutputVariable(pvApiCtx, 1) = nbInputArgument(pvApiCtx)+1;
	//ReturnArguments(pvApiCtx);
	return 0;
}

int returnDoubleMatrixToScilab(int itemPos, int rows, int cols, const double *dest)
{
	SciErr sciErr;
	//same steps as above
	sciErr = createMatrixOfDouble(pvApiCtx, nbInputArgument(pvApiCtx) + itemPos, rows, cols, dest);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}

	AssignOutputVariable(pvApiCtx, itemPos) = nbInputArgument(pvApiCtx)+itemPos;

	return 0;
}

int returnIntegerMatrixToScilab(int itemPos, int rows, int cols, const int *dest)
{
	SciErr sciErr;
	//same steps as above
	sciErr = createMatrixOfInteger32(pvApiCtx, nbInputArgument(pvApiCtx) + itemPos, rows, cols, dest);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}

	AssignOutputVariable(pvApiCtx, itemPos) = nbInputArgument(pvApiCtx)+itemPos;

	return 0;
}

