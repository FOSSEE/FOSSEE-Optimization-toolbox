
// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Keyur Joshi
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "symphony.h"
#include "sci_iofunc.hpp"

extern sym_environment* global_sym_env; //defined in globals.cpp

extern "C" {
#include "api_scilab.h"
#include "Scierror.h"
#include "sci_malloc.h"
#include <localization.h>
#include "BOOL.h"
#include <string.h>

//error management variables
static SciErr sciErr;
static int iRet;

//data declarations
static int *varAddress=NULL,numVars,numConstr,*conMatrixColStart=NULL,*conMatrixRowIndex=NULL,*isIntVarBool=NULL,colIter,rowIter,inputMatrixCols,inputMatrixRows;
static double inputDouble,objSense,*objective=NULL,*lowerBounds=NULL,*upperBounds=NULL,*conLower=NULL,*conUpper=NULL,*conRange=NULL,*conRHS=NULL,*conMatrix=NULL;
static char *conType=NULL,*isIntVar=NULL;

//delete all allocd arrays before exit, and return output argument
static void cleanupBeforeExit()
	{
	if(conMatrixColStart) delete[] conMatrixColStart;
	if(conMatrixRowIndex) delete[] conMatrixRowIndex;
	if(isIntVar) delete[] isIntVar;
	if(conType) delete[] conType;
	if(conRange) delete[] conRange;
	if(conRHS) delete[] conRHS;
	iRet = createScalarDouble(pvApiCtx, nbInputArgument(pvApiCtx)+1,0);
	if(iRet)
	{
		/* If error, no return variable */
		AssignOutputVariable(pvApiCtx, 1) = 0;
		return;
	}
	AssignOutputVariable(pvApiCtx, 1) = nbInputArgument(pvApiCtx)+1;
	//ReturnArguments(pvApiCtx);
}

static int checkNumArgs()
{
	CheckInputArgument(pvApiCtx,10,10);
	CheckOutputArgument(pvApiCtx,1,1);
	return 1;
}

//both basic and advanced loader use this code
static int commonCodePart1()
{
	
	//ensure that environment is active
	if(global_sym_env==NULL)
	{
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		return 1;
	}
	
	//code to check arguments and get them
	if(checkNumArgs()==0)
		return 1;
	
	//get input 1: number of variables
	if(getUIntFromScilab(1,&numVars))
		return 1;
	
	//get input 2: number of constraints
	if(getUIntFromScilab(2,&numConstr))
		return 1;
	
	//allocate and prepare some arrays
	isIntVar=new char[numVars]; //is the variable constrained to be an integer?
	conType=new char[numConstr]; //char representing constraint type
	conRange=new double[numConstr]; //range of each constraint
	conRHS=new double[numConstr]; //RHS to be given to Symphony

	return 0;
}

//both basic and advanced loader use this code
static int commonCodePart2()
{
	//get input 3: lower bounds of variables
	if(getFixedSizeDoubleMatrixFromScilab(3,1,numVars,&lowerBounds))
	{
		cleanupBeforeExit();
		return 1;
	}
	
	//get input 4: upper bounds of variables
	if(getFixedSizeDoubleMatrixFromScilab(4,1,numVars,&upperBounds))
	{
		cleanupBeforeExit();
		return 1;
	}
	
	//get input 5: coefficients of variables in objective function to be minimized
	if(getFixedSizeDoubleMatrixFromScilab(5,1,numVars,&objective))
	{
		cleanupBeforeExit();
		return 1;
	}
	
	//get input 6: array that specifies wether a variable is constrained to be an integer
	sciErr = getVarAddressFromPosition(pvApiCtx, 6, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		cleanupBeforeExit();return 1;
	}
	if ( !isBooleanType(pvApiCtx, varAddress) )
	{
		Scierror(999, "Wrong type for input argument #6: A matrix of booleans is expected.\n");
		cleanupBeforeExit();return 1;
	}
	sciErr = getMatrixOfBoolean(pvApiCtx, varAddress, &inputMatrixRows, &inputMatrixCols, &isIntVarBool);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		cleanupBeforeExit();return 1;
	}
	if(inputMatrixRows!=1 || inputMatrixCols!=numVars)
	{
		Scierror(999, "Wrong type for input argument #6: Incorrectly sized matrix.\n");
		cleanupBeforeExit();return 1;
	}
	for(colIter=0;colIter<numVars;colIter++)
	{
		if(isIntVarBool[colIter])
			isIntVar[colIter]=TRUE;
		else
			isIntVar[colIter]=FALSE;
	}
	
	//get input 7: wether to minimize or maximize objective
	sciErr = getVarAddressFromPosition(pvApiCtx, 7, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		return 1;
	}
	if ( !isDoubleType(pvApiCtx,varAddress) ||  isVarComplex(pvApiCtx,varAddress) )
	{
		Scierror(999, "Wrong type for input argument #7: Either 1 (sym_minimize) or -1 (sym_maximize) is expected.\n");
		return 1;
	}
	iRet = getScalarDouble(pvApiCtx, varAddress, &objSense);
	if(iRet || (objSense!=-1 && objSense!=1))
	{
		Scierror(999, "Wrong type for input argument #7: Either 1 (sym_minimize) or -1 (sym_maximize) is expected.\n");
		return 1;
	}
	iRet=sym_set_obj_sense(global_sym_env,objSense);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY)
	{
		Scierror(999, "An error occured.\n");
		return 1;
	}
	
	//get input 9: constraint lower bound
	if(!(numConstr == 0))
	{	
		if(getFixedSizeDoubleMatrixFromScilab(9,numConstr,1,&conLower) )
		{
			cleanupBeforeExit();
			return 1;
		}
	}	
	//get input 10: constraint upper bound
	if(!(numConstr == 0))
	{
		if(getFixedSizeDoubleMatrixFromScilab(10,numConstr,1,&conUpper) && (numConstr == 0) )
		{
			cleanupBeforeExit();
			return 1;
		}
	}	
	//deduce type of constraint
	for(rowIter=0;rowIter<numConstr;rowIter++)
	{
		if(conLower[rowIter]>conUpper[rowIter])
		{
			Scierror(999, "Error: the lower bound of constraint %d is more than its upper bound.\n",rowIter);
			cleanupBeforeExit();
			return 1;
		}
		//#ifdef _MSC_VER
  	      //  double INFINITY = sym_get_infinity();
		//#endif
		if(conLower[rowIter]==(-INFINITY) && conUpper[rowIter]==INFINITY){
			conType[rowIter]='N';
			conRange[rowIter]=0;
			conRHS[rowIter]=0;
		}else if(conLower[rowIter]==(-INFINITY)){
			conType[rowIter]='L';
			conRange[rowIter]=0;
			conRHS[rowIter]=conUpper[rowIter];
		}else if(conUpper[rowIter]==INFINITY){
			conType[rowIter]='G';
			conRange[rowIter]=0;
			conRHS[rowIter]=conLower[rowIter];
		}else if(conUpper[rowIter]==conLower[rowIter]){
			conType[rowIter]='E';
			conRange[rowIter]=0;
			conRHS[rowIter]=conLower[rowIter];
		}else{
			conType[rowIter]='R';
			conRange[rowIter]=conUpper[rowIter]-conLower[rowIter];
			conRHS[rowIter]=conUpper[rowIter];
		}
	}
	
	/*
	//for debug: show all data
	sciprint("Vars: %d Constr: %d ObjType: %lf\n",numVars,numConstr,objSense);
	for(colIter=0;colIter<numVars;colIter++)
		sciprint("Var %d: upper: %lf lower: %lf isInt: %d ObjCoeff: %lf\n",colIter,lowerBounds[colIter],upperBounds[colIter],isIntVar[colIter],objective[colIter]);
	for(rowIter=0;rowIter<numConstr;rowIter++)
		sciprint("Constr %d: type: %c lower: %lf upper: %lf range: %lf\n",rowIter,conType[rowIter],conLower[rowIter],conRange[rowIter]);
	*/
	
	//call problem loader
	sym_explicit_load_problem(global_sym_env,numVars,numConstr,conMatrixColStart,conMatrixRowIndex,conMatrix,lowerBounds,upperBounds,isIntVar,objective,NULL,conType,conRHS,conRange,TRUE);
	// sciprint("Problem loaded into environment.\n");
	
	//code to give output
	cleanupBeforeExit();
	
	return 0;
}


//advanced problem loader, expects sparse matrix. For use with larger problems (>10 vars)
int sci_sym_loadProblem(char *fname)
{
	int retVal,nonZeros,*itemsPerRow,*colIndex,matrixIter,newPos,*oldRowIndex,*colStartCopy;
	double *data;
	
	if(commonCodePart1())
		return 1;
		
	//get input 8: matrix of constraint equation coefficients
	sciErr = getVarAddressFromPosition(pvApiCtx, 8, &varAddress);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		cleanupBeforeExit();return 1;
	}
	
	if(!(numConstr == 0))
	{

	if ( (!isSparseType(pvApiCtx,varAddress) ||  isVarComplex(pvApiCtx,varAddress)) && (numConstr == 0))
	{
		Scierror(999, "Wrong type for input argument #8: A sparse matrix of doubles is expected.\n");
		cleanupBeforeExit();return 1;
	}
	}
	sciErr = getSparseMatrix(pvApiCtx,varAddress,&inputMatrixRows,&inputMatrixCols,&nonZeros,&itemsPerRow,&colIndex,&data);
	if (sciErr.iErr)
	{
		printError(&sciErr, 0);
		cleanupBeforeExit();return 1;
	}
	if(inputMatrixRows!=numConstr || inputMatrixCols!=numVars)
	{
		Scierror(999, "Wrong type for input argument #8: Incorrectly sized matrix.\n");
		cleanupBeforeExit();return 1;
	}
	
	//convert SciLab format sparse matrix into the format required by Symphony
	conMatrix=new double[nonZeros]; //matrix contents
	conMatrixColStart=new int[numVars+1]; //where each column of the matrix starts
	conMatrixRowIndex=new int[nonZeros]; //row number of each element
	oldRowIndex=new int[nonZeros]; //row number in old matrix
	colStartCopy=new int[numVars+1]; //temporary copy of conMatrixColStart
	for(rowIter=matrixIter=0;rowIter<numConstr;rowIter++) //assign row number to each element in old matrix
		for(colIter=0;colIter<itemsPerRow[rowIter];colIter++,matrixIter++)
			oldRowIndex[matrixIter]=rowIter;
	for(colIter=0;colIter<=numVars;colIter++) //initialize to 0
		conMatrixColStart[colIter]=0;
	for(matrixIter=0;matrixIter<nonZeros;matrixIter++) //get number of elements in each column
		conMatrixColStart[colIndex[matrixIter]]++;
	for(colIter=1;colIter<=numVars;colIter++) //perfrom cumulative addition to get final data about where each column starts
	{
		conMatrixColStart[colIter]+=conMatrixColStart[colIter-1];
		colStartCopy[colIter]=conMatrixColStart[colIter];
	}
	colStartCopy[0]=0;
	for(matrixIter=0;matrixIter<nonZeros;matrixIter++) //move data from old matrix to new matrix
	{
		newPos=colStartCopy[colIndex[matrixIter]-1]++; //calculate its position in the new matrix
		conMatrix[newPos]=data[matrixIter]; //move the data
		conMatrixRowIndex[newPos]=oldRowIndex[matrixIter]; //assign it its row number
	}
	
	retVal=commonCodePart2();
	
	//cleanup some more allocd memory
	if(conMatrix) delete[] conMatrix;
	if(oldRowIndex) delete[] oldRowIndex;
	if(colStartCopy) delete[] colStartCopy;
	
	return retVal;
}

}
