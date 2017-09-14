// Symphony Toolbox for Scilab
// (Declaration of) Functions for input and output from Scilab
// By Keyur Joshi

#ifndef SCI_IOFUNCHEADER
#define SCI_IOFUNCHEADER

//input
int getFunctionFromScilab(int argNum, int **dest);
int getDoubleFromScilab(int argNum, double *dest);
int getUIntFromScilab(int argNum, int *dest);
int getIntFromScilab(int argNum, int *dest);
int getFixedSizeDoubleMatrixFromScilab(int argNum, int rows, int cols, double **dest);
int getDoubleMatrixFromScilab(int argNum, int *rows, int *cols, double **dest);
int getFixedSizeDoubleMatrixInList(int argNum, int itemPos, int rows, int cols, double **dest);
int getStringFromScilab(int argNum,char** dest);
bool getFunctionFromScilab1(int n,char name[], double *x,int posFirstElementOnStackForSF,int nOfRhsOnSF,int nOfLhsOnSF, double **dest);
bool getHessFromScilab(int n,int numConstr_,char name[], double *x,double *obj,double *lambda,int posFirstElementOnStackForSF,int nOfRhsOnSF,int nOfLhsOnSF, double **dest);
int getIntMatrixFromScilab(int argNum, int *rows, int *cols, int **dest);

//output
int return0toScilab();
int returnDoubleToScilab(double retVal);
int returnDoubleMatrixToScilab(int itemPos, int rows, int cols, const double *dest);
int returnIntegerMatrixToScilab(int itemPos, int rows, int cols, const int *dest);

#endif //SCI_IOFUNCHEADER