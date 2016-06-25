#ifdef __cplusplus
extern "C" {
#endif
#include <mex.h> 
#include <sci_gateway.h>
#include <api_scilab.h>
#include <MALLOC.h>
static int direct_gateway(char *fname,void F(void)) { F();return 0;};
extern Gatefunc sci_sym_open;
extern Gatefunc sci_sym_close;
extern Gatefunc sci_sym_isEnvActive;
extern Gatefunc sci_sym_set_defaults;
extern Gatefunc sci_sym_set_int_param;
extern Gatefunc sci_sym_get_int_param;
extern Gatefunc sci_sym_set_dbl_param;
extern Gatefunc sci_sym_get_dbl_param;
extern Gatefunc sci_sym_set_str_param;
extern Gatefunc sci_sym_get_str_param;
extern Gatefunc sci_sym_getInfinity;
extern Gatefunc sci_sym_loadProblemBasic;
extern Gatefunc sci_sym_loadProblem;
extern Gatefunc sci_sym_load_mps;
extern Gatefunc sci_sym_get_num_int;
extern Gatefunc sci_sym_get_num_int;
extern Gatefunc sci_sym_get_num_int;
extern Gatefunc sci_sym_isContinuous;
extern Gatefunc sci_sym_isBinary;
extern Gatefunc sci_sym_isInteger;
extern Gatefunc sci_sym_set_continuous;
extern Gatefunc sci_sym_set_integer;
extern Gatefunc sci_sym_get_dbl_arr;
extern Gatefunc sci_sym_get_dbl_arr;
extern Gatefunc sci_sym_setVarBound;
extern Gatefunc sci_sym_setVarBound;
extern Gatefunc sci_sym_get_dbl_arr;
extern Gatefunc sci_sym_setObjCoeff;
extern Gatefunc sci_sym_getObjSense;
extern Gatefunc sci_sym_setObjSense;
extern Gatefunc sci_sym_get_dbl_arr;
extern Gatefunc sci_sym_get_dbl_arr;
extern Gatefunc sci_sym_get_dbl_arr;
extern Gatefunc sci_sym_get_dbl_arr;
extern Gatefunc sci_sym_setConstrBound;
extern Gatefunc sci_sym_setConstrBound;
extern Gatefunc sci_sym_setConstrType;
extern Gatefunc sci_sym_get_matrix;
extern Gatefunc sci_sym_get_row_sense;
extern Gatefunc sci_sym_addConstr;
extern Gatefunc sci_sym_addVar;
extern Gatefunc sci_sym_delete_cols;
extern Gatefunc sci_sym_delete_rows;
extern Gatefunc sci_sym_getPrimalBound;
extern Gatefunc sci_sym_setPrimalBound;
extern Gatefunc sci_sym_setColSoln;
extern Gatefunc sci_sym_solve;
extern Gatefunc sci_sym_get_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_getVarSoln;
extern Gatefunc sci_sym_getObjVal;
extern Gatefunc sci_sym_get_iteration_count;
extern Gatefunc sci_sym_getRowActivity;
extern Gatefunc sci_linearprog;
extern Gatefunc sci_rmps;
extern Gatefunc sci_solveqp;
extern Gatefunc sci_solveminuncp;
extern Gatefunc sci_solveminbndp;
extern Gatefunc sci_solveminconp;
extern Gatefunc cpp_intqpipopt;
static GenericTable Tab[]={
  {(Myinterfun)sci_gateway,sci_sym_open,"sym_open"},
  {(Myinterfun)sci_gateway,sci_sym_close,"sym_close"},
  {(Myinterfun)sci_gateway,sci_sym_isEnvActive,"sym_isEnvActive"},
  {(Myinterfun)sci_gateway,sci_sym_set_defaults,"sym_resetParams"},
  {(Myinterfun)sci_gateway,sci_sym_set_int_param,"sym_setIntParam"},
  {(Myinterfun)sci_gateway,sci_sym_get_int_param,"sym_getIntParam"},
  {(Myinterfun)sci_gateway,sci_sym_set_dbl_param,"sym_setDblParam"},
  {(Myinterfun)sci_gateway,sci_sym_get_dbl_param,"sym_getDblParam"},
  {(Myinterfun)sci_gateway,sci_sym_set_str_param,"sym_setStrParam"},
  {(Myinterfun)sci_gateway,sci_sym_get_str_param,"sym_getStrParam"},
  {(Myinterfun)sci_gateway,sci_sym_getInfinity,"sym_getInfinity"},
  {(Myinterfun)sci_gateway,sci_sym_loadProblemBasic,"sym_loadProblemBasic"},
  {(Myinterfun)sci_gateway,sci_sym_loadProblem,"sym_loadProblem"},
  {(Myinterfun)sci_gateway,sci_sym_load_mps,"sym_loadMPS"},
  {(Myinterfun)sci_gateway,sci_sym_get_num_int,"sym_getNumConstr"},
  {(Myinterfun)sci_gateway,sci_sym_get_num_int,"sym_getNumVar"},
  {(Myinterfun)sci_gateway,sci_sym_get_num_int,"sym_getNumElements"},
  {(Myinterfun)sci_gateway,sci_sym_isContinuous,"sym_isContinuous"},
  {(Myinterfun)sci_gateway,sci_sym_isBinary,"sym_isBinary"},
  {(Myinterfun)sci_gateway,sci_sym_isInteger,"sym_isInteger"},
  {(Myinterfun)sci_gateway,sci_sym_set_continuous,"sym_setContinuous"},
  {(Myinterfun)sci_gateway,sci_sym_set_integer,"sym_setInteger"},
  {(Myinterfun)sci_gateway,sci_sym_get_dbl_arr,"sym_getVarLower"},
  {(Myinterfun)sci_gateway,sci_sym_get_dbl_arr,"sym_getVarUpper"},
  {(Myinterfun)sci_gateway,sci_sym_setVarBound,"sym_setVarLower"},
  {(Myinterfun)sci_gateway,sci_sym_setVarBound,"sym_setVarUpper"},
  {(Myinterfun)sci_gateway,sci_sym_get_dbl_arr,"sym_getObjCoeff"},
  {(Myinterfun)sci_gateway,sci_sym_setObjCoeff,"sym_setObjCoeff"},
  {(Myinterfun)sci_gateway,sci_sym_getObjSense,"sym_getObjSense"},
  {(Myinterfun)sci_gateway,sci_sym_setObjSense,"sym_setObjSense"},
  {(Myinterfun)sci_gateway,sci_sym_get_dbl_arr,"sym_getRhs"},
  {(Myinterfun)sci_gateway,sci_sym_get_dbl_arr,"sym_getConstrRange"},
  {(Myinterfun)sci_gateway,sci_sym_get_dbl_arr,"sym_getConstrLower"},
  {(Myinterfun)sci_gateway,sci_sym_get_dbl_arr,"sym_getConstrUpper"},
  {(Myinterfun)sci_gateway,sci_sym_setConstrBound,"sym_setConstrLower"},
  {(Myinterfun)sci_gateway,sci_sym_setConstrBound,"sym_setConstrUpper"},
  {(Myinterfun)sci_gateway,sci_sym_setConstrType,"sym_setConstrType"},
  {(Myinterfun)sci_gateway,sci_sym_get_matrix,"sym_getMatrix"},
  {(Myinterfun)sci_gateway,sci_sym_get_row_sense,"sym_getConstrSense"},
  {(Myinterfun)sci_gateway,sci_sym_addConstr,"sym_addConstr"},
  {(Myinterfun)sci_gateway,sci_sym_addVar,"sym_addVar"},
  {(Myinterfun)sci_gateway,sci_sym_delete_cols,"sym_deleteVars"},
  {(Myinterfun)sci_gateway,sci_sym_delete_rows,"sym_deleteConstrs"},
  {(Myinterfun)sci_gateway,sci_sym_getPrimalBound,"sym_getPrimalBound"},
  {(Myinterfun)sci_gateway,sci_sym_setPrimalBound,"sym_setPrimalBound"},
  {(Myinterfun)sci_gateway,sci_sym_setColSoln,"sym_setVarSoln"},
  {(Myinterfun)sci_gateway,sci_sym_solve,"sym_solve"},
  {(Myinterfun)sci_gateway,sci_sym_get_status,"sym_getStatus"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isOptimal"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isInfeasible"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isAbandoned"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isIterLimitReached"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isTimeLimitReached"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isTargetGapAchieved"},
  {(Myinterfun)sci_gateway,sci_sym_getVarSoln,"sym_getVarSoln"},
  {(Myinterfun)sci_gateway,sci_sym_getObjVal,"sym_getObjVal"},
  {(Myinterfun)sci_gateway,sci_sym_get_iteration_count,"sym_getIterCount"},
  {(Myinterfun)sci_gateway,sci_sym_getRowActivity,"sym_getConstrActivity"},
  {(Myinterfun)sci_gateway,sci_linearprog,"linearprog"},
  {(Myinterfun)sci_gateway,sci_rmps,"rmps"},
  {(Myinterfun)sci_gateway,sci_solveqp,"solveqp"},
  {(Myinterfun)sci_gateway,sci_solveminuncp,"solveminuncp"},
  {(Myinterfun)sci_gateway,sci_solveminbndp,"solveminbndp"},
  {(Myinterfun)sci_gateway,sci_solveminconp,"solveminconp"},
  {(Myinterfun)sci_gateway,cpp_intqpipopt,"sci_intqpipopt"},
};
 
int C2F(libFOSSEE_Optimization_Toolbox)()
{
  Rhs = Max(0, Rhs);
  if (*(Tab[Fin-1].f) != NULL) 
  {
     if(pvApiCtx == NULL)
     {
       pvApiCtx = (StrCtx*)MALLOC(sizeof(StrCtx));
     }
     pvApiCtx->pstName = (char*)Tab[Fin-1].name;
    (*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
  }
  return 0;
}
#ifdef __cplusplus
}
#endif
