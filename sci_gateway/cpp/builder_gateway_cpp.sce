// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

mode(-1)
lines(0)

toolbox_title = "FOSSEE_Optimization_Toolbox";

[a, opt] = getversion();
Version = opt(2);

Build_64Bits = %f;

path_builder = get_absolute_file_path('builder_gateway_cpp.sce');

tools_path  = path_builder + "../../thirdparty/linux/";

C_Flags=["-D__USE_DEPRECATED_STACK_FUNCTIONS__ -w -fpermissive -I"+tools_path+"include/coin -Wl,-rpath="+tools_path+"lib/"+Version+filesep()+" "]

Linker_Flag = ["-L"+tools_path+"lib/"+Version+filesep()+"libSym"+" "+"-L"+tools_path+"lib/"+Version+filesep()+"libipopt"+" "+"-L"+tools_path+"lib/"+Version+filesep()+"libClp"+" "+"-L"+tools_path+"lib/"+Version+filesep()+"libOsiClp"+" "+"-L"+tools_path+"lib/"+Version+filesep()+"libCoinUtils" ]


//Name of All the Functions
Function_Names = [
		//for opening/closing environment and checking if it is open/close
		"sym_open","sci_sym_open";
		"sym_close","sci_sym_close";
		"sym_isEnvActive","sci_sym_isEnvActive";
		
		//run time parameters
		"sym_resetParams","sci_sym_set_defaults";
		"sym_setIntParam","sci_sym_set_int_param";
		"sym_getIntParam","sci_sym_get_int_param";
		"sym_setDblParam","sci_sym_set_dbl_param";
		"sym_getDblParam","sci_sym_get_dbl_param";
		"sym_setStrParam","sci_sym_set_str_param";
		"sym_getStrParam","sci_sym_get_str_param";
		"sym_getInfinity","sci_sym_getInfinity";
		
		//problem loaders
		"sym_loadProblemBasic","sci_sym_loadProblemBasic";
		"sym_loadProblem","sci_sym_loadProblem";
		"sym_loadMPS","sci_sym_load_mps";
		
		//basic data
		"sym_getNumConstr","sci_sym_get_num_int";
		"sym_getNumVar","sci_sym_get_num_int";
		"sym_getNumElements","sci_sym_get_num_int";
		
		//variable and objective data
		"sym_isContinuous","sci_sym_isContinuous";
		"sym_isBinary","sci_sym_isBinary";
		"sym_isInteger","sci_sym_isInteger";
		"sym_setContinuous","sci_sym_set_continuous";
		"sym_setInteger","sci_sym_set_integer";
		"sym_getVarLower","sci_sym_get_dbl_arr";
		"sym_getVarUpper","sci_sym_get_dbl_arr";
		"sym_setVarLower","sci_sym_setVarBound";
		"sym_setVarUpper","sci_sym_setVarBound";
		"sym_getObjCoeff","sci_sym_get_dbl_arr";
		"sym_setObjCoeff","sci_sym_setObjCoeff";
		"sym_getObjSense","sci_sym_getObjSense";
		"sym_setObjSense","sci_sym_setObjSense";
		
		//constraint data
		"sym_getRhs","sci_sym_get_dbl_arr";
		"sym_getConstrRange","sci_sym_get_dbl_arr";
		"sym_getConstrLower","sci_sym_get_dbl_arr";
		"sym_getConstrUpper","sci_sym_get_dbl_arr";
		"sym_setConstrLower","sci_sym_setConstrBound";
		"sym_setConstrUpper","sci_sym_setConstrBound";
		"sym_setConstrType","sci_sym_setConstrType";
		"sym_getMatrix","sci_sym_get_matrix";
		"sym_getConstrSense","sci_sym_get_row_sense";
		
		//add/remove variables and constraints
		"sym_addConstr","sci_sym_addConstr";
		"sym_addVar","sci_sym_addVar";
		"sym_deleteVars","sci_sym_delete_cols";
		"sym_deleteConstrs","sci_sym_delete_rows";
		
		//primal bound
		"sym_getPrimalBound","sci_sym_getPrimalBound";
		"sym_setPrimalBound","sci_sym_setPrimalBound";
		
		//set preliminary solution
		"sym_setVarSoln","sci_sym_setColSoln";
		
		//solve
		"sym_solve","sci_sym_solve";
		
		//post solve functions
		"sym_getStatus","sci_sym_get_status";
		"sym_isOptimal","sci_sym_get_solver_status";
		"sym_isInfeasible","sci_sym_get_solver_status";
		"sym_isAbandoned","sci_sym_get_solver_status";
		"sym_isIterLimitReached","sci_sym_get_solver_status";
		"sym_isTimeLimitReached","sci_sym_get_solver_status";
		"sym_isTargetGapAchieved","sci_sym_get_solver_status";
		"sym_getVarSoln","sci_sym_getVarSoln";
		"sym_getObjVal","sci_sym_getObjVal";
		"sym_getIterCount","sci_sym_get_iteration_count";
		"sym_getConstrActivity","sci_sym_getRowActivity";

		//Linprog function
		"linearprog","sci_linearprog"
        "rmps","sci_rmps"

		//QP function
		"solveqp","sci_solveqp"

		//fminunc function and fminbnd function
		"solveminuncp","sci_solveminuncp"
		"solveminbndp","sci_solveminbndp"
		"solveminconp","sci_solveminconp"
		
		//intqpipopt functions
		"sci_intqpipopt", "cpp_intqpipopt"
	];

//Name of all the files to be compiled
Files = [
		"globals.cpp",
		"sci_iofunc.hpp",
		"sci_iofunc.cpp",
		"sci_sym_openclose.cpp",
		"sci_solver_status_query_functions.cpp",
		"sci_sym_solve.cpp",
		"sci_sym_loadproblem.cpp",
		"sci_sym_isenvactive.cpp",
		"sci_sym_load_mps.cpp",
		"sci_vartype.cpp",
		"sci_sym_getinfinity.cpp",
		"sci_sym_solution.cpp",
		"sym_data_query_functions.cpp"
		"sci_sym_set_variables.cpp",
		"sci_sym_setobj.cpp",
		"sci_sym_varbounds.cpp",
		"sci_sym_rowmod.cpp",
		"sci_sym_set_indices.cpp",
		"sci_sym_addrowcol.cpp",
		"sci_sym_primalbound.cpp",
		"sci_sym_setcolsoln.cpp",
		"sci_sym_getrowact.cpp",
		"sci_sym_getobjsense.cpp",
		"sci_sym_remove.cpp",
		"sci_QuadNLP.cpp",
		"QuadNLP.hpp",
		"sci_ipopt.cpp",
		"minuncNLP.hpp",
		"sci_minuncNLP.cpp",
		"sci_ipoptfminunc.cpp",
		"minbndNLP.hpp",
		"sci_minbndNLP.cpp",
		"sci_ipoptfminbnd.cpp",
		"minconNLP.hpp",
		"sci_minconNLP.cpp",
		"sci_ipoptfmincon.cpp",
		"sci_ipopt.cpp",
		"sci_LinProg.cpp",
        "sci_LinCLP.cpp",
        "LinCLP.hpp",
        "read_mps.cpp",
        "QuadTMINLP.hpp",
        "sci_QuadTMINLP.cpp",
        "cpp_intqpipopt.cpp"
	]

tbx_build_gateway(toolbox_title,Function_Names,Files,get_absolute_file_path("builder_gateway_cpp.sce"), [], Linker_Flag, C_Flags, [], "g++");

clear WITHOUT_AUTO_PUTLHSVAR toolbox_title Function_Names Files Linker_Flag C_Flags;
