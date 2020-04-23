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

Build_64Bits = %t;

path_builder = get_absolute_file_path('builder_gateway_cpp.sce');

if getos()=="Windows" then
//Name of All the Functions
Function_Names = [
    	
		//for opening-closing environment and checking if it is open-close
		//"sym_open",//"sci_sym_open", "csci6";
		//"sym_close",//"sci_sym_close", "csci6";
		
		//run time parameters
		//"sym_resetParams",//"sci_sym_set_defaults", "csci6";
		//"sym_setIntParam",//"sci_sym_set_int_param", "csci6";
		//"sym_getIntParam",//"sci_sym_get_int_param", "csci6";
		//"sym_setDblParam",//"sci_sym_set_dbl_param", "csci6";
		//"sym_getDblParam",//"sci_sym_get_dbl_param", "csci6";
		//"sym_setStrParam",//"sci_sym_set_str_param", "csci6";
		//"sym_getStrParam",//"sci_sym_get_str_param", "csci6";

		//problem loaders
		//"sym_loadProblemBasic",//"sci_sym_loadProblemBasic", "csci6";
		////"sym_loadProblem",//"sci_sym_loadProblem", "csci";
		
		//solve
		//"sym_solve",//"sci_sym_solve", "csci6";
		
		//post solve functions
		//"sym_getStatus",//"sci_sym_get_status", "csci6";
		//"sym_isOptimal",//"sci_sym_get_solver_status", "csci6";
		//"sym_isInfeasible",//"sci_sym_get_solver_status", "csci6";
		//"sym_isAbandoned",//"sci_sym_get_solver_status", "csci6";
		//"sym_isIterLimitReached",//"sci_sym_get_solver_status", "csci6";
		//"sym_isTimeLimitReached",//"sci_sym_get_solver_status", "csci6";
		//"sym_isTargetGapAchieved",//"sci_sym_get_solver_status", "csci6";
		//"sym_getVarSoln",//"sci_sym_getVarSoln", "csci6";
		//"sym_getObjVal",//"sci_sym_getObjVal", "csci6";
		//"sym_getIterCount",//"sci_sym_get_iteration_count", "csci6";
  
        //Linprog function
        "linearprog","sci_linearprog", "csci6";
        "rmps","sci_rmps", "csci6";
		"quadprog_CLP","sci_quadprog_CLP", "csci6";

		//QP function
		"solveqp","sci_solveqp", "csci6";

		//Unconstrained Optimization
		"solveminuncp","sci_solveminuncp", "csci6";  

		//Bounded Optimization
		"solveminbndp","sci_solveminbndp", "csci6";

		//Constrained Optimization
		"solveminconp","sci_solveminconp", "csci6";

		//Integer programming functions (CBC)
		'matintlinprog','sci_matintlinprog', 'csci6';
		'mpsintlinprog','sci_mpsintlinprog','csci6';

		//BONMIN Functions
		'solveintqp','sci_solveintqp', 'csci6';
		'inter_fminunc','cpp_intfminunc', 'csci6';
		'inter_fminbnd','cpp_intfminbnd', 'csci6';
		'inter_fmincon','cpp_intfmincon', 'csci6';
		
		 //fotversion
        "fotversion","sci_fotversion", 'csci6';

    ];

//Name of all the files to be compiled
Files = [
        "sci_iofunc.cpp",
        
 		//Symphony
		//"globals.cpp",
		//"sci_sym_openclose.cpp",
		//"sci_solver_status_query_functions.cpp",
		//"sci_sym_solve.cpp",                    
		//"sci_sym_loadproblem.cpp",
		//"sci_sym_solution.cpp",
    	//"sci_sym_get_iteration_count.cpp",
		//"sci_sym_set_variables.cpp",

		// IPOPT
		"sci_QuadNLP.cpp",
		"sci_ipoptquadprog.cpp",
		"sci_QuadNLP.cpp",

		"sci_ipoptfminunc.cpp",
		"sci_minuncNLP.cpp",

		"sci_ipoptfminbnd.cpp",
		"sci_minbndNLP.cpp",
		
		"sci_ipoptfmincon.cpp",
		"sci_minconNLP.cpp",
		
        //CLP
        "sci_LinProg.cpp",
        "read_mps.cpp",
		"sci_quadprogCLP.cpp",
		
		//Bonmin
  		'sci_minuncTMINLP.cpp',
		'cpp_intfminunc.cpp',
		'sci_minbndTMINLP.cpp',
		'cpp_intfminbnd.cpp',		
		'sci_minconTMINLP.cpp',
		'cpp_intfmincon.cpp',
		'sci_intlinprog_matrixcpp.cpp',
		'sci_QuadTMINLP.cpp',
		'sci_intquadprog.cpp'
		'sci_intlinprog_mpscpp.cpp'    

		"sci_fotversion.cpp"    
    
    ]

elseif getos()=="Darwin" then

Function_Names = [
        
		//for opening-closing environment and checking if it is open-close
		//"sym_open",//"sci_sym_open", "csci6";
		//"sym_close",//"sci_sym_close", "csci6";
		
		//run time parameters
		//"sym_resetParams",//"sci_sym_set_defaults", "csci6";
		//"sym_setIntParam",//"sci_sym_set_int_param", "csci6";
		//"sym_getIntParam",//"sci_sym_get_int_param", "csci6";
		//"sym_setDblParam",//"sci_sym_set_dbl_param", "csci6";
		//"sym_getDblParam",//"sci_sym_get_dbl_param", "csci6";
		//"sym_setStrParam",//"sci_sym_set_str_param", "csci6";
		//"sym_getStrParam",//"sci_sym_get_str_param", "csci6";

		//problem loaders
		//"sym_loadProblemBasic",//"sci_sym_loadProblemBasic", "csci6";
		////"sym_loadProblem",//"sci_sym_loadProblem";
		
		//solve
		//"sym_solve",//"sci_sym_solve", "csci6";
		
		//post solve functions
		//"sym_getStatus",//"sci_sym_get_status", "csci6";
		//"sym_isOptimal",//"sci_sym_get_solver_status", "csci6";
		//"sym_isInfeasible",//"sci_sym_get_solver_status", "csci6";
		//"sym_isAbandoned",//"sci_sym_get_solver_status", "csci6";
		//"sym_isIterLimitReached",//"sci_sym_get_solver_status", "csci6";
		//"sym_isTimeLimitReached",//"sci_sym_get_solver_status", "csci6";
		//"sym_isTargetGapAchieved",//"sci_sym_get_solver_status", "csci6";
		//"sym_getVarSoln",//"sci_sym_getVarSoln", "csci6";
		//"sym_getObjVal",//"sci_sym_getObjVal", "csci6";
		//"sym_getIterCount",//"sci_sym_get_iteration_count", "csci6";
        
        //Linprog function
    	"linearprog","sci_linearprog", "csci6";
    	"rmps","sci_rmps","csci6"; 
		"quadprog_CLP","sci_quadprog_CLP", "csci6";

		//QP function
		"solveqp","sci_solveqp", "csci6"; 
		//"sci_quadprogCLP","sci_quadprogCLP", "csci6"; 
	
		//Unconstrained Optimization
		"solveminuncp","sci_solveminuncp", "csci6"; 

		//Bounded Optimization
		"solveminbndp","sci_solveminbndp", "csci6";   

		//Constrained Optimization
		"solveminconp","sci_solveminconp", "csci6";

		//Integer programming functions (CBC)
		'matintlinprog','sci_matintlinprog', 'csci6';
		'mpsintlinprog','sci_mpsintlinprog','csci6';

		//BONMIN Functions
		'solveintqp','sci_solveintqp', 'csci6';
		'inter_fminunc','cpp_intfminunc', 'csci6';
		'inter_fminbnd','cpp_intfminbnd', 'csci6';
		'inter_fmincon','cpp_intfmincon', 'csci6';

		//fotversion
        "fotversion","sci_fotversion", 'csci6';
    ];

//Name of all the files to be compiled
Files = [
        //"sci_iofunc.cpp",

		//Symphony
		//"globals.cpp",
		//"sci_sym_openclose.cpp",
		//"sci_solver_status_query_functions.cpp",
		//"sci_sym_solve.cpp",                    
		//"sci_sym_loadproblem.cpp",
		//"sci_sym_solution.cpp",
    	//"sci_sym_get_iteration_count.cpp",
		//"sci_sym_set_variables.cpp",

        // IPOPT
		"sci_QuadNLP.cpp",
		"sci_ipoptquadprog.cpp",
		"sci_QuadNLP.cpp",

		"sci_ipoptfminunc.cpp",
		"sci_minuncNLP.cpp",

		"sci_ipoptfminbnd.cpp",
		"sci_minbndNLP.cpp",

		"sci_ipoptfmincon.cpp",
		"sci_minconNLP.cpp",

        //CLP
        "sci_LinProg.cpp",
        "read_mps.cpp",

		"sci_quadprogCLP.cpp",


		//Bonmin
  		'sci_minuncTMINLP.cpp',
		'cpp_intfminunc.cpp',
		'sci_minbndTMINLP.cpp',
		'cpp_intfminbnd.cpp',		
		'sci_minconTMINLP.cpp',
		'cpp_intfmincon.cpp',
		'sci_intlinprog_matrixcpp.cpp',
		'sci_QuadTMINLP.cpp',
		'sci_intquadprog.cpp',
		'sci_intlinprog_mpscpp.cpp',

		"sci_fotversion.cpp"
        
    ]

else
//Name of All the Functions
Function_Names = [
        
       //for opening-closing environment and checking if it is open-close
		//"sym_open",//"sci_sym_open", "csci6";
		//"sym_close",//"sci_sym_close", "csci6";
		
		//run time parameters
		//"sym_resetParams",//"sci_sym_set_defaults", "csci6";
		//"sym_setIntParam",//"sci_sym_set_int_param", "csci6";
		//"sym_getIntParam",//"sci_sym_get_int_param", "csci6";
		//"sym_setDblParam",//"sci_sym_set_dbl_param", "csci6";
		//"sym_getDblParam",//"sci_sym_get_dbl_param", "csci6";
		//"sym_setStrParam",//"sci_sym_set_str_param", "csci6";
		//"sym_getStrParam",//"sci_sym_get_str_param", "csci6";

		//problem loaders
		//"sym_loadProblemBasic",//"sci_sym_loadProblemBasic", "csci6";
		////"sym_loadProblem",//"sci_sym_loadProblem", "csci6";
		
		//solve
		//"sym_solve",//"sci_sym_solve", "csci6";
		
		//post solve functions
		//"sym_getStatus",//"sci_sym_get_status", "csci6";
		//"sym_isOptimal",//"sci_sym_get_solver_status", "csci6";
		//"sym_isInfeasible",//"sci_sym_get_solver_status", "csci6";
		//"sym_isAbandoned",//"sci_sym_get_solver_status", "csci6";
		//"sym_isIterLimitReached",//"sci_sym_get_solver_status", "csci6";
		//"sym_isTimeLimitReached",//"sci_sym_get_solver_status", "csci6";
		//"sym_isTargetGapAchieved",//"sci_sym_get_solver_status", "csci6";
		//"sym_getVarSoln",//"sci_sym_getVarSoln", "csci6";
		//"sym_getObjVal",//"sci_sym_getObjVal", "csci6";
		//"sym_getIterCount",//"sci_sym_get_iteration_count", "csci6";
        
        //Linprog function
        "linearprog","sci_linearprog", "csci6";
        "rmps","sci_rmps","csci6";   
		"quadprog_CLP","sci_quadprog_CLP", "csci6";

		//QP function
		"solveqp","sci_solveqp", "csci6";  

		//Unconstrained Optimization
		"solveminuncp","sci_solveminuncp", "csci6"; 

		//Bounded Optimization
		"solveminbndp","sci_solveminbndp", "csci6";   

		//Constrained Optimization
		"solveminconp","sci_solveminconp", "csci6";

		//Integer programming functions (CBC)
		'matintlinprog','sci_matintlinprog', 'csci6';
		'mpsintlinprog','sci_mpsintlinprog','csci6';

		//BONMIN Functions
		'solveintqp','sci_solveintqp', 'csci6';
		'inter_fminunc','cpp_intfminunc', 'csci6';
		'inter_fminbnd','cpp_intfminbnd', 'csci6';
		'inter_fmincon','cpp_intfmincon', 'csci6';

		//fotversion
        "fotversion","sci_fotversion", 'csci6';
    ];

//Name of all the files to be compiled
Files = [
        //"sci_iofunc.cpp",
        
        //Symphony
		//"globals.cpp",
		//"sci_sym_openclose.cpp",
		//"sci_solver_status_query_functions.cpp",
		//"sci_sym_solve.cpp",                    
		//"sci_sym_loadproblem.cpp",
		//"sci_sym_solution.cpp",
    	//"sci_sym_get_iteration_count.cpp",
		//"sci_sym_set_variables.cpp",

        // IPOPT
		"sci_QuadNLP.cpp",
		"sci_ipoptquadprog.cpp",
		"sci_QuadNLP.cpp",

		"sci_ipoptfminunc.cpp",
		"sci_minuncNLP.cpp",

		"sci_ipoptfminbnd.cpp",
		"sci_minbndNLP.cpp",

		"sci_ipoptfmincon.cpp",
		"sci_minconNLP.cpp",

        //CLP
        "sci_LinProg.cpp",
        "read_mps.cpp"
		"sci_quadprogCLP.cpp",

		//Bonmin
  		'sci_minuncTMINLP.cpp',
		'cpp_intfminunc.cpp',
		'sci_minbndTMINLP.cpp',
		'cpp_intfminbnd.cpp',		
		'sci_minconTMINLP.cpp',
		'cpp_intfmincon.cpp',
		'sci_intlinprog_matrixcpp.cpp',
		'sci_QuadTMINLP.cpp',
		'sci_intquadprog.cpp',
		'sci_intlinprog_mpscpp.cpp'

		"sci_fotversion.cpp"
        
    ]

end


[a, opt] = getversion();
Version = opt(2);

//Build_64Bits = %f;

if getos()=="Windows" then
    third_dir = path_builder+filesep()+'..'+filesep()+'..'+filesep()+'thirdparty';
    lib_base_dir = third_dir + filesep() + 'windows' + filesep() + 'lib' + filesep() + Version + filesep();
    //inc_base_dir = third_dir + filesep() + 'windows' + filesep() + 'include' + filesep() + 'coin';
    inc_base_dir = third_dir + filesep() + 'linux' + filesep() + 'include' + filesep() + 'coin';
    threads_dir=third_dir + filesep() + 'linux' + filesep() + 'include' + filesep() + 'pthreads-win32';
    C_Flags=['-D__USE_DEPRECATED_STACK_FUNCTIONS__  -I -w '+path_builder+' '+ '-I '+inc_base_dir+' '+'-I '+threads_dir+' ']   
    Linker_Flag  = [lib_base_dir+"libcoinblas.lib "+lib_base_dir+"libcoinlapack.lib "+lib_base_dir+"libcoinmumps.lib "+lib_base_dir+"libClp.lib "+lib_base_dir+"libipopt.lib "+lib_base_dir+"libOsi.lib "+lib_base_dir+"libOsiClp.lib "+lib_base_dir+"libCoinUtils.lib "+lib_base_dir+"libCgl.lib "+lib_base_dir+"libOsiSym.lib "+lib_base_dir+"libSym.lib "+lib_base_dir+"libCbcSolver.lib "+lib_base_dir+"libCbc.lib "+lib_base_dir+"libbonmin.lib "+lib_base_dir+"pthreadVC2.lib " ]

elseif getos()=="Darwin" then 
	third_dir = path_builder+filesep()+'..'+filesep()+'..'+filesep()+'thirdparty';
    	lib_base_dir = third_dir + filesep() + 'Mac' + filesep() + 'lib' + filesep() + Version + filesep();
    	inc_base_dir = third_dir + filesep() + 'Mac' + filesep() + 'include' + filesep() + 'coin';
    	C_Flags=["-D__USE_DEPRECATED_STACK_FUNCTIONS__ -w -fpermissive -I"+path_builder+" -I"+inc_base_dir+" -Wl,-rpath "+lib_base_dir+" "]
    	Linker_Flag = ["-L"+lib_base_dir+"libSym"+" "+"-L"+lib_base_dir+"libipopt"+" "+"-L"+lib_base_dir+"libClp"+" "+"-L"+lib_base_dir+"libOsiClp"+" "+"-L"+lib_base_dir+"libCoinUtils" + " "+"-L"+lib_base_dir+"libbonmin"]


else
    third_dir = path_builder+filesep()+'..'+filesep()+'..'+filesep()+'thirdparty';
    lib_base_dir = third_dir + filesep() + 'linux' + filesep() + 'lib' + filesep() + Version + filesep();
    inc_base_dir = third_dir + filesep() + 'linux' + filesep() + 'include' + filesep() + 'coin';
    
    C_Flags=["-D__USE_DEPRECATED_STACK_FUNCTIONS__ -w -fpermissive -I"+path_builder+" -I"+inc_base_dir+" -Wl,-rpath="+lib_base_dir+" "+"-std=gnu++11"+" -fno-zero-initialized-in-bss"]
    
    Linker_Flag = ["-L"+lib_base_dir+"libSym"+" "+"-L"+lib_base_dir+"libipopt"+" "+"-L"+lib_base_dir+"libClp"+" "+"-L"+lib_base_dir+"libOsiClp"+" "+"-L"+lib_base_dir+"libCoinUtils" ]
    
end

tbx_build_gateway(toolbox_title,Function_Names,Files,get_absolute_file_path("builder_gateway_cpp.sce"), [], Linker_Flag, C_Flags);

clear toolbox_title Function_Names Files Linker_Flag C_Flags;
