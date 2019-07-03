// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: R.Vidyadhar & Vignesh Kannan
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in


function [xopt,fopt,exitflag,output,lambda] = fminbnd (varargin)
  	// Solves a multi-variable optimization problem on a bounded interval
  	//
  	//   Calling Sequence
  	//   xopt = fminbnd(f,x1,x2)
  	//   xopt = fminbnd(f,x1,x2,options)
  	//   [xopt,fopt] = fminbnd(.....)
  	//   [xopt,fopt,exitflag]= fminbnd(.....)
  	//   [xopt,fopt,exitflag,output]=fminbnd(.....)
  	//   [xopt,fopt,exitflag,output,lambda]=fminbnd(.....)
  	//
  	//   Parameters
  	//   f : a function, representing the objective function of the problem 
  	//   x1 : a vector, containing the lower bound of the variables of size (1 X n) or (n X 1) where 'n' is the number of Variables, where n is number of Variables
  	//   x2 : a vector, containing the upper bound of the variables of size (1 X n) or (n X 1) or (0 X 0) where 'n' is the number of Variables. If x2 is empty it means upper bound is +infinity
  	//   options : a list, containing the option for user to specify. See below for details.
  	//   xopt : a vector of doubles, containing the the computed solution of the optimization problem.
  	//   fopt : a scalar of double, containing the the function value at x.
  	//   exitflag : a scalar of integer, containing the flag which denotes the reason for termination of algorithm. See below for details.
  	//   output : a structure, containing the information about the optimization. See below for details.
 	//   lambda : a structure, containing the Lagrange multipliers of lower bound and upper bound at the optimized point. See below for details.
  	//
  	//   Description
  	//   Search the minimum of a multi-variable function on bounded interval specified by :
  	//   Find the minimum of f(x) such that 
  	//
  	//   <latex>
  	//    \begin{eqnarray}
  	//    &\mbox{min}_{x}
  	//    & f(x)\\
  	//    & \text{subject to} & x1 \ < x \ < x2 \\
  	//    \end{eqnarray}
  	//   </latex>
  	//
  	//   The routine calls Ipopt for solving the Bounded Optimization problem, Ipopt is a library written in C++.
  	//
  	// The options allows the user to set various parameters of the Optimization problem. 
  	// It should be defined as type "list" and contains the following fields.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("MaxIter", [---], "CpuTime", [---], "TolX", [----]);</listitem>
	//   <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
	//   <listitem>CpuTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
	//   <listitem>TolX : a Scalar, containing the Tolerance value that the solver should take.</listitem>
	//   <listitem>Default Values : options = list("MaxIter", [3000], "CpuTime", [600], TolX, [1e-4]);</listitem>
	// </itemizedlist>
	//
	// The exitflag allows to know the status of the optimization which is given back by Ipopt.
	// <itemizedlist>
	//   <listitem>exitflag=0 : Optimal Solution Found </listitem>
	//   <listitem>exitflag=1 : Maximum Number of Iterations Exceeded. Output may not be optimal.</listitem>
	//   <listitem>exitflag=2 : Maximum CPU Time exceeded. Output may not be optimal.</listitem>
	//   <listitem>exitflag=3 : Stop at Tiny Step.</listitem>
	//   <listitem>exitflag=4 : Solved To Acceptable Level.</listitem>
	//   <listitem>exitflag=5 : Converged to a point of local infeasibility.</listitem>
	// </itemizedlist>
	//
	// For more details on exitflag see the ipopt documentation, go to http://www.coin-or.org/Ipopt/documentation/
	//
	// The output data structure contains detailed informations about the optimization process. 
	// It has type "struct" and contains the following fields.
	// <itemizedlist>
	//   <listitem>output.Iterations: The number of iterations performed during the search</listitem>
	//   <listitem>output.Cpu_Time: The total cpu-time spend during the search</listitem>
	//   <listitem>output.Objective_Evaluation: The number of Objective Evaluations performed during the search</listitem>
	//   <listitem>output.Dual_Infeasibility: The Dual Infeasiblity of the final soution</listitem>
	//	 <listitem>output.Message: The output message for the problem</listitem>
	// </itemizedlist>
	//
	// The lambda data structure contains the Lagrange multipliers at the end 
	// of optimization. In the current version the values are returned only when the the solution is optimal. 
	// It has type "struct" and contains the following fields.
	// <itemizedlist>
	//   <listitem>lambda.lower: The Lagrange multipliers for the lower bound constraints.</listitem>
	//   <listitem>lambda.upper: The Lagrange multipliers for the upper bound constraints.</listitem>
	// </itemizedlist>
	//
  	// Examples
  	//	//Find x in R^6 such that it minimizes:
  	//    //f(x)= sin(x1) + sin(x2) + sin(x3) + sin(x4) + sin(x5) + sin(x6)
  	//	//-2 <= x1,x2,x3,x4,x5,x6 <= 2
  	//    //Objective function to be minimised
  	//    function y=f(x)
  	//		y=0
  	//		for i =1:6
  	//			y=y+sin(x(i));
  	//		end	
  	//	endfunction
  	//	//Variable bounds  
  	//	x1 = [-2, -2, -2, -2, -2, -2];
  	//    x2 = [2, 2, 2, 2, 2, 2];
  	//	//Options
  	//	options=list("MaxIter",[1500],"CpuTime", [100],"TolX",[1e-6])
  	//    //Calling Ipopt
  	//	[x,fval] =fminbnd(f, x1, x2, options)
	// // Press ENTER to continue
  	//
  	// Examples
  	//	//Find x in R such that it minimizes:
  	//    //f(x)= 1/x^2
  	//	//0 <= x <= 1000
  	//    //Objective function to be minimised
  	//    function y=f(x)
  	//		y=1/x^2
  	//	endfunction
  	//	//Variable bounds  
  	//	x1 = [0];
  	//    x2 = [1000];
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,output,lambda] =fminbnd(f, x1, x2)
	// // Press ENTER to continue
  	//
  	// Examples
  	//    //The below problem is an unbounded problem:
  	//	//Find x in R^2 such that it minimizes:
  	//    //f(x)= -[(x1-1)^2 + (x2-1)^2]
  	//	//-inf <= x1,x2 <= inf
  	//    //Objective function to be minimised
  	//    function y=f(x)
  	// 		y=-((x(1)-1)^2+(x(2)-1)^2);
  	//	endfunction
  	//	//Variable bounds  
  	//	x1 = [-%inf , -%inf];
  	//    x2 = [];
  	//	//Options
  	//	options=list("MaxIter",[1500],"CpuTime", [100],"TolX",[1e-6])
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,output,lambda] =fminbnd(f, x1, x2, options)  
  	// Authors
  	// R.Vidyadhar , Vignesh Kannan
	

   	//To check the number of input and output arguments
   	[lhs , rhs] = argn();
	
   	//To check the number of arguments given by the user
   	if ( rhs<3 | rhs>4 ) then
    		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 3 or 4"), "fminbnd", rhs);
   		error(errmsg)
   	end
   
   	//Storing the 1st and 2nd Input Parameters  
   	fun = varargin(1);
   	x1 = varargin(2);
   	x2 = varargin(3);
	
	
	//To check whether the 1st Input argument (fun) is a function or not
   	if (type(fun) ~= 13 & type(fun) ~= 11) then
   		errmsg = msprintf(gettext("%s: Expected function for Objective (1st Parameter)"), "fminbnd");
   		error(errmsg);
   	end
   	
   	//Converting the User defined Objective function into Required form (Error Detectable)
   	function [y,check] = f(x)
   		if(execstr('y=fun(x)','errcatch')==32 | execstr('y=fun(x)','errcatch')==27)
			y=zeros(1,1);
			check=1;
		else
			if (isreal(y)==%F) then
				y=zeros(1,1);
				check=1;
  			else
				y=fun(x);
				check=0;
			end
		end
	endfunction
 
   	//To check whether the 2nd Input argument (x1) is a vector/scalar
   	if (type(x1) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Lower Bound Vector (2nd Parameter)"), "fminbnd");
   		error(errmsg);
  	end
  	
  	if (size(x1,2)==0) then
        errmsg = msprintf(gettext("%s: Lower Bound (2nd Parameter) cannot be empty"), "fminbnd");
   		error(errmsg);
    end
    
   	if (size(x1,1)~=1) & (size(x1,2)~=1) then
      errmsg = msprintf(gettext("%s: Lower Bound (2nd Parameter) should be a vector"), "fminbnd");
      error(errmsg); 
   	elseif (size(x1,2)==1) then
   	 	x1=x1;
   	elseif (size(x1,1)==1) then
   		x1=x1';
   	end
   	s=size(x1);
   	
   	//To check the match between f (1st Parameter) and x1 (2nd Parameter)
   	if(execstr('init=fun(x1)','errcatch')==21) then
		errmsg = msprintf(gettext("%s: Objective function and x1 did not match"), "fminbnd");
   		error(errmsg);
	end
   	
   	//To check whether the 3rd Input argument (x2) is a vector/scalar
   	if (type(x2) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Upper Bound Vector (3rd Parameter)"), "fminbnd");
   		error(errmsg);
  	end
   	
   	//To check for the correct size and data of x2 (3rd paramter) and convert it to a column vector as required by Ipopt
    if (size(x2,2)==0) then
		x2 = repmat(%inf,s(1),1);
    end
    
    
    if (size(x2,1)~=1) & (size(x2,2)~=1) then
      errmsg = msprintf(gettext("%s: Upper Bound (3rd Parameter) should be a vector"), "fminbnd");
      error(errmsg); 
    elseif(size(x2,1)~=s(1) & size(x2,2)==1) then
   		errmsg = msprintf(gettext("%s: Upper Bound and Lower Bound are not matching"), "fminbnd");
   		error(errmsg);
   	elseif(size(x2,1)==s(1) & size(x2,2)==1) then
   	 	x2=x2;
   	elseif(size(x2,1)==1 & size(x2,2)~=s(1)) then
   		errmsg = msprintf(gettext("%s: Upper Bound and Lower Bound are not matching"), "fminbnd");
   		error(errmsg);
   	elseif(size(x2,1)==1 & size(x2,2)==s(1)) then
   		x2=x2';
   	end 
    
    //To check the match between f (1st Parameter) and x2 (3rd Parameter)
   	if(execstr('init=fun(x2)','errcatch')==21) then
		errmsg = msprintf(gettext("%s: Objective function and x2 did not match"), "fminbnd");
   		error(errmsg);
	end
	
    //To check the contents of x1 and x2 (2nd & 3rd Parameter)
    
    for i = 1:s(1)
		if (x1(i) == %inf) then
		   	errmsg = msprintf(gettext("%s: Value of Lower Bound can not be infinity"), "fminbnd");
    		error(errmsg); 
  		end	
		if (x2(i) == -%inf) then
		   	errmsg = msprintf(gettext("%s: Value of Upper Bound can not be negative infinity"), "fminbnd");
    		error(errmsg); 
		end	
		if(x2(i)-x1(i)<=1e-6) then
			errmsg = msprintf(gettext("%s: Difference between Upper Bound and Lower bound should be atleast > 10^-6 for variable number=	 %d "), "fminbnd", i);
    		error(errmsg)
    	end
	end
	
   	//To check, whether options has been entered by the user  
   	if ( rhs<4 | size(varargin(4)) ==0 ) then
      		param = list(); 
   	else
      		param =varargin(4); //Storing the 3rd Input Parameter in an intermediate list named 'param'
   	end
   
   	options = list("MaxIter",[3000],"CpuTime",[600],"TolX", [1e-4]);
      
   	//To check the user entry for options and storing it
   	for i = 1:(size(param))/2
       	select convstr(param(2*i-1),'l')
    		case "maxiter" then
          			if (type(param(2*i))~=1) then
          				errmsg = msprintf(gettext("%s: Value for Maximum Iteration should be a Constant"), "fminbnd");
    	      			error(errmsg);
          			else
          				options(2) = param(2*i);    //Setting the maximum number of iterations as per user entry
          			end
       		case "cputime" then
          			if (type(param(2*i))~=1) then
          				errmsg = msprintf(gettext("%s: Value for Maximum Cpu-time should be a Constant"), "fminbnd");
    	      			error(errmsg);
          			else
          				options(4) = param(2*i);    //Setting the maximum CPU time as per user entry
          			end
        	case "tolx" then
          			if (type(param(2*i))~=1) then
          				errmsg = msprintf(gettext("%s: Value for Tolerance should be a Constant"), "fminbnd");
    	      			error(errmsg);
          			else
          				options(6) = param(2*i);    //Setting the tolerance as per user entry
          			end
    		else
    	     	 	errmsg = msprintf(gettext("%s: Unrecognized parameter name %s."), "fminbnd", param(2*i-1));
    	      		error(errmsg)
    	end
   	end

   	
   	//Defining a function to calculate Gradient or Hessian in an error deductable form.
   	function [y,check]=gradhess(x,t)
		if t==1 then
			if(execstr('y=numderivative(fun,x)','errcatch')==10000)
				y=zeros(s(1),1);
				check=1;
			else
				y=numderivative(fun,x);
				if (isreal(y)==%F) then
					y=zeros(s(1),1);
					check=1;
  				else
					check=0;
				end
			end
		else 
			if(execstr('[grad,y]=numderivative(fun,x)','errcatch')==10000)
				y=zeros(s(1),1);
				check=1;
			else
				[grad,y]=numderivative(fun,x);
				if (isreal(y)==%F) then
					y=zeros(s(1),s(1));
					check=1;
  				else
					check=0;
				end	
			end
		end
   	endfunction
   
   	//Calling the Ipopt function for solving the above problem 
	[xopt,fopt,status,iter,cpu,obj_eval,dual,zl,zu] = solveminbndp("f","gradhess",x1,x2,options);
   
   	//Calculating the values for the output
   	xopt = xopt';
   	exitflag = status;
   	output = struct("Iterations", [],"Cpu_Time",[],"Objective_Evaluation",[],"Dual_Infeasibility",[],"Message","");
   	output.Iterations = iter;
    output.Cpu_Time = cpu;
    output.Objective_Evaluation = obj_eval;
    output.Dual_Infeasibility = dual;
    lambda = struct("lower", zl,"upper",zu);

	//In the cases of the problem not being solved, return NULL to the output matrices
	if( status~=0 & status~=1 & status~=2 & status~=3 & status~=4 & status~=7 ) then
		xopt=[]
		fopt=[]
		output = struct("Iterations", [],"Cpu_Time",[],"Message","");
		output.Iterations = iter;
    	output.Cpu_Time = cpu;
		lambda = struct("lower",[],"upper",[]);
	end
	
    //To print output message
    select status
    
    	case 0 then
        	printf("\nOptimal Solution Found.\n");
        	output.Message="Optimal Solution Found";
    	case 1 then
        	printf("\nMaximum Number of Iterations Exceeded. Output may not be optimal.\n");
        	output.Message="Maximum Number of Iterations Exceeded. Output may not be optimal";
    	case 2 then
       		printf("\nMaximum CPU Time exceeded. Output may not be optimal.\n");
       		output.Message="Maximum CPU Time exceeded. Output may not be optimal";
    	case 3 then
        	printf("\nStop at Tiny Step\n");
        	output.Message="Stop at Tiny Step";
    	case 4 then
        	printf("\nSolved To Acceptable Level\n");
        	output.Message="Solved To Acceptable Level";
    	case 5 then
        	printf("\nConverged to a point of local infeasibility.\n");
        	output.Message="Converged to a point of local infeasibility";
    	case 6 then
        	printf("\nStopping optimization at current point as requested by user.\n");
        	output.Message="Stopping optimization at current point as requested by user";
    	case 7 then
        	printf("\nFeasible point for square problem found.\n");
        	output.Message="Feasible point for square problem found";
    	case 8 then 
        	printf("\nIterates diverging; problem might be unbounded.\n");
        	output.Message="Iterates diverging; problem might be unbounded";
    	case 9 then
        	printf("\nRestoration Failed!\n");
        	output.Message="Restoration Failed!";
    	case 10 then
        	printf("\nError in step computation (regularization becomes too large?)!\n");
        	output.Message="Error in step computation (regularization becomes too large?)!";
    	case 11 then
        	printf("\nProblem has too few degrees of freedom.\n");
        	output.Message="Problem has too few degrees of freedom";
    	case 12 then
        	printf("\nInvalid option thrown back by Ipopt\n");
        	output.Message="Invalid option thrown back by Ipopt";
    	case 13 then
        	printf("\nNot enough memory.\n");
        	output.Message="Not enough memory";
    	case 15 then
        	printf("\nINTERNAL ERROR: Unknown SolverReturn value - Notify Ipopt Authors.\n");
        	output.Message="INTERNAL ERROR: Unknown SolverReturn value - Notify Ipopt Authors";
    	else
        	printf("\nInvalid status returned. Notify the Toolbox authors\n");
        	output.Message="Invalid status returned. Notify the Toolbox authors";
        	break;
        end
    

endfunction
