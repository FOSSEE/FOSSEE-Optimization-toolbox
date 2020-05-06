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

function [xopt,resnorm,residual,exitflag,output,lambda,gradient] = lsqnonlin (varargin)
	// Solves a non linear data fitting problems.
	//
	//   Calling Sequence
	//   xopt = lsqnonlin(fun,x0)
	//   xopt = lsqnonlin(fun,x0,lb,ub)
	//   xopt = lsqnonlin(fun,x0,lb,ub,options)
	//   [xopt,resnorm] = lsqnonlin( ... )
	//   [xopt,resnorm,residual] = lsqnonlin( ... )
	//   [xopt,resnorm,residual,exitflag] = lsqnonlin( ... )
	//   [xopt,resnorm,residual,exitflag,output,lambda,gradient] = lsqnonlin( ... )
	//   
	//   Parameters
  	//   fun : a function, representing the objective function and gradient (if given) of the problem
	//   x0 : a vector of double, contains initial guess of variables.
	//   lb : a vector of double, contains lower bounds of the variables.
	//   ub : a vector of double,  contains upper bounds of the variables.
	//   options : a list containing the parameters to be set.
	//   xopt : a vector of double, the computed solution of the optimization problem.
	//   resnorm : a double, objective value returned as the scalar value i.e. sum(fun(x).^2).
	//   residual : a vector of double, solution of objective function i.e. fun(x).
	//   exitflag : The exit status. See below for details.
	//   output : The structure consist of statistics about the optimization. See below for details.
	//   lambda : The structure consist of the Lagrange multipliers at the solution of problem. See below for details.
  	//   gradient : a vector of doubles, containing the Objective's gradient of the solution.
	//   
	//   Description
	//   Search the minimum of a constrained non-linear least square problem specified by :
	//
	//   <latex>
	//    \begin{eqnarray}
	//    &\mbox{min}_{x}
	//    & (f_1(x)^2 + f_2(x)^2 + ... + f_n(x)^2) \\
	//    & lb \leq x \leq ub \\
	//    \end{eqnarray}
	//   </latex>
	//   
	//   The routine calls fmincon which calls Ipopt for solving the non-linear least square problem, Ipopt is a library written in C++.
	//
  	// The options allows the user to set various parameters of the Optimization problem. 
  	// It should be defined as type "list" and contains the following fields.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("MaxIter", [---], "CpuTime", [---],"GradObj", "on");</listitem>
	//   <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
	//   <listitem>CpuTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
	//   <listitem>GradObj : a string, representing the gradient function is on or off.</listitem>
	//   <listitem>Default Values : options = list("MaxIter", [3000], "CpuTime", [600], "GradObj", "off");</listitem>
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
	//   <listitem>output.iterations: The number of iterations performed during the search</listitem>
	//   <listitem>output.constrviolation: The max-norm of the constraint violation.</listitem>
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
	// //A simple non-linear least square example taken from leastsq default present in scilab
	//	function y=yth(t, x)
	//	   y  = x(1)*exp(-x(2)*t)
	//	endfunction
	//	// we have the m measures (ti, yi):
	//	m = 10;
	//	tm = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]';
	//	ym = [0.79, 0.59, 0.47, 0.36, 0.29, 0.23, 0.17, 0.15, 0.12, 0.08]';
	//	// measure weights (here all equal to 1...)
	//	wm = ones(m,1);
	//	// and we want to find the parameters x such that the model fits the given
	//	// data in the least square sense:
	//	//
	//	//  minimize  f(x) = sum_i  wm(i)^2 ( yth(tm(i),x) - ym(i) )^2
	//	// initial parameters guess
	//	x0 = [1.5 ; 0.8];
	//	// in the first examples, we define the function fun and dfun
	//	// in scilab language
	//	function y=myfun(x, tm, ym, wm)
	//	   y = wm.*( yth(tm, x) - ym )
	//	endfunction
	//	// the simplest call
	//	[xopt,resnorm,residual,exitflag,output,lambda,gradient] = lsqnonlin(myfun,x0)
	// // Press ENTER to continue 
	//    
	// Examples
	// //A basic example taken from leastsq default present in scilab with gradient
	//	function y=yth(t, x)
	//	   y  = x(1)*exp(-x(2)*t)
	//	endfunction
	//	// we have the m measures (ti, yi):
	//	m = 10;
	//	tm = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]';
	//	ym = [0.79, 0.59, 0.47, 0.36, 0.29, 0.23, 0.17, 0.15, 0.12, 0.08]';
	//	// measure weights (here all equal to 1...)
	//	wm = ones(m,1);
	//	// and we want to find the parameters x such that the model fits the given
	//	// data in the least square sense:
	//	//
	//	//  minimize  f(x) = sum_i  wm(i)^2 ( yth(tm(i),x) - ym(i) )^2
	//	// initial parameters guess
	//	x0 = [1.5 ; 0.8];
	//	// in the first examples, we define the function fun and dfun
	//	// in scilab language
	//	function [y,dy]=myfun(x, tm, ym, wm)
	//	   y = wm.*( yth(tm, x) - ym )
	//	   v = wm.*exp(-x(2)*tm)
	//     dy = [v , -x(1)*tm.*v]
	//	endfunction
	//	options = list("GradObj", "on")
	//	[xopt,resnorm,residual,exitflag,output,lambda,gradient] = lsqnonlin(myfun,x0,[],[],options)
	// Authors
	// Harpreet Singh


	//To check the number of input and output argument
	[lhs , rhs] = argn();

	//To check the number of argument given by user
	if ( rhs < 2 | rhs == 3 | rhs > 5  ) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [2 4 5]"), "lsqnonlin", rhs);
		error(errmsg)
	end

// Initializing all the values to empty matrix
	

	_fun = varargin(1);
	x0 = varargin(2);
	nbVar = size(x0,'*');

	if(nbVar == 0) then
		errmsg = msprintf(gettext("%s: Cannot determine the number of variables because input initial guess is empty"), "lsqcurvefit");
		error(errmsg);
	end

    if (size(x0,2)~=1) then
		errmsg = msprintf(gettext("%s: Initial Guess should be a column matrix"), "lsqcurvefit");
		error(errmsg);
    end

	if ( rhs<3 ) then
		lb = repmat(-%inf,nbVar,1);
		ub = repmat(%inf,nbVar,1);
	else
		lb = varargin(3);
		ub = varargin(4);
	end

	if ( rhs<5 | size(varargin(5)) ==0 ) then
		param = list();
	else
		param =varargin(5);
	end

	if (size(lb,2)==0) then
		lb = repmat(-%inf,nbVar,1);
	end

	if (size(ub,2)==0) then
		ub = repmat(%inf,nbVar,1);
	end

	//Check type of variables
	fot_Checktype("lsqnonlin", _fun, "fun", 1, "function")
	fot_Checktype("lsqnonlin", x0, "x0", 2, "constant")
	fot_Checktype("lsqnonlin", lb, "lb", 3, "constant")
	fot_Checktype("lsqnonlin", ub, "ub", 4, "constant")
	fot_Checktype("lsqnonlin", param, "param", 5, "list")

	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of parameters should be even"), "lsqnonlin");
		error(errmsg);
	end

	options = list(	"MaxIter"   , [3000], ...
					"CpuTime"   , [600], ...
					"GradObj"	, ["off"]);

	for i = 1:(size(param))/2

		select convstr(param(2*i-1),'l')
			case "maxiter" then
				options(2) = param(2*i);
				fot_Checktype("lsqcurvefit", param(2*i), "maxiter", i, "constant")
			case "cputime" then
				options(4) = param(2*i);
				fot_Checktype("lsqcurvefit", param(2*i), "cputime", i, "constant")
        	case "gradobj" then
        		if (convstr(param(2*i))=="on") then
    				options(6) = "on";
        		else
					errmsg = msprintf(gettext("%s: Unrecognized String [%s] entered for gradobj."), "lsqnonlin",param(2*i));
      				error(errmsg);
        		end
			else
				errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "lsqnonlin", param(2*i-1));
				error(errmsg)
		end
	end

	// Check if the user gives row vector 
	// and Changing it to a column matrix

	if (size(lb,2)== [nbVar]) then
		lb = lb(:);
	end

	if (size(ub,2)== [nbVar]) then
		ub = ub(:);
	end

  	//To check the match between fun (1st Parameter) and x0 (2nd Parameter)
   	if(execstr('init=_fun(x0)','errcatch')==21) then
		errmsg = msprintf(gettext("%s: Objective function and x0 did not match"), "lsqnonlin");
   		error(errmsg);
	end
	
    ierr = execstr('init=_fun(x0)', "errcatch")
    if ierr <> 0 then
    lamsg = lasterror();
    lclmsg = "%s: Error while evaluating the function: ""%s""\n";
    error(msprintf(gettext(lclmsg), "lsqnonlin", lamsg));
    end
        
	//Check the size of Lower Bound which should be equal to the number of variables
	if ( size(lb,1) ~= nbVar) then
		errmsg = msprintf(gettext("%s: The Lower Bound is not equal to the number of variables"), "lsqnonlin");
		error(errmsg);
	end

	//Check the size of Upper Bound which should equal to the number of variables
	if ( size(ub,1) ~= nbVar) then
		errmsg = msprintf(gettext("%s: The Upper Bound is not equal to the number of variables"), "lsqnonlin");
		error(errmsg);
	end

	//Check if the user gives a matrix instead of a vector

	if (size(lb,1)~=1)& (size(lb,2)~=1) then
		errmsg = msprintf(gettext("%s: Lower Bound should be a vector"), "lsqnonlin");
		error(errmsg); 
	end

	if (size(ub,1)~=1)& (size(ub,2)~=1) then
		errmsg = msprintf(gettext("%s: Upper Bound should be a vector"), "lsqnonlin");
		error(errmsg); 
	end

	for i = 1:nbVar
		if(lb(i)>ub(i))
			errmsg = msprintf(gettext("%s: Problem has inconsistent variable bounds"), "lsqnonlin");
			error(errmsg);
		end
	end

	//Converting it into fmincon Problem
	function [y] = __fun(x)
		[yVal] = _fun(x);
		y = sum(yVal.^2);
	endfunction

	if (options(6) == "on")
        ierr = execstr('[initx initdx]=_fun(x0)', "errcatch")
        if ierr <> 0 then
        lamsg = lasterror();
        lclmsg = "%s: Error while evaluating the function: ""%s""\n";
        error(msprintf(gettext(lclmsg), "lsqnonlin", lamsg));
        end
		function [dy] = __fGrad(x)
			[y_user,dy_user] = _fun(x)
			dy = 2*y_user'*dy_user;
			clear y_user dy_user;
		endfunction
	else
		function [dy] = __fGrad(x)
			y_user = _fun(x);
			dy_nd = numderivative(_fun,x);
			dy = 2*y_user'*dy_nd;
			clear y_user dy_nd;
		endfunction
	end
	
    options(6) = __fGrad;
	
    function [c, ceq] = nlc(x)
        c = [];
        ceq = [];
	endfunction

	[xopt,fopt,exitflag,output,lambda_fmincon,gradient] = fmincon(__fun,x0,[],[],[],[],lb,ub,nlc,options);
	
	lambda = struct("lower", lambda_fmincon.lower,"upper",lambda_fmincon.upper);
	
	resnorm = sum(_fun(xopt).^2);
	residual = _fun(xopt);

endfunction
