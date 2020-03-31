// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh, Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function [xopt,fopt,exitflag,gradient,hessian] = intfminunc (varargin)
	// Solves an unconstrainted multi-variable mixed integer non linear programming optimization problem
	//
	//   Calling Sequence
	//   xopt = intfminunc(f,x0)
	//   xopt = intfminunc(f,x0,intcon)
	//   xopt = intfminunc(f,x0,intcon,options)
	//   [xopt,fopt] = intfminunc(.....)
	//   [xopt,fopt,exitflag]= intfminunc(.....)
	//   [xopt,fopt,exitflag,gradient,hessian]= intfminunc(.....)
	//
	//   Parameters
	//   f : a function, representing the objective function of the problem 
	//   x0 : a vector of doubles, containing the starting of variables.
	//	 intcon : a vector of integers, represents which variables are constrained to be integers
	//   options: a list, containing the option for user to specify. See below for details.
	//   xopt : a vector of doubles, the computed solution of the optimization problem.
	//   fopt : a scalar of double, the function value at x. 
	//   exitflag : a scalar of integer, containing the flag which denotes the reason for termination of algorithm. See below for details.
	//   gradient : a vector of doubles, containing the Objective's gradient of the solution.
	//   hessian  : a matrix of doubles, containing the Objective's hessian of the solution.
	//
	//   Description
	//   Search the minimum of a multi-variable mixed integer non linear programming unconstrained optimization problem specified by :
	//   Find the minimum of f(x) such that 
	//
	//   <latex>
	//    \begin{eqnarray}
	//    &\mbox{min}_{x}
	//    & f(x)
	//    & x_i \in \!\, \mathbb{Z}, i \in \!\, I
	//    \end{eqnarray}
	//   </latex>
	//
	//   The routine calls Bonmin for solving the Un-constrained Optimization problem, Bonmin is a library written in C++.
	//
	// The options allows the user to set various parameters of the Optimization problem. 
	// It should be defined as type "list" and contains the following fields.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("IntegerTolerance", [---], "MaxNodes", [---], "CpuTime", [---], "AllowableGap", [---], "MaxIter", [---]);</listitem>
	//   <listitem>IntegerTolerance : a Scalar, containing the Integer tolerance value that the solver should take.</listitem>
	//   <listitem>MaxNodes : a Scalar, containing the maximum nodes that the solver should make.</listitem>
	//	 <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
	//   <listitem>AllowableGap : a Scalar, containing the allowable gap value that the solver should take.</listitem>
	//   <listitem>CpuTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
	//   <listitem>gradobj : a string, to turn on or off the user supplied objective gradient.</listitem>
	//   <listitem>hessian : a Scalar, to turn on or off the user supplied objective hessian.</listitem>
  //   <listitem>Default Values : options = list('integertolerance',1d-06,'maxnodes',2147483647,'cputime',1d10,'allowablegap',0,'maxiter',2147483647,'gradobj',"off",'hessian',"off")</listitem>
  // </itemizedlist>
	//
	// The exitflag allows to know the status of the optimization which is given back by Bonmin.
	// <itemizedlist>
	//   <listitem>exitflag=0 : Optimal Solution Found. </listitem>
	//   <listitem>exitflag=1 : InFeasible Solution.</listitem>
	//   <listitem>exitflag=2 : Output is Continuous Unbounded.</listitem>
	//   <listitem>exitflag=3 : Limit Exceeded.</listitem>
	//   <listitem>exitflag=4 : User Interrupt.</listitem>
	//   <listitem>exitflag=5 : MINLP Error.</listitem>
	// </itemizedlist>
	//
	// For more details on exitflag see the Bonmin page, go to http://www.coin-or.org/Bonmin
	//
	// Examples
	//     //Find x in R^2 such that it minimizes the Rosenbrock function 
	//     //f = 100*(x2 - x1^2)^2 + (1-x1)^2
	//     //Objective function to be minimised
	//     function y= f(x)
	//   	    y= 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
	//     endfunction
	//	 //Starting point  
	//     x0=[-1,2];
  //     intcon = [2]
	//     //Options
	//     options=list("MaxIter", [1500], "CpuTime", [500]);
	//     //Calling 
	//     [xopt,fopt,exitflag,gradient,hessian]=intfminunc(f,x0,intcon,options)
  // // Press ENTER to continue
	//
	// Examples
	//      //Find x in R^2 such that the below function is minimum
	//      //f = x1^2 + x2^2
	//      //Objective function to be minimised
	//      function y= f(x)
	//   	     y= x(1)^2 + x(2)^2;
	//      endfunction
	//	  //Starting point  
	//      x0=[2,1];
  //      intcon = [1];  
	//      [xopt,fopt]=intfminunc(f,x0,intcon)
  // // Press ENTER to continue
	//
	// Examples
	//     //The below problem is an unbounded problem:
	//     //Find x in R^2 such that the below function is minimum
	//     //f = - x1^2 - x2^2
	//     //Objective function to be minimised
	//     function [y,g,h] = f(x)
	//        y = -x(1)^2 - x(2)^2;
  //        g = [-2*x(1),-2*x(2)];
  //        h = [-2,0;0,-2];
	//     endfunction
	//	 //Starting point  
	//     x0=[2,1];
  //     intcon = [1]
  //     options = list("gradobj","ON","hessian","on");
	//    [xopt,fopt,exitflag,gradient,hessian]=intfminunc(f,x0,intcon,options)

	//To check the number of input and output arguments
   	[lhs , rhs] = argn();
	
	//To check the number of arguments given by the user
   	if ( rhs<2 | rhs>4 ) then
    		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be int [2 3 4] "), "intfminunc", rhs);
    		error(errmsg);
   	end
 
	//Storing the 1st and 2nd Input Parameters  
   	fun = varargin(1);
   	x0 = varargin(2);
   //To add intcon
   intcon=[];
   if ( rhs >=3 ) then
			intcon = varargin(3);
	end

  param = list();
  //To check whether options has been entered by user   
 	if ( rhs>=4  ) then
      param =varargin(4);
 	end
  
  nbvar = size(x0,"*");

  ///////////////// To check whether the Input arguments /////////////////
  Checktype("intfminunc", fun, "fun", 1, "function");
  Checktype("intfminunc", x0, "x0", 2, "constant");
  Checktype("intfminunc", intcon, "intcon", 3, "constant");
  Checktype("intfminunc", param, "options", 4, "list");  
  
  ///////////////// To check x0 ///////////////// 
  Checkvector("intfminunc", x0, "x0", 2, nbvar)
  x0 = x0(:);
  Checkvector("intfminbnd", intcon, "intcon", 3, size(intcon,"*"))
  intcon = intcon(:);

  //Error Checks for intcon 
  for i=1:size(intcon,1)
      if(intcon(i)>nbvar) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be less than the number of variables"), "intfminunc");
        error(errmsg);
      end

      if (intcon(i)<0) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be greater than 0 "), "intfminunc");
        error(errmsg);
      end

      if(modulo(intcon(i),1)) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be an integer "), "intfminunc");
        error(errmsg);
      end
  end
   
  //If options has been entered, then check whether an even number of entires has been entered   
  if (modulo(size(param),2)) then
    errmsg = msprintf(gettext("%s: Size of parameters should be even"), "intfminunc");
    error(errmsg);
  end 
   
  intconSize = length(intcon);
  
  options = list('integertolerance',1d-06,'maxnodes',2147483647,'cputime',200,'allowablegap',1d-10,'maxiter',2147483647,'gradobj',"off",'hessian', "off")

  //Pushing param into default value
  
  for i = 1:(size(param))/2
    select convstr(param(2*i-1),'l')
      case 'integertolerance' then
       if (type(param(2*i))~=1) then
			errmsg = msprintf(gettext("%s: Value for integer tolerance should be a Constant"), "intfminunc");
  			error(errmsg);
		else
			options(2) = param(2*i);    //Setting the integer tolerance as per user entry
		end
      case 'maxnodes' then
        if (type(param(2*i))~=1) then
			errmsg = msprintf(gettext("%s: Value for Maximum nodes should be a Constant"), "intfminunc");
  			error(errmsg);
		else
			options(4) = param(2*i);    //Setting the Maximum nodes as per user entry
		end
      case 'cputime' then 
      	if (type(param(2*i))~=1) then
			errmsg = msprintf(gettext("%s: Value for Maximum Cpu-time should be a Constant"), "intfminunc");
  			error(errmsg);
		else
			options(6) = param(2*i);    //Setting the maximum CPU time as per user entry
		end

      case 'allowablegap' then
         if (type(param(2*i))~=1) then
			errmsg = msprintf(gettext("%s: Value for allowable gap should be a Constant"), "intfminunc");
  			error(errmsg);
		else
			options(8) = param(2*i);    //Setting the allowable gap as per user entry
		end
      case 'maxiter' then
        if (type(param(2*i))~=1) then
			errmsg = msprintf(gettext("%s: Value for max. no. of iterations should be a Constant"), "intfminunc");
  			error(errmsg);
		else
			options(10) = param(2*i);    //Setting the max. no. of iterations as per user entry
		end
      case 'gradobj' then
        Checktype("intfminbnd_options", param(2*i), param(2*i-1), 2*i, "string");
        if(convstr(param(2*i),'l') == "on") then
          options(12) = "on"
        elseif(convstr(param(2*i),'l') == "off") then
          options(12) = "off"
        else
          error(999, 'Unknown string passed in gradobj.');
        end
      case 'hessian' then
        Checktype("intfminbnd_options", param(2*i), param(2*i-1), 2*i, "string");
        if(convstr(param(2*i),'l') == "on") then
          options(14) = "on";
        elseif(convstr(param(2*i),'l') == "off") then
          options(14) = "off";
        else
          error(999, 'Unknown string passed in hessian.');
        end
      else
        error(999, 'Unknown string argument passed.');
    end
  end 

  ///////////////// Functions Check /////////////////

  //To check the match between fun (1st Parameter) and x0 (2nd Parameter)
  if(execstr('init=fun(x0)','errcatch')==21) then
    errmsg = msprintf(gettext("%s: Objective function and x0 did not match"), "intfminunc");
    error(errmsg);
  end

  if(options(12) == "on") then
    
      if(execstr('[grad_y,grad_dy]=fun(x0)','errcatch')==59) then
        errmsg = msprintf(gettext("%s: Gradient of objective function is not provided"), "intfminunc");
        error(errmsg);
      end

      Checkvector("intfminunc_options", grad_dy, "dy", 12, nbvar);
  end

  if(options(14) == "on") then
    
      if(execstr('[hessian_y,hessian_dy,hessian]=fun(x0)','errcatch')==59) then
        errmsg = msprintf(gettext("%s: Gradient of objective function is not provided"), "intfminunc");
        error(errmsg);
      end

      if ( ~isequal(size(hessian),[nbvar nbvar]) ) then
        errmsg = msprintf(gettext("%s: Size of hessian should be nbvar X nbvar"), "intfminunc");
        error(errmsg);
      end
  end

  //Converting the User defined Objective function into Required form (Error Detectable)
 	function [y,check] = _f(x)
 		try
      y=fun(x)
      [y,check] =  checkIsreal(y)
	  catch
		  y=0;
		  check=1;
	  end
   endfunction
    
   //Defining an inbuilt Objective gradient function 
    function [dy,check] = _gradf(x) 
      if (options(12) =="on") then
        try
          [y,dy]=fun(x);
          [dy,check] =  checkIsreal(dy);
        catch
          dy = 0;
          check=1;
        end
      else
        try 
          dy=numderivative(fun,x);
          [dy,check] =  checkIsreal(dy);
        catch
          dy=0;
          check=1;
        end
      end
    endfunction

  //Defining a function to calculate Hessian if the respective user entry is OFF 
  function [hessy,check]=_gradhess(x)
    if (options(14) == "on") then
      try
        [obj,dy,hessy] = fun(x)
        [hessy,check] =  checkIsreal(hessy)
      catch
        hessy = 0;
        check=1;
      end
    else
      try
        [dy,hessy]=numderivative(fun,x)
        [hessy,check] =  checkIsreal(hessy)
      catch
        hessy=0;
        check=1;
      end
    end
  endfunction

    //Calling the bonmin function for solving the above problem
	  [xopt,fopt,exitflag] = inter_fminunc("_f","_gradf","_gradhess",x0,intcon,options,intconSize,nbvar);

  	//In the cases of the problem not being solved, return NULL to the output matrices
  	if( exitflag~=0 & exitflag~=3 ) then
  		gradient = [];
        hessian = [];
    else
      [ gradient, hessian] = numderivative(_f, xopt, [], [], "blockmat");
  	end

    //To print output message
    select exitflag
        case 0 then
            printf("\nOptimal Solution Found.\n");
        case 1 then
            printf("\nInFeasible Solution.\n");
        case 2 then
            printf("\nObjective Function is Continuous Unbounded.\n");
        case 3 then
            printf("\Stop at Tiny Step.\n");
        case 4 then
            printf("\nUser Interrupt.\n");
        case 5 then
            printf("\nMINLP Error.\n");
        else
            printf("\nInvalid status returned. Notify the Toolbox authors\n");
            break;
    end
endfunction

function [y, check] = checkIsreal(x)
	if ((~isreal(x))) then
		y = 0
		check=1;
	else
        y = x;
		check=0;
	end	
endfunction
