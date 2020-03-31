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

function [xopt,fopt,exitflag,gradient,hessian] = intfminbnd (varargin)
	// Solves a multi-variable mixed integer non linear programming optimization problem on a bounded interval
	//
	//   Calling Sequence
	//   xopt = intfminbnd(f,intcon,x1,x2)
	//   xopt = intfminbnd(f,intcon,x1,x2,options)
	//   [xopt,fopt] = intfminbnd(.....)
	//   [xopt,fopt,exitflag]= intfminbnd(.....)
	//   [xopt,fopt,exitflag,output]=intfminbnd(.....)
	//   [xopt,fopt,exitflag,gradient,hessian]=intfminbnd(.....)
	//
	//   Parameters
	//   f : a function, representing the objective function of the problem 
	//   x1 : a vector, containing the lower bound of the variables.
	//   x2 : a vector, containing the upper bound of the variables.
  //   intcon : a vector of integers, represents which variables are constrained to be integers
	//   options : a list, containing the option for user to specify. See below for details.
	//   xopt : a vector of doubles, containing the the computed solution of the optimization problem.
	//   fopt : a scalar of double, containing the the function value at x.
	//   exitflag : a scalar of integer, containing the flag which denotes the reason for termination of algorithm. See below for details.
  //   gradient : a vector of doubles, containing the Objective's gradient of the solution.
  //   hessian  : a matrix of doubles, containing the Objective's hessian of the solution.
	//
	//   Description
	//   Search the minimum of a multi-variable mixed integer non linear programming optimization on bounded interval specified by :
	//   Find the minimum of f(x) such that 
	//
	//   <latex>
	//    \begin{eqnarray}
	//    &\mbox{min}_{x}
	//    & f(x)\\
	//    & \text{subject to} & x1 \ < x \ < x2 \\
  //    & x_i \in \!\, \mathbb{Z}, i \in \!\, I
	//    \end{eqnarray}
	//   </latex>
	//
	//   The routine calls Bonmin for solving the Bounded Optimization problem, Bonmin is a library written in C++.
	//
	// The options allows the user to set various parameters of the Optimization problem. 
	// It should be defined as type "list" and contains the following fields.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("IntegerTolerance", [---], "MaxNodes",[---], "MaxIter", [---], "AllowableGap",[---] "CpuTime", [---],"gradobj", "off", "hessian", "off" );</listitem>
  //   <listitem>IntegerTolerance : a Scalar, a number with that value of an integer is considered integer..</listitem>
  //   <listitem>MaxNodes : a Scalar, containing the Maximum Number of Nodes that the solver should search.</listitem>
	//   <listitem>CpuTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
	//   <listitem>AllowableGap : a Scalar, to stop the tree search when the gap between the objective value of the best known solution is reached.</listitem>
  //   <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
  //   <listitem>gradobj : a string, to turn on or off the user supplied objective gradient.</listitem>
  //   <listitem>hessian : a Scalar, to turn on or off the user supplied objective hessian.</listitem>
	//   <listitem>Default Values : options = list('integertolerance',1d-06,'maxnodes',2147483647,'cputime',1d10,'allowablegap',0,'maxiter',2147483647,'gradobj',"off",'hessian',"off")</listitem>
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
	// For more details on exitflag see the Bonmin documentation, go to http://www.coin-or.org/Bonmin
  //
	// Examples
	//	//Find x in R^6 such that it minimizes:
	//    //f(x)= sin(x1) + sin(x2) + sin(x3) + sin(x4) + sin(x5) + sin(x6)
	//	//-2 <= x1,x2,x3,x4,x5,x6 <= 2
	//    //Objective function to be minimised
	//    function y=f(x)
	//		y=0;
	//		for i =1:6
	//			y=y+sin(x(i));
	//		end	
	//	endfunction
	//	//Variable bounds  
	//	x1 = [-2, -2, -2, -2, -2, -2];
	//  x2 = [2, 2, 2, 2, 2, 2];
  //  intcon = [2 3 4];
	//	//Options
	//	options=list("MaxIter",[1500],"CpuTime", [100])
	//	[x,fval] =intfminbnd(f ,intcon, x1, x2, options)
  // // Press ENTER to continue
	//
	// Examples
	//	//Find x in R such that it minimizes:
	//    //f(x)= 1/x^2
	//	//0 <= x <= 1000
	//    //Objective function to be minimised
	//    function y=f(x)
	//		   y=1/x^2;
	//	  endfunction
	//	//Variable bounds  
	//	x1 = [0];
	//  x2 = [1000];
  //  intcon = [1];
	//	[x,fval,exitflag,output,lambda] =intfminbnd(f,intcon , x1, x2)
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
	//  x2 = [ %inf , %inf];
	//	//Options
  //	options=list("MaxIter",[1500],"CpuTime", [100])
  //  intcon = [1 2];
	//	[x,fval,exitflag,output,lambda] =intfminbnd(f,intcon, x1, x2, options)  
	// Authors
	// Harpreet Singh

	//To check the number of input and output arguments
 	[lhs , rhs] = argn();
	
	//To check the number of arguments given by the user
 	if ( rhs<4 | rhs>5 ) then
  		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be int [4 5] "), "intfminbnd", rhs);
  		error(errmsg);
 	end
 
	//Storing the Input Parameters  
 	fun = varargin(1);
 	intcon = varargin(2);	
 	x1 = varargin(3);
 	x2 = varargin(4);
  nbvar = size(x1,"*");

  param = list();
  //To check whether options has been entered by user   
  if ( rhs>=5  ) then
      param =varargin(5);
  end

  //To check whether the Input arguments
  Checktype("intfminbnd", fun, "fun", 1, "function");
  Checktype("intfminbnd", intcon, "intcon", 2, "constant");
  Checktype("intfminbnd", x1, "x1", 3, "constant");
  Checktype("intfminbnd", x2, "x2", 4, "constant");
  Checktype("intfminbnd", param, "options", 5, "list");  


  if(nbvar==0) then
    errmsg = msprintf(gettext("%s: x1 cannot be an empty"), "intfminbnd");
    error(errmsg);    
  end

  ///////////////// To check vectors ///////////////// 
  Checkvector("intfminbnd", x1, "x1", 3, nbvar)
  x1 = x1(:);
  Checkvector("intfminbnd", x2, "x2", 4, nbvar)
  x2 = x2(:);
  Checkvector("intfminbnd", intcon, "intcon", 2, size(intcon,"*"))
  intcon = intcon(:);

  if(~isequal(size(x1),size(x2))) then
    errmsg = msprintf(gettext("%s: x1 and x2 should be of same size"), "intfminbnd");
    error(errmsg);
  end

  for i=1:size(intcon,1)
      if(intcon(i)>nbvar) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be less than the number of variables"), "intfminbnd");
        error(errmsg);
      end

      if (intcon(i)<0) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be greater than 0 "), "intfminbnd");
        error(errmsg);
      end

      if(modulo(intcon(i),1)) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be an integer "), "intfminbnd");
        error(errmsg);
      end
  end   
  
options = list('integertolerance',1d-06,'maxnodes',2147483647,'cputime',1d10,'allowablegap',0,'maxiter',2147483647,'gradobj',"off",'hessian',"off")

  //Pushing param into default value
  
  for i = 1:(size(param))/2
    select convstr(param(2*i-1),'l')
      case 'integertolerance' then
        Checktype("intfminbnd_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(2) = param(2*i);
      case 'maxnodes' then
        Checktype("intfminbnd_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(4) = options(2*i);
      case 'cputime' then 
        Checktype("intfminbnd_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(6) = options(2*i);
      case 'allowablegap' then
        Checktype("intfminbnd_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(8) = options(2*i);
      case 'maxiter' then
        Checktype("intfminbnd_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(10) = options(2*i);
      case 'gradobj' then
        Checktype("intfminbnd_options", param(2*i), param(2*i-1), 2*i, "string");
        if(convstr(options(2*i),'l') == "on") then
          options(12) = "on"
        elseif(convstr(options(2*i),'l') == "off") then
          options(12) = "off"
        else
          error(999, 'Unknown string passed in gradobj.');
        end
      case 'hessian' then
        Checktype("intfminbnd_options", param(2*i), param(2*i-1), 2*i, "string");
        if(convstr(options(2*i),'l') == "on") then
          options(14) = "on";
        elseif(convstr(options(2*i),'l') == "off") then
          options(14) = "off";
        else
          error(999, 'Unknown string passed in hessian.');
        end
      else
        error(999, 'Unknown string argument passed.');
    end
  end 

  ///////////////// Functions Check /////////////////

  //To check the match between f (1st Parameter) and x1 (2nd Parameter)
  if(execstr('init=fun(x1)','errcatch')==21) then
    errmsg = msprintf(gettext("%s: Objective function and x1 did not match"), "intfminbnd");
    error(errmsg);
  end

  if(options(12) == "on") then
      if(execstr('[grad_y,grad_dy]=fun(x1)','errcatch')==59) then
        errmsg = msprintf(gettext("%s: Gradient of objective function is not provided"), "intfminbnd");
        error(errmsg);
      end
      Checkvector("intfminbnd_options", grad_dy, "dy", 12, nbvar);
  end

  if(options(14) == "on") then
      if(execstr('[hessian_y,hessian_dy,hessian]=fun(x1)','errcatch')==59) then
        errmsg = msprintf(gettext("%s: Gradient of objective function is not provided"), "intfminbnd");
        error(errmsg);
      end

      if ( ~isequal(size(hessian) == [nbvar nbvar]) ) then
        errmsg = msprintf(gettext("%s: Size of hessian should be nbvar X nbvar"), "intfminbnd");
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
    
   //Defining an inbuilt Objective gradient function 
    function [dy,check] = _gradf(x) 
      if (options(12) =="on") then
        try
          [y,dy]=fun(x)
          [dy,check] =  checkIsreal(dy)
        catch
          dy = 0;
          check=1;
        end
      else
          try 
            dy=numderivative(fun,x)
            [dy,check] =  checkIsreal(dy)
          catch
            dy=0;
          check=1;
          end
      end
    endfunction
    
    intconsize = size(intcon,"*");

	[xopt,fopt,exitflag] = inter_fminbnd("_f","_gradf","_gradhess",x1,x2,intcon,options,nbvar);

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
        printf("\nnObjective Function is Continuous Unbounded.\n");
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
