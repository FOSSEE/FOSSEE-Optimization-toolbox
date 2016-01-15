// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: Bhanu Priya Sayal, Guru Pradeep Reddy
// Organization: FOSSEE, IIT Bombay
// Email:bhanupriyasayal@gmail.com,gurupradeept@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt


function [xopt,fopt,exitflag,output,lambda] =mps_linprog(varargin)
  // Solves a linear programming problem given in mps format.
  //
  //   Calling Sequence
  //   xopt = mps_linprog(file)
  //   xopt = mps_linprog(file,param)
  //   [xopt, fopt, exitflag, output, lambda] = mps_linprog(file)
  //   [xopt, fopt, exitflag, output, lambda] = mps_linprog(file,param)
  //   
  //  Parameters
  //   file : a string describing the path to the mps file.
  //   param: used to set maximum number of iterations. The default value for maximum iterations is given as 9999.
  //   xopt : a vector of double, the computed solution of the optimization problem.
  //   fopt : a double, the function value at xopt.
  //   exitflag : integer identifying the status of the solved problem.                                                                    0 - optimal solution                                                                                                                      1 - primal infeasible                                                                                                                     2 - dual infeasible                                                                                                                       3 - iteration limit reached                                                                                                               4 - abandoned                                                                                                                             5 - primal objective limit reached                                                                                                        6 - dual objective limit reached
  //   output : structure containing information about the optimization.
  //   lambda : structure containing the Lagrange multipliers at the solution x (separated by constraint type).
  //
  //   Description
  //   OSI-CLP is used for solving the linear programming problem, OSI-CLP is a library written in C++. It reads the data for linear programming problem from the MPS file and solves the problem.
  //   Search the minimum of a constrained linear programming problem specified by : 
  //
  //   <latex>
  //    \begin{eqnarray}
  //    &\mbox{min}_{x}
  //    & f^T*x  \\
  //    & \text{subject to} & A.x \leq b \\
  //    & & Aeq = eq beq \\
  //    & & lb \leq x \leq ub \\
  //    \end{eqnarray}
  //   </latex>
  //   
  //
  // Examples
  // file1_link="https://github.com/bhanusayal/symphony/blob/master/demos/p0033.mps";
  // [xopt,fopt,exitflag,output,lambda] =lp_linprog(file1);
  // Examples
  // file2_link="https://github.com/bhanusayal/symphony/blob/master/demos/exmip1.mps";
  // [xopt,fopt,exitflag,output,lambda] =lp_linprog(file2);
  // Authors
  // Bhanu Priya Sayal, Guru Pradeep Reddy
    
//To check the number of input and output argument
   [lhs , rhs] = argn();
	
//To check the number of argument given by user
   if ( rhs > 2) then
    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [1 2]"), "mps_linprog",   		rhs);
    error(errmsg)
   end
   file = varargin(1);
   if ( rhs<2 | size(varargin(2)) ==0 ) then
      param = list();
   else
       param =varargin(2);
   end 
   if (type(param) ~= 15) then
      errmsg = msprintf(gettext("%s: options should be a list "), "mps_linprog");
      error(errmsg);
   end
   if (modulo(size(param),2)) then
   errmsg = msprintf(gettext("%s: Size of parameters should be even"), "mps_linprog");
   error(errmsg);
   end
   options = list(..
      "MaxIter"     , [3000],);
   for i = 1:(size(param))/2
         select param(2*i-1)
            case "MaxIter" then
        		options(2*i) = param(2*i);
         end 
   end
   //Calling the function by passing the required parameters 
   
   [xopt,fopt,status,iter,Zl,Dl] = rmps(file,options);
   
   xopt = xopt;
   fopt=fopt;
   exitflag = status;
   output = struct("Iterations"      , []);
   output.Iterations = iter;
   lambda = struct("reduced_cost"           , [], ..
                   "dual"           ,[]);
   
    lambda.reduced_cost = Zl;
    lambda.dual =Dl;
  select status

  case 0 then
     printf("\nOptimal Solution.\n");
  case 1 then 
     printf("\nPrimal Infeasible.\n");
  case 2 then 
     printf("\nDual Infeasible.\n");
  case 3 then
     printf("\nIteration limit reached.\n");
  case 4 then 
     printf("\nNumerical Difficulties.\n");
  case 5 then
     printf("\nPrimal Objective limit reached.\n");
  case 6 then
     printf("\nDual Objective limit reached.\n");
  end

endfunction
