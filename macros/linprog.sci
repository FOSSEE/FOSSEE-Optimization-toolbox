// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: Guru Pradeep Reddy, Bhanu Priya Sayal
// Organization: FOSSEE, IIT Bombay
// Email: gurupradeept@gmail.com, bhanupriyasayal@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt


function [xopt,fopt,exitflag,output,lambda] = linprog(varargin)
  // Solves a linear programming problem.
  //
  //   Calling Sequence
  //   xopt = linprog(f,A,b)
  //   xopt = linprog(f,A,b,Aeq,beq)
  //   xopt = linprog(f,A,b,Aeq,beq,lb,ub)
  //   xopt = linprog(f,A,b,Aeq,beq,lb,ub,param)
  //   [xopt,fopt,exitflag,output,lambda] = linprog( ... )
  //   
  //   Parameters
  //   f : a vector of double, represents coefficients of objective function.
  //   A : a matrix of double, represents the linear coefficients in the inequality constraints.
  //   b : a vector of double, represents the inequality constraints.
  //   Aeq : a matrix of double, represents the linear coefficients in the equality constraints.
  //   beq : a vector of double, represents the the equality constraints.
  //   lb : a vector of double, contains lower bounds of the variables.
  //   ub : a vector of double, contains upper bounds of the variables.
  //   param: used to set maximum number of iterations and gives an option to write the problem in MPS format/LP format. The default value for maximum iterations is given as 9999.        
  //   xopt : a vector of double, the computed solution of the optimization problem.
  //   fopt : a double value, the function value at xopt.
  //   exitflag : integer identifying the status of the solved problem.                                                                    0 - optimal solution                                                                                                                      1 - primal infeasible                                                                                                                     2 - dual infeasible                                                                                                                       3 - iteration limit reached                                                                                                               4 - abandoned                                                                                                                             5 - primal objective limit reached                                                                                                        6 - dual objective limit reached
  //   output : structure containing information about the optimization.
  //   lambda : structure containing the Lagrange multipliers at the solution x (separated by constraint type).
  //   
  //   Description
  //   OSI-CLP is used for solving the linear programming problems, OSI-CLP is a library written in C++.
  //   Search the minimum of a constrained linear programming problem specified by :
  //
  //   <latex>
  //    \begin{eqnarray}
  //    &\mbox{min}_{x}
  //    & f^T*x  \\
  //    & \text{subject to} & A.x \leq b \\
  //    & & Aeq.x = beq \\
  //    & & lb \leq x \leq ub \\
  //    \end{eqnarray}
  //   </latex>
  //
  //
  // Examples
  // Optimal problems
  // //Linear program, linear inequality constraints
  // f=[-1,-1/3];
  // A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
  // b=[2,1,2,1,-1,2]   
  // [xopt,fopt,exitflag,output,lambda]=linprog(f, A, b);
  // //Linear program with Linear Inequalities and Equalities
  // f=[-1,-1/3];
  // A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
  // b=[2,1,2,1,-1,2]   
  // Aeq=[1,1/4]
  // beq=[1/2]
  // [xopt,fopt,exitflag,output,lambda]=linprog(f, A, b, Aeq, beq);
  // //Linear program with all constraint types 
  // f=[-1,-1/3];
  // A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
  // b=[2,1,2,1,-1,2]   
  // Aeq=[1,1/4]
  // beq=[1/2]
  // lb=[-1,-0.5]
  // ub=[1.5,1.25]
  // [xopt,fopt,exitflag,output,lambda]=linprog(f, A, b, Aeq, beq, lb, ub);
  // //Linear program with all constraint types, and writing the given problem in MPS/LP file format
  // f=[-1,-1/3];
  // A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
  // b=[2,1,2,1,-1,2]   
  // Aeq=[1,1/4]
  // beq=[1/2]
  // lb=[-1,-0.5]
  // ub=[1.5,1.25]
  // options=list("MaxIter",[200],"WriteMps","ON","WriteLp","ON");
  // mps="../path to the mps file../.."
  // lp="../path to the lp file../.."
  // [xopt,fopt,exitflag,output,lambda]=linprog(f, A, b, Aeq, beq, lb, ub,options,mps,lp);
  // Examples 
  // Primal Infeasible Problem
  // f=[-1,-1,-1];
  // A=[1,2,-1];
  // b=[-4];
  // Aeq=[1,5,3;1,1,0];
  // beq=[10,100];
  // lb=[0,0,0];
  // ub=[%inf,%inf,%inf];
  // [xopt,fopt,exitflag,output,lambda]= linprog(f,A,b,Aeq,beq,lb,ub);
  // Examples
  // Dual Infeasible Problem
  // f=[3,5,-7];
  // A=[-1,-1,4;1,1,4];
  // b=[-8,5];
  // Aeq=[];
  // beq=[];
  // lb=[-%inf,-%inf,-%inf];
  // ub=[%inf,%inf,%inf];
  // [xopt,fopt,exitflag,output,lambda]= linprog(f,A,b,Aeq,beq,lb,ub);
  // Authors
  // Bhanu Priya Sayal, Guru Pradeep Reddy
    
    
//To check the number of input and output argument
   [lhs , rhs] = argn();
	
//To check the number of argument given by user
   if ( rhs < 3 | rhs == 4 | rhs == 6 | rhs >10 ) then
    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [3 5 7 8 9]"), "linprog", rhs);
    error(errmsg)
   end
   
   f = varargin(1);
   if(size(f,1)>size(f,2)) then
   		f=f';
   end
   nbVar = size(f,2);
   A = varargin(2);
   b = varargin(3);
  // nbVar = size(f,2);
   if(size(b,1)>size(b,2)) then
   		b=b';
   end
   if ( rhs<4 ) then
      Aeq = []
      beq = []
   else
      Aeq = varargin(4);
      beq = varargin(5);
  end
  if(size(beq,1)>size(beq,2)) then
   		beq=beq';
   end
  
  if ( rhs<6 ) then
      LB = [];
      UB = [];
   else
      LB = varargin(6);
      UB = varargin(7);
  end
  if(size(LB,1)>size(LB,2)) then
   		LB=LB';
   end
  if(size(UB,1)>size(UB,2)) then
   		UB=UB';
   end

   if ( rhs<8 | size(varargin(8)) ==0 ) then
      param = list();
   else
      param =varargin(8);
   end
   
   if (size(LB,2)==0) then
        LB = repmat(-%inf,1,nbVar);
    end
    
    if (size(UB,2)==0) then
        UB = repmat(%inf,1,nbVar);
    end

    if (size(f,2)==0) then
        f = repmat(0,1,nbVar);
    end

    if (type(param) ~= 15) then
      errmsg = msprintf(gettext("%s: options should be a list "), "linprog");
      error(errmsg);
    end
   

   if (modulo(size(param),2)) then
   errmsg = msprintf(gettext("%s: Size of parameters should be even"), "linprog");
   error(errmsg);
   end

    options = list(..
      "MaxIter"     , [3000],);
   flag_lp = 0;
   flag_mps = 0;
   file_lp = [];
   file_mps = [];
   for i = 1:(size(param))/2
        
      	select param(2*i-1)
    	case "MaxIter" then
       		options(2*i) = param(2*i);
       	case "WriteMps" then
			if(param(2*i)~="ON" & param(2*i)~="OFF")
            	errmsg=msprintf(gettext("%s: The input parameter should be ON or OFF"), "linprog");
                        error(errmsg);
            end
         	if(param(2*i)=="ON")
					//To check if the user has provided the file path
                    if(rhs<9) then
                    	errmsg=msprintf(gettext("%s: File name is missing"), "linprog");
                        error(errmsg);
                    else 
                    //This flag1 is activated(ie. =1) if WriteMps is set
                    flag_mps=1;
                    file_mps = varargin(9);
                    end
             end
        case "WriteLp" then
			if(param(2*i)~="ON" & param(2*i)~="OFF")
            	errmsg=msprintf(gettext("%s: The input parameter should be ON or OFF"), "linprog");
                        error(errmsg);
            end
         	if(param(2*i)=="ON") then
					//To check if the user has provided the file path
                    if(flag_mps==1) then	
                    	if(rhs<10) then
                    		errmsg=msprintf(gettext("%s: File name is missing"), "linprog");
                        	error(errmsg);
                    	else 
                       //This flag1 is activated(ie. =1) if WriteMps is set
                       flag_lp=1;
                       file_lp = varargin(10);
                       end
                    else
                        if(rhs<9) then
                        	errmsg=msprintf(gettext("%s: File name is missing"), "linprog");
                            error(errmsg);
                        else
                            flag_lp=1;
                            file_lp = varargin(9);
                        end
                    end
             end
		end
     end

   nbConInEq = size(A,1);
   nbConEq = size(Aeq,1);

   //Check the size of inequality constraint which should be equal to the number of variables
   if ( size(A,2) ~= nbVar & size(A,2) ~= 0) then
      errmsg = msprintf(gettext("%s: The number of columns in A must be the same as the number of elements of f"), "linprog");
      error(errmsg);
   end

   //Check the size of equality constraint which should be equal to the number of variables
   if ( size(Aeq,2) ~= nbVar & size(Aeq,2) ~= 0 ) then
      errmsg = msprintf(gettext("%s: The number of columns in Aeq must be the same as the number of elements of f"), "linprog");
      error(errmsg);
   end


   //Check the size of Lower Bound which should be equal to the number of variables
   if ( size(LB,2) ~= nbVar) then
      errmsg = msprintf(gettext("%s: The Lower Bound is not equal to the number of variables"), "linprog");
      error(errmsg);
   end

   //Check the size of Upper Bound which should equal to the number of variables
   if ( size(UB,2) ~= nbVar) then
      errmsg = msprintf(gettext("%s: The Upper Bound is not equal to the number of variables"), "linprog");
      error(errmsg);
   end

   //Check the size of constraints of Lower Bound which should equal to the number of constraints
   if ( size(b,2) ~= nbConInEq & size(b,2) ~= 0) then
      errmsg = msprintf(gettext("%s: The number of rows in A must be the same as the number of elements of b"), "linprog");
      error(errmsg);
   end

   //Check the size of constraints of Upper Bound which should equal to the number of constraints
   if ( size(beq,2) ~= nbConEq & size(beq,2) ~= 0) then
      errmsg = msprintf(gettext("%s: The number of rows in Aeq must be the same as the number of elements of beq"), "linprog");
      error(errmsg);
   end

   //Check if the user gives a matrix instead of a vector
   
   if (size(LB,1)~=1)& (size(LB,2)~=1) then
      errmsg = msprintf(gettext("%s: Lower Bound should be a vector"), "linprog");
      error(errmsg); 
   end
   
   if (size(UB,1)~=1)& (size(UB,2)~=1) then
      errmsg = msprintf(gettext("%s: Upper Bound should be a vector"), "linprog");
      error(errmsg); 
   end
   
   if (nbConInEq) then
        if ((size(b,1)~=1)& (size(b,2)~=1)) then
            errmsg = msprintf(gettext("%s: Constraint Lower Bound should be a vector"), "linprog");
            error(errmsg); 
        end
    end
    
    if (nbConEq) then
        if (size(beq,1)~=1)& (size(beq,2)~=1) then
            errmsg = msprintf(gettext("%s: Constraint should be a vector"), "linprog");
            error(errmsg); 
        end
   end
  
	for i = 1:nbConInEq
		if (b(i) == -%inf)
		   	errmsg = msprintf(gettext("%s: Value of b can not be negative infinity"), "linprog");
            error(errmsg); 
        end	
	end
    
	for i = 1:nbConEq
		if (beq(i) == -%inf)
		   	errmsg = msprintf(gettext("%s: Value of beq can not be negative infinity"), "linprog");
            error(errmsg); 
        end	
	end
   
   nbCon = size(A,1);
   nbEqCon= size(Aeq,1);
   
   [xopt,fopt,status,violation,iter,Zl,ineq,eq] = linearprog(f,A,b,Aeq,beq,LB,UB,options,flag_mps,flag_lp,file_mps,file_lp);
   
   xopt = xopt';
   exitflag = status;
   output = struct("Iterations"      , [],..
                   "constrviolation"	, []);
                  
   output.Iterations = iter;
   output.constrviolation = violation
   lambda = struct("reduced_cost"           , [], ..
                   "ineqlin"           , [], ..
                   "eqlin"      , []);
   
   lambda.reduced_cost = Zl;
   lambda.ineqlin =ineq;
   lambda.eqlin = eq;
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
