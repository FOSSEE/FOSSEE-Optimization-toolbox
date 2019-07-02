// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Adarsh Shah
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in


function [xopt,fopt,exitflag,iterations,output,lambda] = quadprogCLP (varargin)
	// Solves a linear quadratic problem.
	//
	//   Syntax
	//   xopt = quadprogCLP(H,f,A,b)
	//   xopt = quadprogCLP(H,f,A,b,Aeq,beq)
	//   xopt = quadprogCLP(H,f,A,b,Aeq,beq,lb,ub)
	//   [xopt,fopt,exitflag,output,lamda] = quadprogCLP( ... )
	//   
	//   Parameters
	//   H : a symmetric matrix of double, represents coefficients of quadratic in the quadratic problem.
	//   f : a vector of double, represents coefficients of linear in the quadratic problem
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
	//   b : a vector of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b.
	//   Aeq : a matrix of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   beq : a vector of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   lb : a vector of double, contains lower bounds of the variables. The default value is 0.
	//   ub : a vector of double, contains upper bounds of the variables. The default value is %inf.
	//   xopt : a vector of double, the computed solution of the optimization problem.
	//   fopt : a double, the value of the function at x.
    //   exitflag : The exit status. See below for details.
    //   iterations: Total number of itereations performed.
	//   output : The structure consist of statistics about the optimization. See below for details.
	//   lambda : The structure consist of the Lagrange multipliers at the solution of problem. See below for details.
	//   
	//   Description
	//   Search the minimum of a constrained linear quadratic optimization problem specified by :
	//
	//   <latex>
	//	\begin{eqnarray}
	//	\begin{align*}
	//  \\min\: 1/2⋅x^T⋅H⋅x + f^T⋅x  
	//	\\ subject\, to
	//	\\ A{}'x \leqslant b
	//	\\ Aeq{}'x = beq
	//	\\ lb \leqslant x \leqslant ub
	//	\end{align*}
	//	\end{eqnarray}
	//   </latex>
	//   
	//   The routine calls CLP for solving the quadratic problem, Clp is a library written in C++.
	//
	// The exitflag allows to know the status of the optimization which is given back by Clp.
	// <itemizedlist>
    // <listitem>exitflag=0 : Optimal Solution Found </listitem>
    // <listitem>exitflag=1 : Primal Infeasible</listitem>
    // <listitem>exitflag=2 : Dual Infeasible</listitem>
    // <listitem>exitflag=3 : Maximum Number of iterations exceeded</listitem>
    // <listitem>exitflag=4 : Solution Abandoned</listitem>
    // <listitem>exitflag=5 : Primal Objective Limit reached</listitem>
	// <listitem>exitflag=6 : Dual Objective Limit reached</listitem>	
	// </itemizedlist>
	// 
	// For more details on exitflag see the CLP documentation, go to http://www.coin-or.org/Clp/documentation/
	//
	// Examples
	//		//Ref : example 14 :
	//		//https://www.me.utexas.edu/~jensen/ORMM/supplements/methods/nlpmethod/S2_quadratic.pdf
	//		// min. -8*x1*x1 -16*x2*x2 + x1 + 4*x2
	//		// such that
	//		//	x1 + x2 <= 5,
	//		//	x1 <= 3,
	//		//	x1 >= 0,
	//		//	x2 >= 0
	//	H = [2 0
	//		 0 8]; 
	//	f = [-8; -16];
	//  A = [1 1;1 0];
	//	b = [5;3];
	//	lb = [0; 0];
	//	ub = [%inf; %inf];
	//	[xopt,fopt,exitflag,output,lambda] = quadprogCLP(H,f,A,b,[],[],lb,ub)
	// // Press ENTER to continue 
	//
	// Examples 
	//  //Find x in R^6 such that:
	//    Aeq= [1,-1,1,0,3,1;
	//         -1,0,-3,-4,5,6;
	//          2,5,3,0,1,0];
	//    beq=[1; 2; 3];
	//    A= [0,1,0,1,2,-1;
	//       -1,0,2,1,1,0];
	//    b = [-1; 2.5];
	//    lb=[-1000; -10000; 0; -1000; -1000; -1000];
	//    ub=[10000; 100; 1.5; 100; 100; 1000];
	//    //and minimize 0.5*x'*H*x + f'*x with
	//    f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
	//    [xopt,fopt,exitflag,output,lambda]=quadprogCLP(H,f,A,b,Aeq,beq,lb,ub)
	//
	// Examples
	//	//Solving Linear Programming Problem
	//	// min -x0 - x1
	//  // subject to 
	//	//	x0+2x1 <= 3
	//	//  2x0+x1 <= 3
	//	f = [-1;-1];
	//	A = [1,2;2,1];
	//	b = [3;3];
	//	[xopt,fopt] = quadprogCLP([],f,A,b);
	// Authors
	// Adarsh Shah
    
    
    //To check the number of input and output argument
    [lhs , rhs] = argn();

	//To check the number of argument given by user
	if ( rhs < 4 | rhs == 5 | rhs > 8 ) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [4 6 7 8]"), "quadprogCLP", rhs);
		error(errmsg)
	end

	H = [];
	f = [];


    H = varargin(1);
	f = varargin(2);
	n = size(H,1);
    

    if(n==0) then
        n= size(f,1);
    end

    if ( rhs<3 ) then
        A = []
        b = []
    else
        A = varargin(3);
        b = varargin(4);
    end

    if ( rhs<5 ) then
        Aeq = []
        beq = []
    else
        Aeq = varargin(5);
        beq = varargin(6);
    end

    if(rhs<7) then
        lb = zeros(n, 1);
    else
        lb = varargin(7); 
    end
    if(rhs<8) then
        ub = ones(n, 1)*%inf;
    else
        ub = varargin(8);
    end
         
    [xopt,fopt,exitflag,iterations,output,lambda] = sci_quadprogCLP(H,f,0,A,b,Aeq,beq,lb,ub);
    
endfunction
