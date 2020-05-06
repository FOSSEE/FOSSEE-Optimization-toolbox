// Copyright (C) 2016 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function [xopt,fopt,status,output] = cbcmatrixintlinprog (varargin)
    // Sci file wrapper for the matrix_cbcintlinprog.cpp file
        
    // To check the number of input and output argument
    [lhs , rhs] = argn();
    
    // To check the number of argument given by user
    if ( rhs < 4 | rhs == 5 | rhs == 7 | rhs > 9 ) then
        errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set [4 6 8 9]"), "cbcintlinprog", rhs);
        error(errmsg);
    end
   
    c = [];
    intcon = [];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    options = list();
    
    c = varargin(1)
    intcon = varargin(2)
    A = varargin(3)
    b = varargin(4)

    if(size(c,2) == 0) then
        errmsg = msprintf(gettext("%s: Cannot determine the number of variables because input objective coefficients is empty"),"cbcintlinprog");
        error(errmsg);
    end

   if (size(c,2)~=1) then
    errmsg = msprintf(gettext("%s: Objective Coefficients should be a column matrix"), "cbcintlinprog");
    error(errmsg);
   end

   nbVar = size(c,1);

   if ( rhs<5 ) then
      Aeq = []
      beq = []
   else
      Aeq = varargin(5);
      beq = varargin(6);
   end
   
   if ( rhs<7 ) then
      lb = repmat(-%inf,1,nbVar);
      ub = repmat(%inf,1,nbVar);
   else
      lb = varargin(7);
      ub = varargin(8);
   end
   
   if (rhs<9|size(varargin(9))==0) then
      options = list();
   else
      options = varargin(9);
   end

    //Check type of variables
    fot_Checktype("cbcintlinprog", c, "c", 1, "constant")
    fot_Checktype("cbcintlinprog", intcon, "intcon", 2, "constant")
    fot_Checktype("cbcintlinprog", A, "A", 3, "constant")
    fot_Checktype("cbcintlinprog", b, "b", 4, "constant")
    fot_Checktype("cbcintlinprog", Aeq, "Aeq", 5, "constant")
    fot_Checktype("cbcintlinprog", beq, "beq", 6, "constant")
    fot_Checktype("cbcintlinprog", lb, "lb", 7, "constant")
    fot_Checktype("cbcintlinprog", ub, "ub", 8, "constant")

    // Check if the user gives empty matrix
    if (size(lb,2)==0) then
        lb = repmat(-%inf,nbVar,1);
    end
    
    if (size(intcon,2)==0) then
        intcon = [];
    end
    
    if (size(ub,2)==0) then
        ub = repmat(%inf,nbVar,1);
    end

    // Calculating the size of equality and inequality constraints
    nbConInEq = size(A,1);
    nbConEq = size(Aeq,1);

    // Check if the user gives row vector 
    // and Changing it to a column matrix

    if (size(lb,2)== [nbVar]) then
        lb = lb';
    end

    if (size(ub,2)== [nbVar]) then
        ub = ub';
    end

    if (size(b,2)== [nbConInEq]) then
        b = b';
    end

    if (size(beq,2)== [nbConEq]) then
        beq = beq';
    end

    for i=1:size(intcon,2)
        if(intcon(i)>nbVar) then
            errmsg = msprintf(gettext("%s: The values inside intcon should be less than the number of variables"), "cbcintlinprog");
            error(errmsg);
        end

        if (intcon(i)<0) then
            errmsg = msprintf(gettext("%s: The values inside intcon should be greater than 0 "), "cbcintlinprog");
            error(errmsg);
        end

        if(modulo(intcon(i),1)) then
            errmsg = msprintf(gettext("%s: The values inside intcon should be an integer "), "cbcintlinprog");
            error(errmsg);
        end
    end

   //Check the size of inequality constraint which should equal to the number of inequality constraints
    if ( size(A,2) ~= nbVar & size(A,2) ~= 0) then
        errmsg = msprintf(gettext("%s: The size of inequality constraint is not equal to the number of variables"), "cbcintlinprog");
        error(errmsg);
    end


   //Check the size of lower bound of inequality constraint which should equal to the number of constraints
    if ( size(b,1) ~= size(A,1)) then
        errmsg = msprintf(gettext("%s: The Lower Bound of inequality constraint is not equal to the number of constraint"), "cbcintlinprog");
        error(errmsg);
    end

   //Check the size of equality constraint which should equal to the number of inequality constraints
    if ( size(Aeq,2) ~= nbVar & size(Aeq,2) ~= 0) then
        errmsg = msprintf(gettext("%s: The size of equality constraint is not equal to the number of variables"), "cbcintlinprog");
        error(errmsg);
    end

   //Check the size of upper bound of equality constraint which should equal to the number of constraints
    if ( size(beq,1) ~= size(Aeq,1)) then
        errmsg = msprintf(gettext("%s: The equality constraint upper bound is not equal to the number of equality constraint"), "cbcintlinprog");
        error(errmsg);
    end

   //Check the size of Lower Bound which should equal to the number of variables
   if ( size(lb,1) ~= nbVar) then
      errmsg = msprintf(gettext("%s: The Lower Bound is not equal to the number of variables"), "cbcintlinprog");
      error(errmsg);
   end

   //Check the size of Upper Bound which should equal to the number of variables
   if ( size(ub,1) ~= nbVar) then
      errmsg = msprintf(gettext("%s: The Upper Bound is not equal to the number of variables"), "cbcintlinprog");
      error(errmsg);
   end
   
   if (type(options) ~= 15) then
      errmsg = msprintf(gettext("%s: Options should be a list "), "cbcintlinprog");
      error(errmsg);
   end
   

   if (modulo(size(options),2)) then
    errmsg = msprintf(gettext("%s: Size of parameters should be even"), "cbcintlinprog");
    error(errmsg);
   end

   //Check if the user gives a matrix instead of a vector
   
   if (((size(intcon,1)~=1)& (size(intcon,2)~=1))&(size(intcon,2)~=0)) then
      errmsg = msprintf(gettext("%s: intcon should be a vector"), "cbcintlinprog");
      error(errmsg); 
   end
   
   if (size(lb,1)~=1)& (size(lb,2)~=1) then
      errmsg = msprintf(gettext("%s: Lower Bound should be a vector"), "cbcintlinprog");
      error(errmsg); 
   end
   
   if (size(ub,1)~=1)& (size(ub,2)~=1) then
      errmsg = msprintf(gettext("%s: Upper Bound should be a vector"), "cbcintlinprog");
      error(errmsg); 
   end
   
   if (nbConInEq) then
        if ((size(b,1)~=1)& (size(b,2)~=1)) then
            errmsg = msprintf(gettext("%s: Constraint Lower Bound should be a vector"), "cbcintlinprog");
            error(errmsg); 
        end
    end
    
    if (nbConEq) then
        if (size(beq,1)~=1)& (size(beq,2)~=1) then
            errmsg = msprintf(gettext("%s: Constraint Upper Bound should be a vector"), "cbcintlinprog");
            error(errmsg); 
        end
   end
  

    //Changing the inputs in symphony's format 
    conMatrix = [A;Aeq]
    nbCon = int32(size(conMatrix,1));
    conLB = [repmat(-%inf,size(A,1),1);beq];
    conUB = [b;beq] ;
    
    isInt = repmat(%f,1,nbVar);
    // Changing intcon into column vector
    intcon = intcon(:);
    for i=1:size(intcon,1)
        isInt(intcon(i)) = %t
    end
    
    objSense = 1.0;

    //Changing into row vector
    lb = lb';
    ub = ub';
    c = c';
    
    //Pusing options as required to a double array
    optval = [];
    if length(options) == 0 then
        optval = [0 0 0 0 0];
    else
        optval = [0 0 0 0 0];
        for i=1:2:length(options)
            select options(i)
                case 'IntegerTolerance' then
                    optval(1) = options(i+1);
                case 'MaxNodes' then
                    optval(2) = options(i+1);
                case 'MaxTime' then
                    optval(3) = options(i+1);
                case 'AllowableGap' then
                    optval(4) = options(i+1);
                case 'ThreadsNumber' then
                    optval(5) = options(i+1);
                else
                    error(999, 'Unknown string argument passed.');
                end
        end
    end
	
    [xopt,fopt,status,nodes,nfpoints,L,U,niter] = matintlinprog(int32(nbVar),nbCon,c,intcon,conMatrix,conLB,conUB,lb,ub,objSense,optval);

    //Debugging Prints
    //disp(c);disp(intcon);disp(conMatrix);disp(conLB);disp(conUB);disp(lb);disp(ub);disp(Aeq);disp(beq);disp(xopt);
    //disp(L);disp(U);
    //disp(options);
    
    output = struct("relativegap"       , [],..
                    "absolutegap"       , [],..
                    "numnodes"          , [],..
                    "numfeaspoints"     , [],..
                    "numiterations"		, [],..
                    "numberthreads"     , [],..
                    "constrviolation"   , [],..
                    "message"           , '');
                    
    output.numnodes=[nodes];
    output.numfeaspoints=[nfpoints];
    output.numiterations=[niter];
    output.relativegap=(U-L)/(abs(U)+1);
    output.absolutegap=(U-L);
    output.numberthreads=optval(5);
    output.constrviolation = max([0;norm(Aeq*xopt-beq, 'inf');(lb'-xopt);(xopt-ub');(A*xopt-b)]);
    
    select status

    case 0 then
        output.message="Optimal Solution"
    case 1 then 
        output.message="Primal Infeasible"
    case 2 then 
        output.message="Solution Limit is reached"
    case 3 then
        output.message="Node Limit is reached"
    case 4 then 
        output.message="Numerical Difficulties"
    case 5 then 
        output.message="Time Limit Reached"
    case 6 then 
        output.message="Continuous Solution Unbounded"
    case 7 then 
        output.message="Dual Infeasible"
    else
        output.message="Invalid status returned. Notify the Toolbox authors"
        break;
    end

endfunction
