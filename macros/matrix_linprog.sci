// Copyright (C) 2015 - IIT Bombay - FOSSEE
// Author: Guru Pradeep Reddy, Bhanu Priya Sayal
// Organization: FOSSEE, IIT Bombay
// Email: gurupradeept@gmail.com, bhanupriyasayal@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt


function [xopt,fopt,exitflag,output,lambda] = matrix_linprog (varargin)
//To check the number of input and output argument
   [lhs , rhs] = argn();
	
	c = [];
	A = [];
	b = [];
	Aeq = [];
	beq = [];
	lb = [];
	ub = [];
	options = list();

//To check the number of argument given by user
   if ( rhs < 3 | rhs == 4 | rhs == 6 | rhs >8 ) then
    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [3 5 7 8]"), "linprog", rhs);
    error(errmsg)
   end
   
   c = varargin(1);

	if(size(c,2) == 0) then
		errmsg = msprintf(gettext("%s: Cannot determine the number of variables because input objective coefficients is empty"), "linprog");
		error(errmsg);
	end

   if (size(c,2)~=1) then
	errmsg = msprintf(gettext("%s: Objective Coefficients should be a column matrix"), "linprog");
	error(errmsg);
   end

   nbVar = size(c,1);
   A = varargin(2);
   b = varargin(3);

	if ( rhs<4 ) then
	  	Aeq = []
	  	beq = []
	else
	  	Aeq = varargin(4);
	  	beq = varargin(5);
	end

	if ( rhs<6 ) then
		lb = [];
		ub = [];
	else
		lb = varargin(6);
		ub = varargin(7);
	end

   	if ( rhs<8 | size(varargin(8)) ==0 ) then
    	param = list();
   	else
	  	param =varargin(8);
   	end

	//Check type of variables
	Checktype("linprog", c, "c", 1, "constant")
	Checktype("linprog", A, "A", 2, "constant")
	Checktype("linprog", b, "b", 3, "constant")
	Checktype("linprog", Aeq, "Aeq", 4, "constant")
	Checktype("linprog", beq, "beq", 5, "constant")
	Checktype("linprog", lb, "lb", 6, "constant")
	Checktype("linprog", ub, "ub", 7, "constant")

   nbConInEq = size(A,1);
   nbConEq = size(Aeq,1);

    lb = lb(:);
    ub = ub(:)
    b = b(:);
    beq = beq(:);

   	if (size(lb,2)==0) then
        lb = repmat(-%inf,nbVar,1);
    end
    
    if (size(ub,2)==0) then
        ub = repmat(%inf,nbVar,1);
    end

    if (type(param) ~= 15) then
      errmsg = msprintf(gettext("%s: options should be a list "), "linprog");
      error(errmsg);
    end
   

   if (modulo(size(param),2)) then
   errmsg = msprintf(gettext("%s: Size of parameters should be even"), "linprog");
   error(errmsg);
   end

   options = list("MaxIter", [3000]);

   for i = 1:(size(param))/2
        
      	select convstr(param(2*i-1),'l')
    	case "maxiter" then
       		options(2*i) = param(2*i);
		else
			  errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "linprog", param(2*i-1));
			  error(errmsg)
		end
     end

   //Check the size of inequality constraint which should be equal to the number of variables
   if ( size(A,2) ~= nbVar & size(A,2) ~= 0) then
      errmsg = msprintf(gettext("%s: The number of columns in A must be the same as the number of elements of c"), "linprog");
      error(errmsg);
   end

   //Check the size of equality constraint which should be equal to the number of variables
   if ( size(Aeq,2) ~= nbVar & size(Aeq,2) ~= 0 ) then
      errmsg = msprintf(gettext("%s: The number of columns in Aeq must be the same as the number of elements of c"), "linprog");
      error(errmsg);
   end
	
   //Check the size of Lower Bound which should be equal to the number of variables
   if ( size(lb,1) ~= nbVar) then
      errmsg = msprintf(gettext("%s: The Lower Bound is not equal to the number of variables"), "linprog");
      error(errmsg);
   end

   //Check the size of Upper Bound which should equal to the number of variables
   if ( size(ub,1) ~= nbVar) then
      errmsg = msprintf(gettext("%s: The Upper Bound is not equal to the number of variables"), "linprog");
      error(errmsg);
   end
   //Check the size of constraints of Lower Bound which should equal to the number of constraints
   if ( size(b,1) ~= nbConInEq & size(b,2) ~= 0) then
      errmsg = msprintf(gettext("%s: The number of rows in A must be the same as the number of elements of b"), "linprog");
      error(errmsg);
   end

   //Check the size of constraints of Upper Bound which should equal to the number of constraints
   if ( size(beq,1) ~= nbConEq & size(beq,2) ~= 0) then
      errmsg = msprintf(gettext("%s: The number of rows in Aeq must be the same as the number of elements of beq"), "linprog");
      error(errmsg);
   end

   //Check if the user gives a matrix instead of a vector
   if (size(lb,1)~=1)& (size(lb,2)~=1) then
      errmsg = msprintf(gettext("%s: Lower Bound should be a vector"), "linprog");
      error(errmsg); 
   end
   
   if (size(ub,1)~=1)& (size(ub,2)~=1) then
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

	nbVar = int32(size(c,1));
	c = c';
	conMatrix = [Aeq;A];
	nbCon = int32(size(conMatrix,1));
	conlb = [beq; repmat(-%inf,nbConInEq,1)];
	conub = [beq;b];
	lb = lb';
	ub = ub';

   [xopt,fopt,status,iter,Zl,dual] = linearprog(nbVar,nbCon,c,conMatrix,conlb,conub,lb,  ub,options);
   
   xopt = xopt';
   exitflag = status;
   output = struct("Iterations"      , [],..
                   "constrviolation"	, []);
                  
   output.Iterations = iter;
   output.constrviolation = []//max([0;norm(Aeq*xopt-beq, 'inf');(lb-xopt);(xopt-ub);(A*xopt-b)]);
   lambda = struct("lower"           , [], ..
                   "upper"           , [], ..
                   "ineqlin"           , [], ..
                   "eqlin"      , []);
   
    rc = Zl; //The reduced cost vector
    for i= 1:length(rc)
    
        if abs(rc(i)) < 10^-6 then
            rc(i) = 0;
        end
        
        if rc(i) == 0 then
            lambda.lower(i) = 0;
            lambda.upper(i) = 0;
        elseif rc(i) > 0 then
            lambda.lower(i) = rc(i);
            lambda.upper(i) = 0;
        else
            lambda.lower(i) = 0;
            lambda.upper(i) = rc(i);
        end
             
    end
	lambda.eqlin =dual(1:nbConEq);
	lambda.ineqlin =dual(nbConEq+1:nbCon);
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
	else
    	printf("\nInvalid status returned. Notify the Toolbox authors\n");
        break;
    end

endfunction
