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

	//To check the number of input and output argument
   	[lhs , rhs] = argn();
	
	//To check the number of argument given by user
   	if ( rhs < 1 | rhs > 2) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [1 2]"),"linprog",rhs);
		error(errmsg)
   	end
   	mpsFile = varargin(1);

	[sizeFile isFile] = fileinfo(mpsFile);
	if(isFile ~= 0) then
		errmsg = msprintf(gettext("%s: File is not present at the given location"),"linprog",rhs);
		error(errmsg)
	end

   if ( rhs<2 | size(varargin(2)) ==0 ) then
      param = list();
   else
       param =varargin(2);
   end 

   if (type(param) ~= 15) then
      errmsg = msprintf(gettext("%s: options should be a list "), "linprog");
      error(errmsg);
   end

	//Check type of variables
	Checktype("linprog", mpsFile, "mpsFile", 1, "string")

   if (modulo(size(param),2)) then
   errmsg = msprintf(gettext("%s: Size of parameters should be even"), "linprog");
   error(errmsg);
   end
   options = list("MaxIter"     , [3000],);

   for i = 1:(size(param))/2
        select convstr(param(2*i-1),'l')
            case "maxiter" then
        		options(2*i) = param(2*i);
        	else
			  errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "linprog", param(2*i-1));
			  error(errmsg)
		end		
   end
   //Calling the function by passing the required parameters 
   
   [xopt,fopt,status,iter,Zl,dual] = rmps(mpsFile,options);
   
   xopt = xopt';
   exitflag = status;
   output = struct("Iterations"      , [],..
                   "constrviolation"	, []);
                  
   output.Iterations = iter;
   output.constrviolation =[] //max([0;norm(Aeq*xopt'-beq, 'inf');(lb-xopt');(xopt'-ub);(A*xopt'-b)]);
   lambda = struct("lower"           , [], ..
                   "upper"           , [], ..
                   "dual"           , []..
                  );
   
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
	
	
	lambda.dual =dual;
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
