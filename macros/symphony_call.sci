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

function [xopt,fopt,status,output] = symphony_call(nbVar,nbCon,objCoef,isInt,lb,ub,A,conLB,conUB,objSense,options)

    xopt = [];
    fopt = [];
    status = [];
    output = [];
    
    //Opening Symphony environment 
    sym_open();

    //Setting Options for the Symphpony
    setOptions(options);
    
   //Choosing to launch basic or advanced version
    if(~issparse(A)) then
        sym_loadProblemBasic(int32(nbVar), int32(nbCon),lb,ub,objCoef,isInt,objSense,A,conLB,conUB);
    else
        // Changing to Constraint Matrix into sparse matrix
        A_advanced=sparse(A);
        sym_loadProblem(nbVar,nbCon,lb,ub,objCoef,isInt,objSense,A_advanced,conLB,conUB);
    end

    op = sym_solve();
    
    status = sym_getStatus();
     if (status == 228 | status == 227 | status == 229 | status == 230 | status == 231 | status == 232 | status == 233) then
            xopt = sym_getVarSoln();
            // Symphony gives a row matrix converting it to column matrix
            xopt = xopt';
            
            fopt = sym_getObjVal();
    end
    status = sym_getStatus();
	output = struct("Iterations", []);
    output.Iterations = sym_getIterCount();


endfunction
