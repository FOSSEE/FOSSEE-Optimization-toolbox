//// Copyright (C) 2015 - IIT Bombay - FOSSEE
////
//// This file must be used under the terms of the CeCILL.
//// This source file is licensed as described in the file COPYING, which
//// you should have received as part of this distribution.  The terms
//// are also available at
//// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
//// Author: Harpreet Singh
//// Organization: FOSSEE, IIT Bombay
//// Email: toolbox@scilab.in
//
//function setOptions(varargin)
//
//    options = varargin(1);
//    nbOpt = size(options);
//    
//
//    if (nbOpt~=0) then
//        for i = 1:(nbOpt/2)
//        
//        //Setting the parameters
//          
//            //Check if the given parameter is String
//            if (type(options(2*i)) == 10 ) then
//                sym_setStrParam(options(2*i - 1),options(2*i));
//          
//            //Check if the given parameter is Double
//            elseif(type(options(2*i))==1) then
//                sym_setDblParam(options(2*i - 1),options(2*i));
//             
//            //Check if the given parameter is Integer
//            elseif(type(options(2*i))==8)
//                sym_setIntParam(options(2*i - 1),options(2*i)); 
//            end
//         end
//    end
//    
//endfunction
//
