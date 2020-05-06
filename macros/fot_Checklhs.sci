// Copyright (C) 2010 - DIGITEO - Michael Baudin
//
// This file must be used under the terms of the GNU LGPL license.

function errmsg = fot_Checklhs ( funname , lhs , lhsset )
  // Generates an error if the number of LHS is not in given set.
  //
  // Calling Sequence
  //   errmsg = Checklhs ( funname , lhs , lhsset )
  //
  // Parameters
  //   funname : a 1 x 1 matrix of strings, the name of the calling function.
  //   lhs : a 1 x 1 matrix of floating point integers, the actual number of output arguments
  //   lhsset : a 1 x n or n x 1 matrix of floating point integers, the authorized number of output arguments
  //   errmsg : a 1 x 1 matrix of strings, the error message. If there was no error, the error message is the empty matrix.
  //
  // Description
  // This function is designed to be used to design functions with 
  // variable number of output arguments.
  // Notice that it is useless to call this function if the 
  // function definition does not use the varargout statement.
  // Notice that a function as a minimum of 1 output argument.
  //   Last update : 29/07/2010.
  //
  // Examples
  // // The function takes 3 input arguments and 1/2 output arguments
  // function varargout = myfunction ( x1 , x2 , x3 )
  //   [lhs, rhs] = argn()
  //   Checkrhs ( "myfunction" , rhs , 3 : 3 )
  //   Checklhs ( "myfunction" , lhs , 1 : 2 )
  //   y1 = x1 + x2
  //   y2 = x2 + x3
  //   varargout(1) = y1
  //   if ( lhs == 2 ) then
  //     varargout(2) = y2
  //   end
  // endfunction
  // // Calling sequences which work
  // myfunction ( 1 , 2 , 3 )
  // y1 = myfunction ( 1 , 2 , 3 )
  // [ y1 , y2 ] = myfunction ( 1 , 2 , 3 )
  // // Calling sequences which generate an error
  // [ y1 , y2 , y3 ] = myfunction ( 1 , 2 , 3 )
  //
  // // The function takes 1 or 3 output arguments, but not 2
  // function varargout = myfunction ( x1 , x2 , x3 )
  //   [lhs, rhs] = argn()
  //   Checkrhs ( "myfunction" , rhs , 3 : 3 )
  //   Checklhs ( "myfunction" , lhs , [1 3] )
  //   y1 = x1 + x2
  //   y2 = x2 + x3
  //   y3 = x1 + x3
  //   varargout(1) = y1
  //   if ( lhs == 3 ) then
  //     varargout(2) = y2
  //     varargout(3) = y3
  //   end
  // endfunction
  // // Calling sequences which work
  // myfunction ( 1 , 2 , 3 )
  // y1 = myfunction ( 1 , 2 , 3 )
  // [ y1 , y2 , y3 ] = myfunction ( 1 , 2 , 3 )
  // // Calling sequences which generate an error
  // [y1 , y2] = myfunction ( 1 , 2 , 3 )
  //
  // Authors
  //   Michael Baudin - 2010 - DIGITEO
  //
  
  errmsg = []
  if ( and ( lhs <> lhsset ) ) then
    lhsstr = strcat(string(lhsset)," ")
    errmsg = msprintf(gettext("%s: Unexpected number of output arguments : %d provided while the expected number of output arguments should be in the set [%s]."), funname , lhs , lhsstr );
    error(errmsg)
  end
endfunction



