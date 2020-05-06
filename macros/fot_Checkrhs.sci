// Copyright (C) 2010 - DIGITEO - Michael Baudin
//
// This file must be used under the terms of the GNU LGPL license.

function errmsg = fot_Checkrhs ( funname , rhs , rhsset )
  // Generates an error if the number of RHS is not in given set.
  //
  // Calling Sequence
  //   errmsg = Checkrhs ( funname , rhs , rhsset )
  //
  // Parameters
  //   funname : a 1 x 1 matrix of strings, the name of the calling function.
  //   rhs : a 1 x 1 matrix of floating point integers, the actual number of input arguments
  //   rhsset : a 1 x n or n x 1 matrix of floating point integers, the authorized number of input arguments
  //   errmsg : a 1 x 1 matrix of strings, the error message. If there was no error, the error message is the empty matrix.
  //
  // Description
  // This function is designed to be used to design functions with 
  // variable number of input arguments.
  // Notice that it is useless to call this function if the 
  // function definition does not use the varargin statement.
  //   Last update : 05/08/2010.
  //   Last update : 29/07/2010.
  //
  // Examples
  // // The function takes 2/3 input arguments and 1 output arguments
  // function y = myfunction ( varargin )
  //   [lhs, rhs] = argn()
  //   Checkrhs ( "myfunction" , rhs , 2:3 )
  //   Checklhs ( "myfunction" , lhs , 1 )
  //   x1 = varargin(1)
  //   x2 = varargin(2)
  //   if ( rhs >= 3 ) then
  //     x3 = varargin(3)
  //   else
  //     x3 = 2
  //   end
  //   y = x1 + x2 + x3
  // endfunction
  // // Calling sequences which work
  // y = myfunction ( 1 , 2 )
  // y = myfunction ( 1 , 2 , 3 )
  // // Calling sequences which generate an error
  // y = myfunction ( 1 )
  // y = myfunction ( 1 , 2 , 3 , 4 )
  //
  // // The function takes 2 or 4 input arguments, but not 3
  // function y = myfunction ( varargin )
  //   [lhs, rhs] = argn()
  //   Checkrhs ( "myfunction" , rhs , [2 4] )
  //   Checklhs ( "myfunction" , lhs , 1 )
  //   x1 = varargin(1)
  //   x2 = varargin(2)
  //   if ( rhs >= 3 ) then
  //     x3 = varargin(3)
  //     x4 = varargin(4)
  //   else
  //     x3 = 2
  //     x4 = 3
  //   end
  //   y = x1 + x2 + x3 + x4
  // endfunction
  // // Calling sequences which work
  // y = myfunction ( 1 , 2 )
  // y = myfunction ( 1 , 2 , 3 , 4 )
  // // Calling sequences which generate an error
  // y = myfunction ( 1 )
  // y = myfunction ( 1 , 2 , 3 )
  // y = myfunction ( 1 , 2 , 3 , 4, 5 )
  //
  // // The function takes 2 input arguments and 0/1 output arguments.
  // // Notice that if the checkrhs function is not called,
  // // the variable x2 might be used from the user's context,
  // // that is, if the caller has defined the variable x2, it 
  // // is used in the callee.
  // // Here, we want to avoid this.
  // function y = myfunction ( x1 , x2 )
  //   [lhs, rhs] = argn()
  //   Checkrhs ( "myfunction" , rhs , 2 )
  //   Checklhs ( "myfunction" , lhs , [0 1] )
  //   y = x1 + x2
  // endfunction
  // // Calling sequences which work
  // y = myfunction ( 1 , 2 )
  // // Calling sequences which generate an error
  // y = myfunction ( 1 )
  // y = myfunction ( 1 , 2 , 3 )
  //
  // Authors
  //   Michael Baudin - 2010 - DIGITEO
  //

  errmsg = []
  if ( and(rhs <> rhsset) ) then
    rhsstr = strcat(string(rhsset)," ")
    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while the number of expected input arguments should be in the set [%s]."), funname , rhs , rhsstr );
    error(errmsg)
  end
endfunction



