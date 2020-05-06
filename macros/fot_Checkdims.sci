
// Copyright (C) 2010 - DIGITEO - Michael Baudin
//
// This file must be used under the terms of the GNU LGPL license.

function errmsg = fot_Checkdims ( funname , var , varname , ivar , matdims )
  // Generates an error if the variable has not the required size.
  //
  // Calling Sequence
  //   errmsg = Checkdims ( funname , var , varname , ivar , matdims )
  //
  // Parameters
  //   funname : a 1 x 1 matrix of strings, the name of the calling function.
  //   var : a 1 x 1 matrix of valid Scilab data type, the variable
  //   varname : a 1 x 1 matrix of string, the name of the variable
  //   ivar : a 1 x 1 matrix of floating point integers, the index of the input argument in the calling sequence
  //   matdims : 1 x 2 matrix of floating point integers, the number of rows, columns for the variable #ivar
  //   errmsg : a 1 x 1 matrix of strings, the error message. If there was no error, the error message is the empty matrix.
  //
  // Description
  // This function is designed to be used to design functions where 
  // the input argument has a known shape.
  // This function cannot be use when var is a function, or more
  // generally, for any input argument for which the size function
  // does not work.
  //   Last update : 05/08/2010.
  //
  // Examples
  // // The function takes a 2 x 3 matrix of doubles.
  // function y = myfunction ( x )
  //   Checkdims ( "myfunction" , x , "x" , 1 , [2 3] )
  //   y = x
  // endfunction
  // // Calling sequences which work
  // y = myfunction ( ones(2,3) )
  // y = myfunction ( zeros(2,3) )
  // // Calling sequences which generate an error
  // y = myfunction ( ones(1,3) )
  // y = myfunction ( zeros(2,4) )
  //
  // Authors
  //   Michael Baudin - 2010 - DIGITEO
  //

  [lhs,rhs]=argn()
  fot_Checkrhs ( funname , rhs , 5 )
  fot_Checklhs ( funname , lhs , [0 1] )

  errmsg = []
  if ( or ( size(var) <> matdims ) ) then
    strexp = strcat(string(matdims)," ")
    strcomp = strcat(string(size(var))," ")
    errmsg = msprintf(gettext("%s: Expected size [%s] for input argument %s at input #%d, but got [%s] instead."), funname, strexp, varname , ivar , strcomp );
    error(errmsg)
  end
endfunction
