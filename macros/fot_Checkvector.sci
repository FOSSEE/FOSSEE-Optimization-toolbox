// Copyright (C) 2010 - DIGITEO - Michael Baudin
//
// This file must be used under the terms of the GNU LGPL license.

function errmsg = fot_Checkvector ( funname , var , varname , ivar , nbval )
  // Generates an error if the variable is not a vector.
  //
  // Calling Sequence
  //   errmsg = Checkvector ( funname , var , varname , ivar )
  //
  // Parameters
  //   funname : a 1 x 1 matrix of strings, the name of the calling function.
  //   var : a 1 x 1 matrix of valid Scilab data type, the variable
  //   varname : a 1 x 1 matrix of string, the name of the variable
  //   ivar : a 1 x 1 matrix of floating point integers, the index of the input argument in the calling sequence
  //   nbval : a 1 x 1 matrix of floating point integers, the number of entries in the vector.
  //   errmsg : a 1 x 1 matrix of strings, the error message. If there was no error, the error message is the empty matrix.
  //
  // Description
  // This function is designed to be used to design functions where 
  // the input argument is a vector, that is, a matrix for which 
  // nrows == 1 or ncols == 1.
  // This function cannot be use when var is a function, or more
  // generally, for any input argument for which the size function
  // does not work.
  //
  // Examples
  // // The function takes a vector of 3 doubles.
  // function y = myfunction ( x )
  //   Checkvector ( "myfunction" , x , "x" , 1 , 3 )
  //   y = x
  // endfunction
  // // Calling sequences which work
  // y = myfunction ( ones(1,3) )
  // y = myfunction ( zeros(3,1) )
  // // Calling sequences which generate an error
  // // The following are not vectors
  // y = myfunction ( ones(2,3) )
  // y = myfunction ( zeros(3,2) )
  // // The following have the wrong number of entries
  // y = myfunction ( ones(1,3) )
  //
  // Authors
  //   Michael Baudin - 2010 - DIGITEO
  //

  errmsg = []
  nrows = size(var,"r")
  ncols = size(var,"c")
  if ( nrows <> 1 & ncols <> 1 ) then
    strcomp = strcat(string(size(var))," ")
    errmsg = msprintf(gettext("%s: Expected a vector matrix for input argument %s at input #%d, but got [%s] instead."), funname, varname , ivar , strcomp );
    error(errmsg)
  end
  if ( ( nrows == 1 & ncols <> nbval ) | ( ncols == 1 & nrows <> nbval ) ) then
    strcomp = strcat(string(size(var))," ")
    errmsg = msprintf(gettext("%s: Expected %d entries for input argument %s at input #%d, but current dimensions are [%s] instead."), funname, nbval , varname , ivar , strcomp );
    error(errmsg)
  end
endfunction



