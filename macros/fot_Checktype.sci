// Copyright (C) 2010 - DIGITEO - Michael Baudin
//
// This file must be used under the terms of the GNU LGPL license.

function errmsg = fot_Checktype ( funname , var , varname , ivar , expectedtype )
  // Generates an error if the given variable is not of expected type.
  //
  // Calling Sequence
  //   errmsg = Checktype ( funname , var , varname , ivar , expectedtype )
  //
  // Parameters
  //   funname : a 1 x 1 matrix of strings, the name of the calling function.
  //   var : a 1 x 1 matrix of valid Scilab data type, the variable
  //   varname : a 1 x 1 matrix of string, the name of the variable
  //   ivar : a 1 x 1 matrix of floating point integers, the index of the input argument in the calling sequence
  //   expectedtype : a n x 1 or 1 x n matrix of strings, the available types for the variable #ivar
  //   errmsg : a 1 x 1 matrix of strings, the error message. If there was no error, the error message is the empty matrix.
  //
  // Description
  // This function is designed to be used to design functions with 
  // input arguments with variable type.
  // We use the typeof function to compute the type of the variable:
  // see help typeof to get the list of all available values for expectedtype.
  //   Last update : 29/07/2010.
  //
  // Examples
  // // The function takes a string argument.
  // function myfunction ( x )
  //   Checktype ( "myfunction" , x , "x" , 1 , "string" )
  //   disp("This is a string")
  // endfunction
  // // Calling sequences which work
  // myfunction ( "Scilab" )
  // // Calling sequences which generate an error
  // myfunction ( 123456 )
  //
  // // The function takes a string or a matrix of doubles argument.
  // function myfunction ( x )
  //   Checktype ( "myfunction" , x , "x" , 1 , [ "string" "constant" ] )
  //   if ( typeof(x) == "string" ) then
  //     disp("This is a matrix of strings")
  //   else
  //     disp("This is a matrix of doubles")
  //   end
  // endfunction
  // // Calling sequences which work
  // myfunction ( "Scilab" )
  // myfunction ( 123456 )
  // // Calling sequences which generate an error
  // myfunction ( uint8(2) )
  //
  // Authors
  //   Michael Baudin - 2010 - DIGITEO
  //

  errmsg = []
  if ( and ( typeof ( var ) <> expectedtype ) ) then
    strexp = """" + strcat(expectedtype,""" or """) + """"
    errmsg = msprintf(gettext("%s: Expected type [%s] for input argument %s at input #%d, but got ""%s"" instead."),funname, strexp, varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction



