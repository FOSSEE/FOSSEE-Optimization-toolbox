// to check if the number of input arguments is equal to the numbers mentioned in the set, which is either 2,4,6,8,9 or 10
function errmsg = fgoalattain_checkrhs ( funname , rhs , rhsset )  
 errmsg = []
  if ( and(rhs <> rhsset) ) then
    rhsstr = strcat(string(rhsset)," ")
    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while the number of expected input arguments should be in the set [%s]."), funname , rhs , rhsstr );
    error(errmsg)
  end
endfunction

// to check if the number of output arguments is equal to the numbers mentioned in the set, which is any number between 1 to 6 (both inclusive)
function errmsg = fgoalattain_checklhs ( funname , lhs , lhsset )
errmsg = []
  if ( and ( lhs <> lhsset ) ) then
    lhsstr = strcat(string(lhsset)," ")
    errmsg = msprintf(gettext("%s: Unexpected number of output arguments : %d provided while the expected number of output arguments should be in the set [%s]."), funname , lhs , lhsstr );
    error(errmsg)
  end
endfunction

// fmincon library of scilab gives error when 'c' component of nonlinfun empty
// add a trivial case of -5 <= 0 to c to bypass this error
function [c,ceq] = nonlinfun(x)
    c = [-5; objfun(x)-weight(x)*x0(numvar+1)-goal(x)]
    
    ceq = []
  endfunction

// Generates an error if the given variable is not of expected type.
function errmsg = fgoalattain_checktype ( funname , var , varname , ivar , expectedtype )
errmsg = []
  if ( and ( typeof ( var ) <> expectedtype ) ) then
    strexp = """" + strcat(expectedtype,""" or """) + """"
    errmsg = msprintf(gettext("%s: Expected type [%s] for input argument %s at input #%d, but got ""%s"" instead."),funname, strexp, varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction

// Generates an error if the variable is not a vector.
function errmsg = fgoalattain_checkvector ( funname , var , varname , ivar , nbval )
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

// Generates an error if the variable has not the required size.
function errmsg = fgoalattain_checkdims ( funname , var , varname , ivar , matdims )
[lhs,rhs]=argn()
  fgoalattain_checkrhs ( "fgoalattain_checkdims" , rhs , 5 )
  fgoalattain_checklhs ( "fgoalattain_checkdims" , lhs , [0 1] )
  errmsg = []
  if ( or ( size(var) <> matdims ) ) then
    strexp = strcat(string(matdims)," ")
    strcomp = strcat(string(size(var))," ")
    errmsg = msprintf(gettext("%s: Expected size [%s] for input argument %s at input #%d, but got [%s] instead."), funname, strexp, varname , ivar , strcomp );
    error(errmsg)
  end
endfunction


