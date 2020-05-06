function [x,y,s,z,info,status] = ecos(varargin)
    //   Solves conic optimization problems with integer and boolean constraints.
    //
    //   Calling Sequence
    //   x = ecos(c,G,h,dims)
    //   x = ecos(c,G,h,dims,A,b)
    //   x = ecos(c,G,h,dims,A,b,options)
    //   [x,y,s,z,info,status] = ecos( ... )
    //   
    //   Parameters
    //   c : a vector of double, contains coefficients of the variables in the objective.
    //   G : a matrix of double, represents the linear coefficients in the inequality constraints G⋅x ≤ h with respect to cone K. 
    //   h : a vector of double, represents the linear coefficients in the inequality constraints G⋅x ≤ h with respect to cone K.
    //  dims : a list containing the number of n-dimensional positive orthant(l),the second order cones(q) and exponential cones(e).
    //   A : a matrix of double, represents the linear coefficients in the equality constraints A⋅x = b.
    //   b : a vector of double, represents the linear coefficients in the equality constraints A⋅x = b.
    //   options : a list containing the parameters to be set.
    //   x : a vector of double, the primal solution variable of the optimization problem.
    //   y : a vector of double,  the dual variables fo equality constraints of the optimization problem.
    //   s : a vector of double,  the slack for inequality constraints of the optimization problem.
    //   z : a vector of double,  the dula varibale for inequality constraints of the optimization problem.
    //   info : a struct defining different parameters of the ecos solver. See below for details.
    //   status : The status returns the solver status after the optimization process. See below for details.
    //   
    //   Description
    //   Search the minimum of a conic constrained mixed integer programming optimization problem specified by :
    //
    //   <latex>
    //  \begin{eqnarray}
    // &\mbox{min}_{x}
    // & C^T⋅x \\
    // & \text{subject to}& A⋅x = b \\
    // & & G⋅x \preceq_K h \\
    // & & x_i \in \{0,1\}, i \in \!\, bool_vars\
    // & & x_j \in \!\, \mathbb{Z}, i \in \!\, int_vars\\
    //  \end{eqnarray}
    //   </latex>
    //
    //    Examples 
    // // Objective function
    // c = [-750 -1000];
    // // Constraint inequality matrix LHS
    // G = [
    //   1 1;
    //   1 2;
    //   4 3;
    //     ];
    // // Constraint inequality matrix RHS
    // h = [10 15 25]';
    // // Constraint equality matrix LHS
    // A = [
    //   0.5 1
    // ];
    // // Constraint equality matrix RHS
    // b=[7.5];
    // // Dimension of positive orthant
    // l = [3];
    // q = [];
    // e = [];
    // dims=list("l",l,"q",q,"e",e)
    //  //Calling ecos
    // [x,y,s,z,info,status] =ecos(c,G,h,dims,A,b);

    // Examples
    //
    // // Objective function
    // c = [0 0 0 0 1];
    //  //Constraint inequality matrix LHS
    // G = [
    //   0.4167578    0.0562668    0.           0.           0.  
    //     2.1361961   -1.6402708    0.           0.           0.  
    //     1.7934356    0.8417474    0.           0.           0.  
    //     0.           0.           0.4167578    0.0562668    0.  
    //     0.           0.           2.1361961   -1.6402708    0.  
    //     0.           0.           1.7934356    0.8417474    0.  
    //     0.           0.           0.           0.          -1.  
    //    -1.           0.           0.           0.           0.  
    //     0.          -1.           0.           0.           0.  
    //     0.           0.          -1.           0.           0.  
    //     0.           0.           0.          -1.           0.
    //     ];
    //  //Constraint inequality matrix RHS
    // h = [0 0  0  0  0  0  0  0  0  0  0]';
    // // Dimension of positive orthant
    // l = [6];
    // q = [5];
    // e = [0]
    // dims=list("l",l,"q",q,"e",e)
    // [x,y,s,z,info,status] =ecos(c,G,h,dims);
  // Author
  // Georgey John

  function [A1,b1,s0]=linconcheck(A,b,inputs_name)
    //Function to check the linear inputs A,b and Aeq,beq

    if(size(b,2)>1) then
      errmsg = msprintf(gettext("%s: Expected Column vector of size (Number of constraints) for %s"), "ecos",inputs_name(2));
      error(errmsg);
    end

    s0=size(A);
    //To check for correct size of A
    if(s0(2)==0) then 
      if(size(b,2)~=0) then
        errmsg = msprintf(gettext("%s: As Linear Inequality Constraint coefficient Matrix %s is empty, %s should also be empty"), "ecos",inputs_name(1),inputs_name(2));
        error(errmsg);
      end
    else
      if((size(b,1)~=1) & (size(b,2)~=1)) then
        errmsg = msprintf(gettext("%s: Expected Non empty Row/Column Vector for %s for your Inputs "), "ecos",inputs_name(2));
        error(errmsg);
      elseif(size(b,1)~=s0(1) & size(b,2)==1) then
        errmsg = msprintf(gettext("%s: Expected Column Vector (number of linear inequality constraints X 1) for %s"), "ecos",inputs_name(2));
        error(errmsg);
      // elseif(size(b,1)==1 & size(b,2)~=s0(1)) then
      //   errmsg = msprintf(gettext("%s: Expected Row Vector (1 X number of linear inequality constraints) for %s"), "ecos",inputs_name(2));
      //   error(errmsg);
      end 
    end
    b1=b;

    //To check for corrcet size of A
    if(size(A,1)~=size(b,1) & size(A,2)~=0) then
      errmsg = msprintf(gettext("%s: Expected Matrix of size (No of linear inequality constraints X No of Variables) or an Empty Matrix for Linear Inequality Constraint coefficient Matrix %s"), "ecos",inputs_name(1));
      error(errmsg);
    end
    A1=A;
    
  endfunction

	[lhs ,rhs]=argn() 
	
	if ( rhs<4 | rhs>9 ) then
	    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while it should be 4,6,7,8 or 9"), "ecos", rhs);
	    error(errmsg);
  	end

  	if (rhs==5) then
	   	errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while it should be 4,6,8,10 or 11"), "ecos", rhs);
	   	error(errmsg);
 	end

 	c = varargin(1);
 	G = varargin(2);
 	h =	varargin(3);
 	dims = varargin(4);
 	
  A=[];Apr=[];Ajc=[];Air=[];
 	b=[];
 	param=list();
 	m=0,n=0,p=0;

 	if (rhs>4 & rhs<7) then
	 	A = varargin(5);
	 	b = varargin(6);
	else
		A=[];
		b=[];
	end

	if (rhs ==7) then
 		param = varargin(7);
 	end

 	fot_Checktype("ecos", c, "c", 1, "constant");

  // Error check for objective matrix(c)
  if (size(c,1)~=1) then
    errmsg = msprintf(gettext("%s: c should be a row vector"), "ecos");
    error(errmsg);
  elseif (size(c,1)==0) then
    errmsg = msprintf(gettext("%s: c cannot be empty"), "ecos");
    error(errmsg);
  end

 	s=size(c,2);


 	fot_Checktype("ecos", G, "G", 2, "constant");
 	fot_Checktype("ecos", h, "h", 3, "constant");

  [m,n]=size(G);

  // Error check and converting inequaity matrix(G) to CCS format required by ecos
  if (n~=size(c,2)) then
    errmsg = msprintf(gettext("%s: Number of columns of G and c do not match"), "ecos");
    error(errmsg);
  end

  [G,h,s1]=linconcheck(G,h,["G","h"]);
  [Gjc,Gir,Gpr] = sp2adj(sparse(G));

  // Extracting values and Error checks for dims
 	fot_Checktype("ecos", dims, "dims", 4, "list");

 	if (isempty(dims) | modulo(size(dims),2)) then
 		errmsg = msprintf(gettext("%s: dims cannot be empty and should be list of even size"), "ecos");
   			error(errmsg);
	end

	for i = 1:(size(dims))/2
    select convstr(dims(2*i-1),'l')
      case "l" then
        fot_Checktype("ecos", dims(2*i), "l", 4, "constant");
        if (isempty(dims(2*i))) then
        	l=0;
        elseif (dims(2*i)<0 | modulo(dims(2*i),1)) then
        	errmsg = msprintf(gettext("%s: l in dims should be a positive integer"), "ecos");
          	error(errmsg);
        else
        	l=dims(2*i);
          // if (prod(size(l))>1)) then
          //   errmsg = msprintf(gettext("%s: l in dims should be a positive integer"), "ecos");
          //   error(errmsg);
          // end
        end
      case "q" then
        fot_Checktype("ecos", dims(2*i), "q", 4, "constant");
        if (isempty(dims(2*i))) then
        	q=[];
        elseif (dims(2*i)<0 | modulo(dims(2*i),1)) then
        	errmsg = msprintf(gettext("%s: q in dims should be a positive vector"), "ecos");
          	error(errmsg);
        else
        	q=dims(2*i);
        end
       case "e" then
        fot_Checktype("ecos", dims(2*i), "e", 4, "constant");
        if (isempty(dims(2*i))) then
        	e=0;
        elseif (dims(2*i)<0 | modulo(dims(2*i),1)) then
        	errmsg = msprintf(gettext("%s: e in dims should be a positive integer"), "ecos");
          	error(errmsg);
        else
        	e=dims(2*i);
          // if (prod(size(e))>1)) then
          //   errmsg = msprintf(gettext("%s: e in dims should be a positive integer"), "ecos");
          //   error(errmsg);
          // end
        end
    else
      errmsg = msprintf(gettext("%s: Unrecognized parameter name %s."), "ecos", dims(2*i-1));
      error(errmsg);
    end                  
  end

  // Error check and converting equaity matrix(A) to CCS format required by ecos
  if ((size(A,1)*size(b,1)==0) & (size(A,1)+size(b,1)~=0)) then
    errmsg = msprintf(gettext("%s: One of %s and %s is an empty matrix"), "ecos","A","b");
    error(errmsg);
  end

	if (size(A,1)~=0) then
		fot_Checktype("ecos", A, "A", 5, "constant");
    fot_Checktype("ecos", b, "b", 6, "constant");

    [A,b,s2]=linconcheck(A,b,["A","b"]);
    [Ajc,Air,Apr] = sp2adj(sparse(A));
    [m2,p]= size(A);

		if (n~=p) then
			errmsg = msprintf(gettext("%s: Number of columns of G and A do not match"), "ecos");
      error(errmsg);
    end

    if (p~=size(c,2)) then
      errmsg = msprintf(gettext("%s: Number of columns of A and c do not match"), "ecos");
      error(errmsg);
    end
	end

  if (size(G,1)==0 & (size(h,1)==0)) then
    if (size(A,1)==0 & (size(b,1)==0)) then
      errmsg = msprintf(gettext("%s: At most one of the pair (G, h) or (A, b) is allowed to be absent"), "ecos");
      error(errmsg);
    end
  end

  //  Extracting values and Error checks for options(param)
 	fot_Checktype("ecos", param, "param", 1, "list");
 	if (modulo(size(param),2)) then
    	errmsg = msprintf(gettext("%s: Size of Options (list) should be even"), "ecos");
    	error(errmsg);
  	end

  	option = list("maxiter", [100], "feastol", [1e-8],"reltol",[1e-8],"abstol",[1e-8],"feastol_inacc",[1e-4],"abstol_inacc",[5e-5],"reltol_inacc",[5e-5],"verbose",[0],"mi_max_iters",[1000],"mi_int_tol",[1e-4],"mi_abs_eps",[1e-6],"mi_rel_eps",[1e-6]);

 	for i = 1:(size(param))/2
	    select convstr(param(2*i-1),'l')
	      case "maxit" then
	        if (type(option(2*i))~=1 | modulo(option(2*i),1)) then
	          errmsg = msprintf(gettext("%s: Value for Maximum Iteration should be a Constant integer"), "ecos");
	          error(errmsg);
	        else
	          option(2) = param(2*i);    
	        end

	      case "feastol" then
	        if (type(option(2*i))~=1) then
	          errmsg = msprintf(gettext("%s: Value for tolerance should be a Constant"), "ecos");
	          error(errmsg);
	        else
	          option(4) = param(2*i);    
	        end

	      case "reltol" then
	        if (type(option(2*i))~=1) then
	          errmsg = msprintf(gettext("%s: Value for relative tolerance should be a Constant"), "ecos");
	          error(errmsg);
	        else
	          option(6) = param(2*i);    
	        end

	      case "abstol" then
	        if (type(option(2*i))~=1) then
	          errmsg = msprintf(gettext("%s: Value for absolute tolerance should be a Constant"), "ecos");
	          error(errmsg);
	        else
	          option(8) = param(2*i);    
	        end

	      case "feastol_inacc" then
	        if (type(option(2*i))~=1) then
	          errmsg = msprintf(gettext("%s: Value for tolerance with reduced precision should be a Constant"), "ecos");
	          error(errmsg);
	        else
	          option(10) = param(2*i);    
	        end

	      case "reltol_inacc" then
	        if (type(option(2*i))~=1) then
	          errmsg = msprintf(gettext("%s: Value for relative tolerance with reduced precision should be a Constant"), "ecos");
	          error(errmsg);
	        else
	          option(12) = param(2*i);
	        end

	      case "abstol_inacc" then
	        if (type(option(2*i))~=1) then
	          errmsg = msprintf(gettext("%s: Value for absolute tolerance with reduced precision should be a Constant"), "ecos");
	          error(errmsg);
	        else
	          option(14) = param(2*i);
	        end

	      case "verbose" then
	        if (type(option(2*i))~=1 | modulo(option(2*i),1)) then
	          errmsg = msprintf(gettext("%s: Value for verbose level should be a Constant integer"), "ecos");
	          error(errmsg);
	        else
	          option(16) = param(2*i);
	        end

	      case "mi_max_iters" then
	        if (type(option(2*i))~=1 | modulo(option(2*i),1)) then
	          errmsg = msprintf(gettext("%s: Value for maximum branch and bound iterations should be a Constant integer"), "ecos");
	          error(errmsg);
	        else
	          option(18) = param(2*i);
	        end

	      case "mi_int_tol" then
	        if (type(option(2*i))~=1) then
	          errmsg = msprintf(gettext("%s: Value for integer tolerance should be a Constant"), "ecos");
	          error(errmsg);
	        else
	          option(20) = param(2*i);
	        end

	      case "mi_abs_eps" then
	        if (type(option(2*i))~=1) then
	          errmsg = msprintf(gettext("%s: Value for absolute tolerance between upper and lower bounds should be a Constant"), "ecos");
	          error(errmsg);
	        else
	          option(22) = param(2*i);
	        end

	      case "mi_rel_eps" then
	        if (type(option(2*i))~=1) then
	          errmsg = msprintf(gettext("%s: Value for relative tolerance between upper and lower bounds should be a Constant"), "ecos");
	          error(errmsg);
	        else
	          option(24) = param(2*i);
	        end

	    else
	      errmsg = msprintf(gettext("%s: Unrecognized parameter name %s."), "ecos", param(2*i-1));
	      error(errmsg);
	    end
  	end

    // number of second order cones is size of q
  	ncones = size(q,2);

    //Converting form 1-based(Scilab) indexing to 0-based indexing(C)
  	Gjc=Gjc-1;
    Gir=Gir-1;
  	if (p ~=0)
  		Ajc=Ajc-1;
  		Air=Air-1;
  	end

    // Calling ecos
  	[x,y,info1,s,z]=solveecos(c,Gpr,int32(Gjc),int32(Gir),h,Apr,int32(Ajc),int32(Air),b,l,int32(q),e,option,m,n,p,ncones)
  	
    // Assigning output parameters 
  	info=struct();
  	info.setup_time=info1(1);
  	info.solve_time = info1(2);
  	info.primal_objective_cost =info1(3);
  	info.dual_objective_cost =info1(4);
  	info.Normalized_primal_residual = info1(5);
  	info.Normalized_dual_residual = info1(6);
  	info.primal_infeasibile = info1(7);
  	info.dual_infeasibile = info1(8);
  	info.primal_infeasibility_measure = info1(9);
  	info.dual_infeasibility_measure = info1(10);
  	info.Complementarity_gap = info1(11);
  	info.Normalized_complementarity_gap = info1(12);
  	info.Iterations = info1(13);
  	info.exitflag = info1(14);

    // Status
  	select info.exitflag
  		case 0 then
  			status='Optimal solution found';
  		case 1 then
  			status='Certificate of primal infeasibility found';
  		case 2 then
  			status='Certificate of dual infeasibility found';
  		case 10 then
  			status='Optimal solution found subject to reduced tolerances';
  		case 11 then
  			status='Certificate of primal infeasibility found subject to reduced tolerances';
  		case 12 then
  			status='Certificate of dual infeasibility found subject to reduced tolerances';
  		case -1 then
  			status='Maximum number of iterations reached';
  		case -2 then
  			status='Numerical problems (unreliable search direction)';
  		case -3 then
  			status='Numerical problems (slacks or multipliers outside cone)';
  		case -4 then
  			status='Interrupted by signal or CTRL-C';
  		case -7 then
  			status='Unknown problem in solver';
      case -8 then
        status='ecos setup error';
  		else
  			status='Unknown problem in toolbox,contact toolbox authors';
      end
  	disp(status);
endfunction
