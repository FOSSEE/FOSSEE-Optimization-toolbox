# FOSSEE Optimization Toolbox for Scilab 6.0.x

A toolbox that provides mixed integer linear programming, quadratic programming and non linear programming tools in Scilab through the Symphony and Ipopt libraries.

NOTE: On linux systems with gfortran8 as the default version, the user will need to install libgfortran4 for the toolbox to load. This can be done by sudo apt-get install libgfortran4 (for systems with aptitude as the package manager). 


## To Download
1. [Visit the link
   `http://atoms.scilab.org/toolboxes/FOT/`]
2. Select the linux, macOS, or windows version as per your platform.
3. Extract the files.

## To use
1. In Scilab, change the working directory to the root directory of the repository
2. Run `exec loader.sce` in the scilab console.
3. The Toolbox is now ready, to see help type `help FOSSEE Optimization Toolbox` in console.
4. The demos are in `Demos folder`.
5. To run a demo type `exec <name of function>.dem.sce`
6. Test cases are in `tests folder`.

## To build
1. If you have updated the source code then you need to build it again to see the changes.
2. To build it first unlink the toolbox by executing the following command `ulink`.
3. Then type `exec builder.sce` to run the builder. {Prerequisites: In windows you need Visual Studio}
4. Now just run `exec loader.sce` in the scilab console. It will be ready to use.

   This toolbox consists of open-source solvers for a variety of optimization
problems: CLP for linear optimization, CBC and Symphony for integer linear
optimization, IPOPT (with MUMPS) for nonlinear optimization, and BONMIN for
integer nonlinear optimization.

Features
---------
* linprog: Solves a linear optimization problem.
 	
* intlinprog: Solves a mixed-integer linear optimization problem in intlinprog
format with CBC.
  
* symphony: Solves a mixed-integer linear optimization problem.
  
* symphonymat: Solves a mixed-integer linear optimization problem (with input
  in Matlab format).
  
* quadprog: Solves a quadratic optimization problem.
  
* quadprogmat: Solves a quadratic optimization problem (with input in Matlab
  format).
  
* lsqnonneg: Solves a nonnegative linear least squares optimization problem.
  
* lsqlin: Solves a linear least squares optimization problem.
  
* lsqnonlin: Solves a nonlinear least squares optimization problem.
  
* fminunc: Solves an unconstrained optimization problem.
  
* fminbnd: Solves a nonlinear optimization problem on bounded variables.
 
* fmincon: Solves a general nonlinear optimization problem.
  
* fgoalattain: Solves a multiobjective goal attainment problem.
  
* fminimax: Solves a minimax optimization problem.
  
* intfminunc: Solves an unconstrained mixed-integer nonlinear optimization
  problem.
  
* intfminbnd: Solves a mixed-integer nonlinear optimization
  problem on bounded variables.
  
* intfmincon: Solves a constrained mixed-integer nonlinear optimization
problem.
  
* intfminimax: Solves a mixed-integer minimax optimization problem.
  
* intquadprog: Solves an integer quadratic optimization problem.

* quadprogCLP: Solves a quadratic optimization problem.

* qcqp: Solves a quadratic constrained quadratic optimization problem.
