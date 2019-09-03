# FOSSEE Optimization Toolbox for Scilab 6.0.x

A toolbox that provides mixed-integer programming, quadratic programming and nonlinear programming tools in Scilab through various open-source libraries available from Coin-OR.

NOTE: On linux systems with gfortran8 as the default version, the user will need to install libgfortran4 for the toolbox to load. This can be done, for example in Ubuntu, by executing: sudo apt-get install libgfortran4


## To Download
1. [Visit the link
   `http://atoms.scilab.org/toolboxes/FOT/`]
2. Select the linux or windows version as per your platform.
3. Extract the files.

## To use
1. In Scilab, change the working directory to the root directory of the repository
2. Run `exec loader.sce` in the scilab console.
3. The Toolbox is now ready, to see help type `help` in console.
4. The demos are available in `Demos folder`.
5. To run a demo type `exec <name of function>.dem.sce`
6. Test cases are available in `tests folder`.

## To build
1. If you have updated the source code, you need to build it again to see the changes.
2. To build it first unlink the toolbox by executing the command `ulink`.
3. Then type `exec builder.sce` to run the builder. {Prerequisites: In windows you need Visual Studio}
4. Now run `exec loader.sce` in the scilab console. The toolbox will be ready
   to use.

   This toolbox consists of open-source solvers for a variety of optimization
problems: CLP for linear and quadratic optimization, CBC for integer linear
optimization, IPOPT (with MUMPS) for nonlinear optimization, and BONMIN for
integer nonlinear optimization.

Features
---------
* linprog: Solves a linear optimization problem.
 	
* intlinprog: Solves a mixed-integer linear optimization problem in intlinprog
format with CBC.
  
* quadprog: Solves a quadratic optimization problem.
  
* quadprogmat: Solves a quadratic optimization problem (with input in Matlab
  format).
  
* quadprogCLP: Solves a quadratic optimization problem.

* intquadprog: Solves an integer quadratic optimization problem.

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
  
