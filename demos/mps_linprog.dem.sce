mode(1)
//
// Demo of mps_linprog.sci
//

file1_link="https://github.com/bhanusayal/symphony/blob/master/demos/p0033.mps";
[xopt,fopt,exitflag,output,lambda] =lp_linprog(file1);
file2_link="https://github.com/bhanusayal/symphony/blob/master/demos/exmip1.mps";
[xopt,fopt,exitflag,output,lambda] =lp_linprog(file2);
//========= E N D === O F === D E M O =========//
