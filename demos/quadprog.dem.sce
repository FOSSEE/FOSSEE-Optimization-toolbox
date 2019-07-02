mode(1)
//
// Demo of quadprog.sci
//

//Ref : example 14 :
//https://www.me.utexas.edu/~jensen/ORMM/supplements/methods/nlpmethod/S2_quadratic.pdf
// min. -8*x1*x1 -16*x2*x2 + x1 + 4*x2
// such that
//    x1 + x2 <= 5,
//    x1 <= 3,
//    x1 >= 0,
//    x2 >= 0
H = [2 0
0 8];
f = [-8; -16];
A = [1 1;1 0];
conUB = [5;3];
conLB = [-%inf; -%inf];
lb = [0; 0];
ub = [%inf; %inf];
nbVar = 2;
nbCon = 2;
[xopt,fopt,exitflag,output,lambda] = quadprog(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB)
//Press ENTER to continue
halt()   // Press return to continue
 
//Find x in R^6 such that:
A= [1,-1,1,0,3,1;
-1,0,-3,-4,5,6;
2,5,3,0,1,0
0,1,0,1,2,-1;
-1,0,2,1,1,0];
conLB=[1;2;3;-%inf;-%inf];
conUB = [1;2;3;-1;2.5];
lb=[-1000;-10000; 0; -1000; -1000; -1000];
ub=[10000; 100; 1.5; 100; 100; 1000];
//and minimize 0.5*x'⋅H⋅x + f'⋅x with
f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
nbVar = 6;
nbCon = 5;
x0 = repmat(0,nbVar,1);
param = list("MaxIter", 300, "CpuTime", 100);
[xopt,fopt,exitflag,output,lambda]=quadprog(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB,x0,param)
//========= E N D === O F === D E M O =========//
