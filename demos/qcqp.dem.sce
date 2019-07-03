mode(1)
//
// Demo of qcqp.sci
//

//Reference : http://www.minlplib.org/nvs10.html
// var i1 := 1, >= 0, <= 200;
// var i2 := 1, >= 0, <= 200;
// minimize obj: 7*i1^2 + 6*i2^2 - 35*i1 - 80.4*i2;
// subject to
// e1: (-9*i1^2) - 10*i1*i2 - 8*i2^2 >= -583;
// e2: (-6*i1^2) - 8*i1*i2 - 6*i2^2 >= -441;
H = [7 0;0 6];
f = [-35 ; -80.4];
Q = [[9 5;5 8];[3 4;4 3]];
c = [0 0;0 0];
r = [583;441];
x = [1;1];
lb = [0;0];
ub = [200;200];
[xopt,fopt] = qcqp(x,H,f,Q,c,r,[],[],[],[],lb,ub);
// Press ENTER to continue
halt()   // Press return to continue
 
// var i1 integer := 1, >= 0, <= 200;
// var i2 integer := 1, >= 0, <= 200;
// var i3 integer := 1, >= 0, <= 200;
// minimize obj: 7*i1^2 + 6*i2^2 - 15.8*i1 - 93.2*i2 + 8*i3^2 - 6*i3*i1 + 4*i3*i2
//      - 63*i3;
// subject to
// e1: (-9*i1^2) - 10*i1*i2 - 8*i2^2 - 5*i3^2 - 6*i3*i1 - 10*i3*i2 >= -1000;
// e2: (-6*i1^2) - 8*i1*i2 - 6*i2^2 - 4*i3^2 - 2*i3*i1 - 2*i3*i2 >= -550;
// e3: (-9*i1^2) - 6*i2^2 - 8*i3^2 + 2*i2*i1 + 2*i3*i2 >= -340;
x = [1 ; 1; 1];
lb = [0 ; 0; 0];
ub = [200 ; 200 ; 200];
H = [7 0 -3;
0 6 2;
-3 2 8];
f = [-15.8 ; -93.2 ; -63];
Q = [
[9 5 3;
5 8 5;
3 5 8];
[3 4 1;
4 6 1;
1 1 4];
[9 -1 0;
-1 6 -1;
0 -1 8];
];
c = zeros(3,3);
r = [1000;550;340];
[xopt,fopt,lambda,exitflag] = qcqp(x,H,f,Q,c,r,[],[],[],[],lb,ub);
// Press Enter to continue
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
