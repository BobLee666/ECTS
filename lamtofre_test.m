clear;
syms An Bn x;
An = 3.1432e9;
Bn = 1e18;
meiS = 1;
f0 = 2e10;
%a = 8.108*10^(-12);
a = 8.108*10^(-12);


strA = num2str(An);
strB = num2str(Bn);
stra = num2str(-a);
%x = solve(-a + An/(x^2) + 2*Bn/(x^3) == 0);
%x = solve(-8.108*10^(-12) + 3.1432e9/(x^2) + 2*1e18/(x^3) == 0);
equal = strcat([stra,'+',strA,'/(x^2)+2*',strB,'/(x^3)=0']);
x = solve(equal,'x');

y = x(3);
ans = double(y);
lam = callambda(meiS,An,Bn,f0);

