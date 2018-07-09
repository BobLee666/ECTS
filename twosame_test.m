clear;
N = 3;
An = 3.1432e9*ones(N,1);
Bn = 1e18*ones(N,1);
meiS = (1:N);
f0 = 2e10;
a = 8.108*10^(-12)*ones(N,1);
lam = callambda(meiS,An,Bn,f0);

