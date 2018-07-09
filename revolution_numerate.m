clear;
clc;
N = 1;
R = 100;


betaT = 0.25+(0.75-0.25)*rand(N,1);  %user preference for T
betaE = 0.25+(0.75-0.25)*rand(N,1);	 %user preference for E
Fl = (0.5+(1.5-0.5)*rand(N,1))*1e9;
%L0 = R*rand(N,1);          %distance to base station
L0 = 0.8*R*ones(N,1);
v = 20*rand(N,1);
%v = 20*ones(N,1);
%theta = 180*rand(N,1);
theta = 180*ones(N,1);

tpro = 1;
C = 1000e6;
f0 = 20e9;
D = 420*1024*8; %bit
B = 20e6;
W = 1e6;
rou = 1;        %user payment
p = 10e-3*ones(N,1);
Tl = C./Fl;
El = 1e-20*Fl*C;
tao = rou*betaT;



%caculate the transtime 
for i=1:N
	Ttrans(i) = transtime(L0(i),v(i),theta(i),tpro,D,W,p(i));
end
% enumerate
for i = 2:2^N
    location = zeros(N,1);
    bi = dec2bin(i-1);
    st = num2str(bi);
    for j = 1:length(st)
        location(j) = (st(length(st)-j+1)=='1');
    end
    
    if sum(location==1)>B/W
        continue;
    end
    %translate location to schedule
    meiS = find(location==1);
    % calculate fi according to schedule
    su = 0; 
    for j=1:length(meiS)
        su = su + sqrt(tao(meiS(j))*Fl(meiS(j)));
    end
    for j=1:length(meiS)
        k = meiS(j);
        fi(k) = sqrt(tao(k)*Fl(k))*f0/su;
    end
	% calculate myf according to schedule
	An = zeros(N,1);
	Bn = zeros(N,1);
	
	for j=1:length(meiS)
        t = meiS(j);
		An(t) = (rou*betaT(t)*C/Tl(t)) + 2*(tpro+Ttrans(t))*C;
		Bn(t) = C^2;
    end
	lam = callambda(meiS,An,Bn,f0);
	for j=1:length(meiS)
		t = meiS(j);
        strB = num2str(B(t));
        stra = num2str(-a(t));
        equal = strcat([stra,'+',strA,'/(x^2)+2*',strB,'/(x^3)=0']);
        x = solve(equal,'x');

        y = x(3);
        myf(t) = double(y);
    end
end



