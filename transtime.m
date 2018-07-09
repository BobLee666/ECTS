function [time]=transtime(L0,v,theta,tpro,D,W,p)
sum=0;
t=0;
L0 = sqrt(L0^2+(v*tpro)^2-2*L0*v*tpro*cosd(theta));
while(sum<D)
	L = L0^2+(v*t)^2-2*L0*v*t*cosd(theta);
	delta = 0.001*W*log2(1+p*L^(-2)/1e-13);
	t = t+0.001;
	sum = sum+delta;
end
time = t;
end
