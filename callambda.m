function [lam]=callambda(S,An,Bn,f0)
a = 2*Bn/(f0^3)+An/(f0^2);
zeng = 1e-12;
delta = 1e-13;
%a = a + zeng;
while(true)
	fsum = fresum(S,An,Bn,a);
	if(fsum > f0)
		a = a+zeng;
		zeng = zeng*2;
		continue;
    elseif(fsum == f0)
        lam = a;
        break;    
    else
		astart = a-zeng;
		aend = a;
	end
	while(aend-astart > delta)
		amid = (astart+aend)/2;
		fsum = fresum(S,An,Bn,amid);
		if(fsum > f0)
			astart = amid;
		else
			aend = amid;
		end
    end
	lam = aend;
    break;
end
end




