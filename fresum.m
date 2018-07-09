function [fsum]=fresum(S,A,B,a)
	fsum = 0;
	for j=1:length(S)
		t = S(j);
        strA = num2str(A(t));
        strB = num2str(B(t));
        stra = num2str(-a(t));
        equal = strcat([stra,'+',strA,'/(x^2)+2*',strB,'/(x^3)=0']);
        x = solve(equal,'x');

        y = x(3);
        ftem = double(y);
		fsum = fsum + ftem;
    end
end
