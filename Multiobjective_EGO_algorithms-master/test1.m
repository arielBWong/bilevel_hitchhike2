function [f,g] = test1(x)
if nargin == 0
	prob.nx = 4;
	prob.nf = 1;
	prob.ng = 1;
	prob.range = cell(1, prob.nx);
	for i = 1:prob.nx
		prob.range{i} = Range('range', [-10,10]);
	end
	f = prob;
else
	[f,g] = test_true(double(x));
end
return


function [f,g] = test_true(x)
f(1)=x(1)^3+x(2)^2-(x(3)-3)^2+100;
g(1) = 5-(x(1)+x(2));
return

