function [f,g] = g8(x)
if nargin == 0
	prob.nx = 2;
	prob.nf = 1;
	prob.ng = 2;
	prob.range = {Range('range', [0,10]), Range('range', [0,10])};
	f = prob;
else
	[f,g] = g8_true(double(x));
end
return


function [f,g] = g8_true(x)
f = - sin(2*pi*x(1))^3 * sin(2*pi*x(2)) / (x(1)^3 * (x(1)+x(2)));
g(1) = -x(1)^2 + x(2) - 1;
g(2) = -1 + x(1) - (x(2)-4)^2;
return
