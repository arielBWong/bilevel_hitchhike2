%%

handle =  @(x)ackley(x);
num_vari = 2;
design_space = [ones(1,num_vari).*(-32); ones(1,num_vari)*32];
[bestmem1, bestval1, nfeval1] = DE(handle, num_vari, design_space(1,:), design_space(2,:),100, 100);

function y =  ackley(xx)
% assume x is two variables  

d = length(xx);

if (nargin < 4)
    c = 2*pi;
end
if (nargin < 3)
    b = 0.2;
end
if (nargin < 2)
    a = 20;
end

sum1 = 0;
sum2 = 0;
for ii = 1:d
	xi = xx(ii);
	sum1 = sum1 + xi^2;
	sum2 = sum2 + cos(c*xi);
end

term1 = -a * exp(-b*sqrt(sum1/d));
term2 = -exp(sum2/d);

y = term1 + term2 + a + exp(1);

end
