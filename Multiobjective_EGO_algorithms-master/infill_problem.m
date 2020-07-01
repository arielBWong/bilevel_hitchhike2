function [f,g] = infill_problem(x, fun_handle)


if nargin == 0
	prob.nx = 6;
	prob.nf = 1;
	prob.ng = 0;
	prob.range = cell(1, prob.nx);
	for i = 1:prob.nx
		prob.range{i} = Range('range', [0,1]);
	end
	f = prob;
else
	[f,g] = test_true(double(x), fun_handle);
end
return


function [f,g] = test_true(x, fun_handle)
f = feval(fun_handle, x);
g = [];
return
