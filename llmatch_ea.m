function[match_xl] = llmatch_ea(xu, prob, num_pop, num_gen)
% this function use gsolver to find a matching xl for given xu
%   usage
%       input
%           prob:  bilevel problem
%           num_pop: ea parameter population size
%           num_gen: ea parameter generation size
%           xu: upper x to be matched
%       ouput
%           match_xl: upper level x for xu
%           n_fev: number of function evaluation in local search
%----------
funh_obj = @(x)llobj(x, xu, prob);
funh_con = @(x)llcon(x, xu,  prob);

param = stuct;
param.gen = num_gen;
param.popsize = num_pop;

[bestx, bestf, bestc, archive] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param);
% check bestc is feasible or not
num_con = size(bestc, 2);
bestc(bestc <= 0) = 0;
s = sum(bestc, 2) == 0;

[match_xl, n_fev] = losolver(best_x, s, prob,  xu);

end

%objective wrapper for both local and global search
function f = llobj(xl, xu, prob)
[f, ~] = prob.evaluate_l(xu, xl);
end

%constraint wrapper for global search
function [c]  = llcon(xl, xu, prob)
[~, c] = prob.evaluate_l(xu, xl);
end



%constraint wrapper for local search
function [c, ceq]  = llconstraint(xl, xu, prob)
[~, c] = prob.evaluate_l(xu, xl);
ceq = [];
end


