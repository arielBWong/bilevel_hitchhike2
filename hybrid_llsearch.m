function [newxl, n_feval, flag] = hybrid_llsearch(xu, xl_start, prob, hy_pop, hy_gen)
% This function uses a starting point of xl
% to further improve xu's match on lower level
% usage:    
% input     xu             : the final xu for re-evaluation
%           xl_start       : xl matched from upper ego
%           prob           : problem instance
%           hy_pop         : lower global search parameter
%           hy_gen         : lower global search parameter
% output    newxl          : xl returned from re-evaluation
%           n_feval        : number of lower level function evaluation
%           flag           : indicating whether newxl is feasible (true/false)
%--------------------------------------------------------------------------
% use xl_start as initial point

%-global search
obj = @(x)hyobj(x, xu, prob);
con = @(x)hycon(x, xu, prob);

opts = optimoptions('ga');
opts.MaxGenerations = hy_gen;
opts.InitialPopulationMatrix = xl_start;
opts.PopulationSize = hy_pop;

[newxl_g, newfl_g] = ga(obj, prob.n_lvar, [],[],[],[], prob.xl_bl, prob.xl_bu,con, opts);

%-local search
opts = optimset('fmincon');
opts.Algorithm = 'sqp';
opts.MaxFunctionEvaluations = 100;
[newxl_l, newfl_l, ~, output] = fmincon(obj, newxl_g, [], [],[], [],  ...
    prob.xl_bl, prob.xl_bu, con, opts);

%-check and return results
newxl = newxl_l;
n_feval = hy_pop * hy_gen + output.funcCount;
flag = true;
if output.constrviolation > 0
    flag = False;
    warning('hybrid search did not find any feasible lower match');
end
end

function [f] = hyobj(x, xu, prob)
    [f, ~] = prob.evaluate_l(xu, x);
end

function [c, ceq] = hycon(x, xu, prob)
    [~, c] = prob.evaluate_l(xu, x);
    ceq = [];
end