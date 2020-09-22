function [match_xl, flag, num_eval] = ll_localsearch(best_x, best_f, best_c, s, xu, prob)
% This method is part of local search algorithm
% receives starting point
%
% input:
%       best_x: starting point
%       best_f: starting point's objective value
%       best_c: starting point's constraint value
%       s: flag, whether this solution is feasible
% output:
%       match_xl: the final matching xl for a given xu
%       flag: whether match_xl is feasible
%       num_eval: number of function evaluations consumed in local search
%-----------------------------------------------------------


% give starting point to local search
fmin_obj = @(x)llobjective(x, xu, prob);
fmin_con = @(x)llconstraint(x, xu, prob);
opts = optimset('fmincon');
opts.Algorithm = 'sqp';
opts.Display = 'off';
opts.MaxFunctionEvaluations = 100;
[newxl, newfl, ~, output] = fmincon(fmin_obj, best_x, [], [],[], [],  ...
    prob.xl_bl, prob.xl_bu, fmin_con,opts);

% decide which xl to return back to upper level
% compatible with unconstraint problem
flag = true;
if s  % sao return feasible or unconstraint problem
    match_xl = newxl;
    if best_f < newfl % if local search performance is not as good as ego
        match_xl = best_x;
    end
else % sao did not find feasible
    match_xl = newxl;
    if output.constrviolation > 1e-6% local solver also fails
        flag = false;
        % neither global search or local search found feasible, return by smaller
        % constraint
        [~, newfc] = prob.evaluate_l(xu, newxl);
        if sum(newfc) > sum(best_c)
            match_xl = best_x;
        end
    end
end
num_eval = output.funcCount; 
end


%objective wrapper for true evaluation
function f = llobjective(xl, xu, prob)
[f, ~] = prob.evaluate_l(xu, xl);
end

%constraint wrapper
function [c, ceq]  = llconstraint(xl, xu, prob)
[~, c] = prob.evaluate_l(xu, xl);
ceq = [];
end
