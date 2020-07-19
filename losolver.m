function [match_xl, n_fev] = losolver(best_x, s, prob,  xu)
%
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
if s  % ego return feasible or unconstraint problem
    match_xl = newxl; 
    if best_f < newfl % if local search performance is not as good as ego
        match_xl = best_x;
    end
else % ego did not find feasible 
    match_xl = newxl;
    if output.constrviolation > 1e-6% local solver also fails
        flag = false;
        % neither ego or local search found feasible, return by smaller
        % constraint
        [~, newfc] = prob.evaluate_l(xu, newxl);
        if sum(newfc) > sum(best_c)
            match_xl = best_x;
        end
    end
end

n_fev = output.funcCount;

end

%objective wrapper
function f = llobjective(xl, xu, prob)
[f, ~] = prob.evaluate_l(xu, xl);
end

%constraint wrapper
function [c, ceq]  = llconstraint(xl, xu, prob)
[~, c] = prob.evaluate_l(xu, xl);
ceq = [];
end

