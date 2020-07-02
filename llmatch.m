function[match_xl, n_fev, flag] = llmatch(xu, prob, num_pop, num_gen, init_size, iter_size)
% method of searching for a match xl for xu
% usage:
%
% input: xu             - upper level variable, to be matched 
%        prob           - problem instance
%        num_pop        - DE parameter
%        num_gen        - DE parameter
%        init_size      -surrogate parameter: number of initiliazation samples
%        iter_size      -surrogate parameter: number of iterations
%
% output: matching_xl   - found xl for xu
%         n_fev         - total number of function evaluation on lower
%                           level
%         flag          - whether xl is a feasible solution
%--------------------------------------------------------------------------

l_nvar = prob.n_lvar;
upper_bound = prob.xl_bu;
lower_bound = prob.xl_bl;
xu_init = repmat(xu, init_size, 1);
train_xl = lhsdesign(init_size,l_nvar,'criterion','maximin','iterations',1000);
train_xl = repmat(lower_bound, init_size, 1) + repmat((upper_bound - lower_bound), init_size, 1) .* train_xl;

% test only---
testname = 'trainx.csv';
csvwrite(testname,train_xl);
% test only-----

% evaluate training fl from xu_init and train_xl
% compatible with non-constriant problem
[train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);
l_ncons= size(train_fc, 2);

% call EIM next
for iter = 1:iter_size
    % eim propose next
    [new_xl, ~] = EIMnext_znorm(train_xl, train_fl, upper_bound, lower_bound, ...
        num_pop, num_gen, train_fc);
    % evaluate next
    [new_fl, new_fc] = prob.evaluate_l(xu, new_xl);
    
    % add to training
    train_xl = [train_xl; new_xl];
    train_fl = [train_fl; new_fl];
    train_fc = [train_fc; new_fc];  %compatible with nonconstraint
end

% provide local search with best solution
if ~isempty(train_fc)
    ind_feas = sum(train_fc<=0, 2)== l_ncons;
    if sum(ind_feas) ~= 0  % feasible solution exists
        train_xl = train_xl(ind_feas, :);
        train_fl = train_fl(ind_feas, :);
        train_fc = train_fc(ind_feas, :);
    end
    % if there is no feasible select with smallest fl
    % same as non constraint problems, compatibility checked
end

% local search starting point selelction
% lower level is considered as single objective
[~, best_ind] = min(train_fl);
if length(best_ind) > 1
    best_ind = best_ind(0); end
best_x = train_xl(best_ind, :);
% best_f = train_fl(best_ind, :);
% compatible with constraint problem

% give starting point to local search
fmin_obj = @(x)llobjective(x, xu, prob);
fmin_con = @(x)llconstraint(x, xu, prob);
opts = optimset('fmincon');
opts.Algorithm = 'sqp';
opts.MaxFunctionEvaluations = 100;
[newxl, newfl, flag, output] = fmincon(fmin_obj, best_x, [], [],[], [],  ...
    prob.xl_bl, prob.xl_bu, fmin_con,opts);
newxl
newfl
% decide which x to return
% compatible with constraint problem
flag = true;
if ~isempty(train_fc)&& sum(ind_feas) == 0 && output.constrviolation > 0
% for constraint problems
% surrogate has no feasible andlocal search also return no feasible
flag = false; 
end         
    
match_xl = newxl; 
% count local search number
n_fev = iter_size + output.funcCount;
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




