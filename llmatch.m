function[match_xl, n_fev, flag] = llmatch(xu, prob, num_pop, num_gen, propose_nextx, iter_size, llfit_hn,  varargin)
% method of searching for a match xl for xu. 
% Problem(Prob) definition require certain formation for bilevel problems
% evaluation method should be of  form 'evaluation_l(xu, xl)'
% usage
% input: 
%        xu:  upper level variable, to be matched
%        prob: problem instance, require certin evaluation
%                      method name defintion-- check problems
%                      folder
%        num_pop: DE parameter
%        num_gen : DE parameter
%        propose_nextx  : str, function handle to generate next x
%        iter_size : surrogate parameter: number of iterations
%        llfit_hn :  str, lower level seach fitness for next x method, used
%                          by proposed_nextx method 
%
% output: 
%        matching_xl : found xl for xu 
%         n_fev : total number of function evaluation on lower level
%         flag : whether xl is a feasible solution(true/false)
%--------------------------------------------------------------------------

l_nvar = prob.n_lvar;
% init_size = 11 * l_nvar -1;
init_size = 20;
upper_bound = prob.xl_bu;
lower_bound = prob.xl_bl;
xu_init = repmat(xu, init_size, 1);
train_xl = lhsdesign(init_size,l_nvar,'criterion','maximin','iterations',1000);
train_xl = repmat(lower_bound, init_size, 1) ...
                + repmat((upper_bound - lower_bound), init_size, 1) .* train_xl;

% evaluate/get training fl from xu_init and train_xl
% compatible with non-constriant problem
[train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);

% lower level is considered as single objective
if size(train_fl, 2)>1
    error('lower level problem is multiple objective, algorithm is not compatible');
end

% call EIM/Ehv to expand train xl one by one
fithn = str2func(llfit_hn);
nextx_hn = str2func(propose_nextx);
for iter = 1:iter_size
    % eim propose next xl
    % lower level is single objective so no normalization method is needed
    % tic;
    [new_xl, ~] = nextx_hn(train_xl, train_fl, upper_bound, lower_bound, ...
        num_pop, num_gen, train_fc, fithn);
   %  toc;
    
    % evaluate next xl with xu
    [new_fl, new_fc] = prob.evaluate_l(xu, new_xl);
    
    % add to training
    train_xl = [train_xl; new_xl];
    train_fl = [train_fl; new_fl];
    train_fc = [train_fc; new_fc];  %compatible with nonconstraint
end

% connect a local search to ego
% local search starting point selection
[best_x, best_f, best_c, s] =  localsolver_startselection(train_xl, train_fl, train_fc);

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

% count local search number
n_fev = init_size + iter_size + output.funcCount;

% save lower level
llcmp = true;
if llcmp
    method = 'llmatcheim';
    seed = varargin{1};
    % add local search result
    train_xl = [train_xl; match_xl];
    [local_fl, local_fc]  = prob.evaluate_l(xu, match_xl);
    train_fl = [train_fl; local_fl];
    train_fc = [train_fc; local_fc];
    
    perfrecord_umoc(train_xl, train_fl, train_fc, prob, seed, method, 0, 0)
end

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




