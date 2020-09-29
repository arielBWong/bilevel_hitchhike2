function[match_xl, n_fev, flag] = llmatch_eimfix(xu, prob, num_pop, num_gen, propose_nextx, iter_size, llfit_hn,  varargin)
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


% use believer to deal with the last one
norm = str2func('normalization_z');
[krg_obj, krg_con, ~] = update_surrogate(train_xl, train_fl, train_fc, norm);
funh_obj = @(x)llobj(x, krg_obj);
funh_con = @(x)llcon(x, krg_con);

param.gen       = num_gen;
param.popsize   = num_pop;
[~,~,~, archive] = gsolver(funh_obj, l_nvar,  prob.xl_bl, prob.xl_bu, [], funh_con, param);
[train_xl, train_fl, train_fc, growflag] = ulego_sao_updateArchiveL(xu,archive.pop_last.X, ...
     prob, train_xl, train_fl, train_fc, true);    

 
 
% connect a local search to ego
% local search starting point selection
[best_x, best_f, best_c, s] =  localsolver_startselection(train_xl, train_fl, train_fc);
nolocalsearch = true;
if nolocalsearch
    match_xl = best_x;
    n_fev = size(train_xl, 1);
    flag = s;
else
    [match_xl, flag, num_eval] = ll_localsearch(best_x, best_f, best_c, s, xu, prob);
    n_global = size(train_xl, 1);
    n_fev = n_global +num_eval;       % one in a population is evaluated
end




% save lower level
llcmp = true;
%llcmp = false;
if llcmp
    method = 'llmatcheimfix';
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



% believer objectives 
function  f = llobj(x, kriging_obj)
num_obj = length(kriging_obj);   % krg cell array?
num_x = size(x, 1);
f = zeros(num_x, num_obj);
for ii =1:num_obj
    [f(:, ii), ~] = dace_predict(x, kriging_obj{ii});
end
end


% believer constraints
function c = llcon(x, krging_con)
num_con = length(krging_con);
num_x = size(x, 1);
c = zeros(num_x, num_con);
for ii =1:num_con 
    [c(:, ii), ~] = dace_predict(x, krging_con{ii});
end
end




