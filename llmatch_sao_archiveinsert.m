function[match_xl, n_fev, flag] = llmatch_sao_archiveinsert(xu, prob, num_pop, itersize, iter_freq, varargin)
% this lower level matching method uses population based sao to find a
% matching xl for upper level xu
% input:
%   xu: upper level variable to be matched
%   prob: real problem instance for true evaluation
%   num_pop: EA evolution parameter
%   num_gen: EA evolution parameter
%   iter_freq: Every how many generations when true evaluation and krg update  is applied
% output:
%   matching_xl : found xl for xu
%   n_fev : total number of function evaluation on lower level
%   flag : whether xl is a feasible solution(true/false)
%----
% steps:
%    (1) initialize xl population, evaluate and train kriging
%    (2) use kriging as objective function, evolve xl population
%    (3) when update frequence is met, evaluate current population,
%       expand krg training data with UNSEEN data, update gsolver objective function
%    (4) continue to evolve xl population, until num_gen is met
%-----------------------------------------------------------------

% (1) initialize xl population, evaluate and train kriging
l_nvar = prob.n_lvar;
upper_bound = prob.xl_bu;
lower_bound = prob.xl_bl;
% init_sizel = 2 * l_nvar + 1;
init_sizel = 11 * l_nvar - 1;
% init_sizel =7;
xu_init = repmat(xu, init_sizel, 1);
train_xl = lhsdesign(init_sizel,l_nvar,'criterion','maximin','iterations',1000);
train_xl = repmat(lower_bound, init_sizel, 1) ...
    + repmat((upper_bound - lower_bound), init_sizel, 1) .* train_xl;

% evaluate/get training fl from xu_init and train_xl
% compatible with non-constriant problem
[train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);

% lower level is considered as single objective
if size(train_fl, 2)>1
    error('lower level problem is multiple objective, algorithm is not compatible');
end

% (2) use kriging as objective function, evolve xl population
[krg_obj, krg_con, ~] = update_surrogate(train_xl, train_fl, train_fc);

% (3) when update frequence is met, evaluate current population,
%       expand krg training data, update gsolver objective function
param = struct();
initmatrix = train_xl;
initmatrix = [];
n = itersize;
for g = 1: n
    % fprintf('lower gen %d\n', g);
    funh_obj = @(x)llobj(x, krg_obj);
    funh_con = @(x)llcon(x, krg_con);

    param.gen=iter_freq;
    param.popsize = num_pop;
    % (4) continue to evolve xl population, until num_gen is met
    [~,~,~, archive] = gsolver(funh_obj, l_nvar,  prob.xl_bl, prob.xl_bu, initmatrix, funh_con, param);
       
    % last population is re-evaluated with real evaluation, only those
    % unseen in archive (training data)
    [train_xl, train_fl, train_fc, growflag] = ulego_sao_updateArchiveL(xu,archive.pop_last.X, prob, train_xl, train_fl, train_fc, true);    
    
    if growflag % there is unseen data in evolution
        [krg_obj, krg_con, ~] = update_surrogate(train_xl, train_fl, train_fc);
        initmatrix = [];
    else % there is no unseen data in evolution
            % re-introduce random individual 
        fprintf('no unseen data in last population, introduce randomness');
        initmatrix = [];
    end
end

% n = size(train_xl, 1);
% fprintf('ble lower evaluation size %s \n', num2str(n));

% local search for best xl
% connect a local search to sao
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


%----
% save lower level
llcmp = true;
%llcmp = false;
if llcmp
    method = 'llmatchble';
    seed = varargin{1};
    % add local search result
    train_xl = [train_xl; match_xl];
    [local_fl, local_fc]  = prob.evaluate_l(xu, match_xl);
    train_fl = [train_fl; local_fl];
    train_fc = [train_fc; local_fc];
    
    perfrecord_umoc(xu, train_xl, train_fl, train_fc, prob, seed, method, 0, 0);
end 

end


function [sf, sx, sc] = initmatrix_pick(x, f, c)
    [sf, sx, sc] = pop_sort(f, x, c);
end

function  f = llobj(x, kriging_obj)
num_obj = length(kriging_obj);   % krg cell array?
num_x = size(x, 1);
f = zeros(num_x, num_obj);
for ii =1:num_obj
    [f(:, ii), ~] = dace_predict(x, kriging_obj{ii});
end
end

function c = llcon(x, krging_con)
num_con = length(krging_con);
num_x = size(x, 1);
c = zeros(num_x, num_con);
for ii =1:num_con 
    [c(:, ii), ~] = dace_predict(x, krging_con{ii});
end
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



%-----auxiliary function ---
function [krg_obj, krg_con, info] = update_surrogate(trainx, trainy, trainc)
% this function updates train kriging model from train x and train y
%
%
if size(trainy,2)>1
    error(' following zscore norm only applies to single objective problems');
end
[train_y_norm, y_mean, y_std] = zscore(trainy, 1, 1);
num_obj = size(trainy, 2);
krg_obj = cell(1, num_obj);
num_vari = size(trainx, 2);
for ii = 1:num_obj
    % kriging_obj{ii} = dace_train(train_x_norm,train_y_norm(:,ii));
    krg_obj{ii} = dacefit(trainx,train_y_norm(:,ii),...
        'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  % for test
end

info = struct();
info.ymean = y_mean;
info.ystd = y_std;

% deal with constraints
if ~isempty(trainc)
    num_con = size(trainc, 2);
    krg_con = cell(1, num_con);
    [train_c_norm, c_mean, c_std] = zscore(trainc, 1, 1);
    
    for ii = 1:num_con
        krg_con{ii} = dacefit(trainx, train_c_norm(:,ii),...
            'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  % for test
    end
    info.cmean = c_mean;
    info.cstd = c_std;
else
    krg_con = [];
    info.cmean = [];
    info.cstd = [];
end

end