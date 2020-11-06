function[match_xl, n_fev, flag] = llmatch_sao_archiveinsert(xu, prob, num_pop, init_size, itersize, num_gen, varargin)
% this lower level matching method uses population based sao to find a
% matching xl for upper level xu
% input:
%   xu: upper level variable to be matched
%   prob: real problem instance for true evaluation
%   num_pop: EA evolution parameter
%   num_gen: EA evolution parameter
%   num_gen: Every how many generations when true evaluation and krg update  is applied
% output:
%   matching_xl : found xl for xu
%   n_fev : total number of function evaluation on lower level
%   flag : whether xl is a feasible solution(true/false)
%-----------------------------------------------------------------

% (1) initialize xl population, evaluate and train kriging
norm_str = varargin{2};
norm_hn = str2func(norm_str);
l_nvar = prob.n_lvar;
upper_bound = prob.xl_bu;
lower_bound = prob.xl_bl;

% init_size = 11 * l_nvar - 1;
xu_init = repmat(xu, init_size, 1);
train_xl = lhsdesign(init_size,l_nvar,'criterion','maximin','iterations',1000);
train_xl = repmat(lower_bound, init_size, 1) ...
    + repmat((upper_bound - lower_bound), init_size, 1) .* train_xl;

initx = train_xl;

% evaluate/get training fl from xu_init and train_xl
% compatible with non-constriant problem
[train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);

% (2) use kriging as objective function, evolve xl population
[krg_obj, krg_con, ~] = update_surrogatedace(train_xl, train_fl, train_fc, norm_hn);

% (3) when update frequence is met, evaluate current population,
%       expand krg training data, update gsolver objective function
param = struct();
initmatrix = train_xl;
initmatrix = [];
n = itersize;
process_upper = false;
process_believer = true;

fighn = figure(1);


for g = 1: n    
    dace_stable  =  stability_check(train_xl, train_fl, train_fc, krg_obj, krg_con);
    if ~dace_stable
        fprintf('not stable at iteration %d\n', iter);
        break
    end
    %--test demo
    processplot(fighn, train_xl, train_fl, krg_obj, prob, initx);
    %--test demo
    
    
    % fprintf('lower gen %d\n', g);
    funh_obj = @(x)llobj(x, krg_obj);
    funh_con = @(x)llcon(x, krg_con);
    
    param.gen=num_gen;
    param.popsize = num_pop;
    % (4) continue to evolve xl population, until num_gen is met
    [~,~,~, archive] = gsolver(funh_obj, l_nvar,  prob.xl_bl, prob.xl_bu, initmatrix, funh_con, param);    
    [newx, growflag] = believer_select(archive.pop_last.X, train_xl, prob, process_upper, process_believer);
    
    if growflag % there is unseen data in evolution
        [new_xl] = surrogate_localsearch(xu, newx, prob, train_xl, train_fl, train_fc, norm_str);
        % new_xl = newx;
        [new_fl, new_fc] = prob.evaluate_l(xu, new_xl);
        
        % inprocess_plotsearch(fighn, prob, cons_hn, new_xl, train_xl);
        
        train_xl = [train_xl; new_xl];
        train_fl = [train_fl; new_fl];
        train_fc = [train_fc; new_fc];
        [krg_obj, krg_con, ~] = update_surrogatedace(train_xl, train_fl, train_fc, norm_hn);
        initmatrix = [];
    else % there is no unseen data in evolution
        % re-introduce random individual
        fprintf('no unseen data in last population, introduce randomness \n');
        initmatrix = [];
    end
    
    
end

% local search for best xl
% connect a local search to sao
% local search starting point selection
% MO problem no process
[best_x, best_f, best_c, s] =  localsolver_startselection(train_xl, train_fl, train_fc);
nolocalsearch = true;
if nolocalsearch
    match_xl = best_x;
    n_fev = size(train_xl, 1);
    flag = s;
else
    if size(train_fl, 2)> 1
        error('local search does not apply to MO');
    end
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
    % only for SO
    if size(train_fl, 2) ==  1
        train_xl = [train_xl; match_xl];
        [local_fl, local_fc]  = prob.evaluate_l(xu, match_xl);
        train_fl = [train_fl; local_fl];
        train_fc = [train_fc; local_fc];
    end    
    perfrecord_umoc(xu, train_xl, train_fl, train_fc, prob, seed, method, 0, 0, init_size);
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

%-----------paper demo-------------
function  [f, sig] = llobjdemo(x, kriging_obj)
num_obj = length(kriging_obj);   % krg cell array?
num_x = size(x, 1);
f = zeros(num_x, num_obj);
for ii =1:num_obj
    [f(:, ii), sig] = dace_predict(x, kriging_obj{ii});
end
end

% ---paper demo-
function[] = processplot(fighn, trainx, trainy, krg, prob, initx)
clf(fighn);

% (1) create test
testdata = linspace(prob.xl_bl, prob.xl_bu, 2000);
testdata = testdata';

% (2) predict
[fpred, sig] = llobjdemo(testdata, krg);
fpred = denormzscore(trainy, fpred);
plot(testdata, fpred, 'r--', 'LineWidth', 2); hold on;

% y1 = fpred + sig * 1.5;
% y2 = fpred - sig * 1.5;
% y = [y1', fliplr(y2')];
% x = [testdata', fliplr(testdata')];
% fill(x, y, 'r', 'FaceAlpha', 0.1, 'EdgeColor','none'); hold on;

% (3) real
[freal, sig]= prob.evaluate_l([], testdata);
plot(testdata, freal, 'b', 'LineWidth', 2);hold on;


% (4) scatter train
scatter(trainx, trainy, 80, 'ko', 'LineWidth', 2);

inity = prob.evaluate_l([], initx);
scatter(initx, inity, 40, 'ro', 'filled');


pause(1);

end

 function f = denormzscore(trainy, fnorm)
[train_y_norm, y_mean, y_std] = zscore(trainy);
f = fnorm * y_std + y_mean;
 end