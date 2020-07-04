function[best_x, info] = EIMnext_cv(train_x, train_y, xu_bound, xl_bound, ...
    num_pop, num_gen, train_c)
% method of using EIM to generate next point with crossvalidation on krg
% normalization is zscore on all variables
% usage:
%
% input: train_x    - design variables
%                       1/2d array: (num_samples, num_varibles)
%        train_y    - objective values
%                       1/2d array: (num_samples, num_objectives)
%        xu_bound   - upper bound of train_x
%                       1d array
%        xl_bound   - lower bound of train_x
%                       1d array
%        num_pop    - EIM optimization parameter
%        num_gen    - EIM optimization parameter
%        train_c    - constraints values
%                       1/2d array: (num_samples, num_constraints)
% output: best_x    - proposed next x to be evaluated by EIM
%         info      - returned information for functor caller to recreate
%                   - or check information
%                   - info.krg
%                   - info.krgc
%                   - info.train_xmean
%                   - info.train_ymean
%                   - info.train_xstd
%                   - info.train_ystd
%                   - info.info.eim_normf
%--------------------------------------------------------------------------

% number of objective
num_obj = size(train_y, 2);
kriging_obj = cell(1,num_obj);

% number of input designs
num_x = size(train_x,1);
num_vari = size(train_x, 2);

% normalise train data
[train_x_norm, x_mean, x_std] = zscore(train_x, 1, 1);  % flag 0 sample zscore

% adjust x bounds according to zscore normalization
ub = (xu_bound - x_mean)./ x_std;
lb = (xl_bound - x_mean)./ x_std;


% train objective kriging model
if num_obj > 1
    [train_y_norm, y_mean, y_std] = zscore(train_y, 1, 1);
    % further normalization ??
    train_y_norm = (train_y_norm -repmat(min(train_y_norm),num_x,1))./...
        repmat(max(train_y_norm)-min(train_y_norm),num_x,1);
else
    [train_y_norm, y_mean, y_std] = zscore(train_y, 1, 1);
end

% prepare f_best for EIM, first only consider non-constraint situation
if num_obj > 1
    index_p = Paretoset(train_y_norm);
    f_best = train_y_norm(index_p, :); % constraint problem has further process
else
    f_best = min(train_y_norm, [], 1);
end

% compatibility with constraint problems
% if nargin>6
if ~isempty(train_c)
    % constraints should not be normalised
    % version did not scale train_c
    train_c_norm = train_c;
    num_con = size(train_c, 2);
    % adjust f_best according to feasibility
    % feasibility needs to be valued in original range
    index_c = sum(train_c <= 0, 2) == num_con;
    if sum(index_c) == 0 % no feasible
        f_best = [];
    else
        feasible_trainy_norm = train_y_norm(index_c, :);
        % still needs nd front
        % also works for single objective constraint problem
        index_p = Paretoset(feasible_trainy_norm);
        f_best = feasible_trainy_norm(index_p, :);
    end
    % there needs a model selection before dive into EIM
    [kriging_obj, kriging_con] = surrogate_build(train_x_norm, train_y_norm, train_c_norm);
    fitness_val = @(x)EIM_eval(x,f_best, kriging_obj, kriging_con);
else
    [kriging_obj, ~] = surrogate_build(train_x_norm, train_y_norm, []);
    fitness_val = @(x)EIM_eval(x, f_best, kriging_obj);
end

% call DE evolution
[best_x, eim_f] = DE(fitness_val, num_vari, lb, ub, num_pop, num_gen);

% convert best_x to denormalized value
best_x = best_x .* x_std + x_mean;  % commit for test

% fix bound violation
best_x = fixbound_violation(best_x, xu_bound, xl_bound);

%--
info = struct();
info.eim_normf = eim_f;
info.krg = kriging_obj;
info.train_xmean = x_mean;
info.train_xstd = x_std;
info.train_ymean = y_mean;
info.train_ystd = y_std;
end

function x = fixbound_violation(x, xu, xl)
% local function for znorm correctness, if x fall out of boundary
% move them on boudary
l_vio = x < xl;
x(l_vio) = xl(l_vio);

u_vio = x > xu;
x(u_vio) = xu(u_vio);
end

function [fit] = EIM_eval(x, f, kriging_obj, kriging_con)
% function of using EIM as fitness evaluation
% usage:
%
% input: x            - pop to evaluate
%        f            - best f so far/feasible pareto front
%                           in multi-objective probs
%        kriging_obj  - kriging model for objectives
%        kriging_con  - kriging model for constraints
% output: fit         - pop fitness
%--------------------------------------------------------------------------


% number of input designs
num_x = size(x,1);
num_obj = size(f, 2);

% the kriging prediction and varince
u = zeros(num_x,num_obj);
mse = zeros(num_x,num_obj);

% pof init
pof = 1;

if length(f) == 0 && nargin > 3   % refer no feasible solution
    % the number of constraints
    num_con = length(kriging_con);
    % the kriging prediction and varince
    mu_g = zeros(size(x,1), num_con);
    mse_g = zeros(size(x,1), num_con);
    for ii = 1: num_con
        [mu_g(:, ii), mse_g(:, ii)] = dace_predict(x, kriging_con{ii});
        % [mu_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_con{ii}); %test
    end
    pof  = prod(Gaussian_CDF((0-mu_g)./mse_g), 2);
    fit = -pof;
    return
end

for ii = 1:num_obj
    [u(:, ii), mse(:, ii)] = dace_predict(x, kriging_obj{ii});
    % [u(:, ii), mse(:, ii)] = predictor(x, kriging_obj{ii});
end
% mse = sqrt(max(0,mse));

% calcualate eim for objective
if num_obj == 1
    f = repmat(f, num_x, 1);
    imp = f - u;
    z = imp./mse;
    ei1 = imp .* Gaussian_CDF(z);
    ei1(mse==0)=0;
    ei2 = mse .* Gaussian_PDF(z);
    EIM = (ei1 + ei2);
else
    % for multiple objective problems
    % f - refers to pareto front
    r = 1.1*ones(1, num_obj);  % reference point
    num_pareto = size(f, 1);
    r_matrix = repmat(r,num_pareto,1);
    EIM = zeros(num_x,1);
    for ii = 1 : num_x
        u_matrix = repmat(u(ii,:),num_pareto,1);
        s_matrix = repmat(mse(ii,:),num_pareto,1);
        eim_single = (f - u_matrix).*Gaussian_CDF((f - u_matrix)./s_matrix) + s_matrix.*Gaussian_PDF((f - u_matrix)./s_matrix);
        EIM(ii) =  min(prod(r_matrix - f + eim_single,2) - prod(r_matrix - f,2));
    end
end

% calculate probability of feasibility for constraints
if nargin>3
    % the number of constraints
    num_con = length(kriging_con);
    % the kriging prediction and varince
    mu_g = zeros(size(x,1), num_con);
    mse_g = zeros(size(x,1), num_con);
    for ii = 1: num_con
        [mu_g(:, ii), mse_g(:, ii)] = dace_predict(x, kriging_con{ii});
        % [mu_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_con{ii}); %test
    end
    % mse_g = sqrt(max(0,mse_g));  %only compare with predictor not
    % dace_predict
    pof  = prod(Gaussian_CDF((0-mu_g)./mse_g), 2);
end
fit = -EIM .* pof;
end

function[krg_obj, krg_con] = surrogate_build(trainx, trainy, trainc)
% this function checks closeness of x
% use cross validation results to choose which elimination distance to uses

nd = 3;  %nd number of digits
n_obj  =  size(trainy,2);
n_con =  size(trainc,2);  %compatible with trainc=[]
krg_obj  =  cell(1, n_obj);
krg_con = cell(1, n_con);

for i=1:n_obj
    %decide on krg obj
    krg_obj{i} = cvtrain_withdatacheck(trainx, trainy(:, i), nd);
end

for i = 1:n_con
    %decide on krg con
    krg_con{i} = cvtrain_withdatacheck(trainx, trainy(:, i), nd);
end
end


function krg =  cvtrain_withdatacheck(trainx, trainy, nd)
% trainy is one column
if size(trainy, 2)>1
    error('response variable should have only one variable');
end

mse_crossnd = zeros(1, nd);
krg_crossnd = cell(1, nd);
k=5;
num_vari =  size(trainx,2);

for i = 1:nd % nd number of digits
    % --elimination scheme
    % --r potentially reduced data size
    [trainx_r,  trainy_r]  = close_elimination(trainx, trainy, i);
    
    % --create krg with cross validation
    cv_par = cvpartition(trainx_r, 'k', k);
    mse_perfolder = zeros(cv_par.NumTestSets,1);
    krg_perfolder  =  cell(1, k);
    for j = 1:k
        %---train
        trIdx = cv_par.training(i);
        teIdx = cv_par.test(i);
        krg_perfolder {j}= dacefit(trainx_r(trIdx , :), trainy_r(trIdx, ii),...
            'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
        %--test
        [~, mse_test]  = dace_predict(trainx_r(teIdx, :), krg_perfolder {j});
        mse_perfolder(j) =  sum(mse_test);
    end
    %-- pick krg
    [mse_min, ind] = min(mse_perfolder);
    mse_crossnd(i) = mse_min;
    krg_crossnd(i) = krg_perfolder{ind};
end
% -pick across nd
[~, ind] = min(mse_crossnd);
krg = krg_crossnd(ind);
end


function [x_after, y_after] = close_elimination(x, y, ndigi)
x = round(x, ndigi);
[x_after, ia] = unique(x, 'row');
y_after = y(ia, :);
end