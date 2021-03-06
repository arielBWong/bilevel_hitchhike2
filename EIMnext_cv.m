function[best_x, info] = EIMnext_cv(train_x, train_y, xu_bound, xl_bound, ...
    num_pop, num_gen, train_c)
% method of using EIM to generate next point with crossvalidation on krg
% normalization is zscore on all variables
% usage:
%
% input:        train_x                  - design variables
%                                                           1/2d array: (num_samples, num_varibles)
%                    train_y                 - objective values
%                                                           1/2d array: (num_samples, num_objectives)
%                    xu_bound           - upper bound of train_x
%                                                           1d array
%                    xl_bound             - lower bound of train_x
%                                                           1d array
%                    num_pop            - EIM optimization parameter
%                    num_gen             - EIM optimization parameter
%                    train_c                   - constraints values
%                                                            1/2d array: (num_samples, num_constraints)
% output:     best_x                   - proposed next x to be evaluated by EIM
%                    info                       - returned information for functor caller to recreate
%                                                    - or check information
%                                                   - info.krg
%                                                   - info.krgc
%                                                   - info.train_xmean
%                                                    - info.train_ymean
%                                                   - info.train_xstd
%                                                   - info.train_ystd
%                                                    - info.info.eim_normf
%--------------------------------------------------------------------------

% number of objective
num_obj = size(train_y, 2);
kriging_obj = cell(1,num_obj);

% number of input designs
num_x = size(train_x,1);
num_vari = size(train_x, 2);

% normalise train data x
[train_x_norm, x_mean, x_std] = zscore(train_x, 1, 1); 

% adjust x bounds according to zscore normalization
ub = (xu_bound - x_mean)./ x_std;
lb = (xl_bound - x_mean)./ x_std;

%  normalise train data y 
if num_obj > 1
    [train_y_norm, y_mean, y_std] = zscore(train_y, 1, 1);
    % further normalization ??
    train_y_norm = (train_y_norm -repmat(min(train_y_norm),num_x,1))./...
        repmat(max(train_y_norm)-min(train_y_norm),num_x,1);
else
    [train_y_norm, y_mean, y_std] = zscore(train_y, 1, 1);
end

% prepare f_best for EIM, first only consider non-constraint situation
% for contrainted problems update is done later
if num_obj > 1
    index_p = Paretoset(train_y_norm);
    f_best = train_y_norm(index_p, :); % constraint problem has further update
else
    f_best = min(train_y_norm, [], 1);
end

% compatibility with constraint problems
if ~isempty(train_c)
    % constraints should not be normalised
    % lazy in changing names
    train_c_norm = train_c;
    num_con = size(train_c, 2);
    % update f_best according to feasibility
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
    % Generate krging model  before dive into EIM
    [kriging_obj, kriging_con] = surrogate_build(train_x_norm, train_y_norm, train_c_norm);
    fitness_val = @(x)EIM_eval(x,f_best, kriging_obj, kriging_con);
else
    [kriging_obj, ~] = surrogate_build(train_x_norm, train_y_norm, []);
    fitness_val = @(x)EIM_eval(x, f_best, kriging_obj);
end

% call DE evolution
[best_x, eim_f] = DE(fitness_val, num_vari, lb, ub, num_pop, num_gen);

% convert best_x to denormalized value
best_x = best_x .* x_std + x_mean;  

% fix bound violation
best_x = fixbound_violation(best_x, xu_bound, xl_bound);

% form output 
info = struct();
info.eim_normf = eim_f;
info.krg = kriging_obj;
info.train_xmean = x_mean;
info.train_xstd = x_std;
info.train_ymean = y_mean;
info.train_ystd = y_std;
info.krg_con = [];
if  ~isempty(train_c)
    info.krgc = kriging_con;
end
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
% input: 
%        x                            - pop to evaluate
%        f                             - best f so far/feasible pareto front
%                                           in multi-objective probs
%        kriging_obj         - kriging model for objectives
%        kriging_con        - kriging model for constraints
% output: 
%       fit                           - pop fitness
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

%-----test  needed---
function[krg_obj, krg_con] = surrogate_build(trainx, trainy, trainc)
% this function checks closeness of x
% and use cross validation results to choose 
% which elimination distance to uses
% usage:
% input
%           trainx
%           trainy
%           trainc
% output
%           krg_obj
%           krg_con
%-----------------------------------------------------------------

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

