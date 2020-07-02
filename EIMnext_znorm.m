function[best_x, info] = EIMnext_znorm(train_x, train_y, xu_bound, xl_bound, ...
    num_pop, num_gen, train_c)
% method of using EIM to generate next point
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
%                   - info.train_cmean
%                   - info.train_xstd
%                   - info.train_ystd
%                   - info.train_cstd
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
%     train_y_norm = (train_y -repmat(min(train_y),num_x,1))./...
%          repmat(max(train_y)-min(train_y),num_x,1);
%     y_mean = NaN;
%     y_std = NaN;
    [train_y_norm, y_mean, y_std] = zscore(train_y, 1, 1);
    % further normalization ??
    train_y_norm = (train_y_norm -repmat(min(train_y_norm),num_x,1))./...
         repmat(max(train_y_norm)-min(train_y_norm),num_x,1);
else
    [train_y_norm, y_mean, y_std] = zscore(train_y, 1, 1);
end

for ii = 1:num_obj
    % kriging_obj{ii} = dace_train(train_x_norm,train_y_norm(:,ii));
    kriging_obj{ii} = dacefit(train_x_norm,train_y_norm(:,ii),...
        'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  % for test

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
    num_con = size(train_c, 2);
    kriging_con = cell(1,num_con);

    [train_c_norm, c_mean, c_std] = zscore(train_c, 1, 1);
    % version did not scale train_c
    % train_c_norm = train_c;
    % c_mean = NaN;
    % c_std = NaN;
    for ii = 1:num_con
        % kriging_con{ii} = dace_train(train_x_norm,train_c_norm(:,ii));
        kriging_con{ii} = dacefit(train_x_norm,train_c_norm(:,ii),...
             'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  %test

    end
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
    fitness_val = @(x)EIM_eval(x,f_best, kriging_obj, kriging_con);
else
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

%if nargin>6
if ~isempty(train_c)
    info.train_cmean = c_mean;
    info.train_cstd = c_std;
end

end

function x = fixbound_violation(x, xu, xl)
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
    mse_g = sqrt(max(0,mse_g));
    pof  = prod(Gaussian_CDF((0-mu_g)./mse_g), 2);
end
fit = -EIM .* pof;
end
