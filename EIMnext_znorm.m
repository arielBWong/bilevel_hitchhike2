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
[train_x_norm, x_mean, x_std] = zscore(train_x, 0, 1);  % flag 0 sample zscore

% adjust x bounds according to zscore normalization
ub = (xu_bound - x_mean)./ x_std;
lb = (xl_bound - x_mean)./ x_std;

% train objective kriging model
if num_obj > 1
    [train_y_norm, y_mean, y_std] = zscore(train_y, 0, 1);
    % further normalization ??
     train_y_norm = (train_y_norm -repmat(min(train_y_norm),num_x,1))./...
         repmat(max(train_y_norm)-min(train_y_norm),num_x,1);
else
    [train_y_norm, y_mean, y_std] = zscore(train_y, 0, 1);
end

for ii = 1:num_obj
    kriging_obj{ii} = dace_train(train_x_norm,train_y_norm(:,ii));
end

% prepare f_best for EIM, first only consider non-constraint situation
if num_obj > 1
    index_p = Paretoset(train_y_norm);
    f_best = train_y_norm(index_p, :); % constraint problem has further process
else
    f_best = min(train_y_norm, [], 1);
end

% compatibility with constraint problems
if nargin>6
    num_con = size(train_c, 2);
    kriging_con = cell(1,num_con);
    [train_c_norm, c_mean, c_std] = zscore(train_c, 0, 1);
    train_c_norm = train_c;
    c_mean = NaN;
    c_std = NaN;
    for ii = 1:num_con
        kriging_con{ii} = dace_train(train_x_norm,train_c_norm(:,ii));
    end
    % adjust f_best according to feasibility
    % feasibility needs to be valued in original range
    index_c = sum(train_c <= 0, 2) == num_con;
    if sum(index_c) == 0 % no feasible
        f_best = [];
    else
        feasible_trainy_norm = train_y_norm(index_c, :);
        % still needs nd front
        index_p = Paretoset(feasible_trainy_norm);
        f_best = feasible_trainy_norm(index_p, :);
    end
    fitness_val = @(x)EIM_eval(x,f_best, kriging_obj, kriging_con);
else
    fitness_val = @(x)EIM_eval(x, f_best, kriging_obj);
end

% call DE evolution
[best_x, eim_f] = DE(fitness_val, num_vari,lb, ub, num_pop, num_gen);

% convert best_x to denormalized value
best_x = best_x .* x_std + x_mean;  % commit for test

info = struct();
info.eim_normf = eim_f;
info.krg = kriging_obj;
info.train_xmean = x_mean;
info.train_xstd = x_std;
info.train_ymean = y_mean;
info.train_ystd = y_std;

if nargin>6
    info.train_cmean = c_mean;
    info.train_cstd = c_std;
end

end


