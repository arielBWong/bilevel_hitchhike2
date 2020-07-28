function [best_x, info] = Believer_next(train_x, train_y, xu_bound, xl_bound, ...
    num_pop, num_gen, train_c, fitnesshn, normhn)
% This function is also a surrogate assisted method, but it does not
% consider variances in prediction
% normalization is zscore on all variables
% usage:
%
% input:
%        train_x                - design variables
%                                           1/2d array: (num_samples, num_varibles)
%        train_y                - objective values
%                                           1/2d array: (num_samples, num_objectives)
%        xu_bound          - upper bound of train_x
%                                           1d array
%        xl_bound           - lower bound of train_x
%                                           1d array
%        num_pop          - EIM optimization parameter
%        num_gen          - EIM optimization parameter
%        train_c                - constraints values
%                                           1/2d array: (num_samples, num_constraints)
%        fitnesshn           -fitness valuation handle (callback) for ea process of
%                                           proposing next poing
%         normhn             -normalization handle on y
% output:
%       best_x                 - proposed next x to be evaluated by EIM
%       info                     - returned information for functor caller to recreate
%                                   - or check information
%                                   - info.krg
%                                   - info.krgc
%                                   - info.train_xmean
%                                   - info.train_ymean
%                                   - info.train_xstd
%                                   - info.train_ystd
%                                   - info.info.eim_normf
%--------------------------------------------------------------------------

% number of objective
info = struct();
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
    info.train_ynormmin = min(train_y_norm);
    info.train_ynormmax = max(train_y_norm);
    % further normalization
    train_y_norm = normhn(train_y_norm);
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
    % constraints should not be normalised
    % version did not scale train_c
    train_c_norm = train_c;
    for ii = 1:num_con
        % kriging_con{ii} = dace_train(train_x_norm,train_c_norm(:,ii));
        kriging_con{ii} = dacefit(train_x_norm, train_c_norm(:,ii),...
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
    fitness_val = @(x)fitnesshn(x,f_best, kriging_obj, kriging_con);
else
    fitness_val = @(x)fitnesshn(x, f_best, kriging_obj);
end

% call global solver evolution
funh_obj = @(x)obj(x,kriging_obj);
if ~isempty(train_c)
    funh_con = @(x)con(x,kriging_con);
else
    funh_con = @(x)con(x, []);
end
param = struct();
param.gen= num_gen;
param.popsize = num_pop;
[best_x, best_f, ~, ~] = gsolver(funh_obj, num_vari, lb, ub, [], funh_con, param);

% convert best_x to denormalized value
best_x = best_x .* x_std + x_mean;  % commit for test

% fix bound violation
best_x = fixbound_violation(best_x, xu_bound, xl_bound);

%--info
info.eim_normf = best_f;
info.krg = kriging_obj;
if  ~isempty(train_c)
    info.krgc = kriging_con;
end
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

function [c] = con(x, krg)

if ~isempty(krg)
    % number of input designs
    num_x = size(x, 1);
    num_con = length(krg);
    
    % the kriging prediction and varince
    c = zeros(num_x,num_con);
    for ii = 1:num_con
        [c(:, ii), ~] = dace_predict(x, krg{ii});
    end
else
    c = [];
end
end


function [f] = obj(x, krg)
% number of input designs
num_x = size(x, 1);
num_obj = length(krg);

% the kriging prediction and varince
f = zeros(num_x,num_obj);
for ii = 1:num_obj
    [f(:, ii), ~] = dace_predict(x, krg{ii});
end
end


