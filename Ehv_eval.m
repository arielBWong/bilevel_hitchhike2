function [fit] = Ehv_eval(x, f, kriging_obj, kriging_con)
% function of using expected hyper volume as fitness evaluation
% usage:
%
% input:
%        x                          - pop to evaluate
%        f                           - best f so far/feasible pareto front
%                                           in multi-objective probs
%        kriging_obj       - kriging model for objectives
%        kriging_con      - kriging model for constraints
% output:
%       fit                         - pop fitness
%--------------------------------------------------------------------------

% number of input designs/population size
num_x = size(x,1);
num_obj = size(f, 2);
if num_obj == 1
    error('this evaluation is not for single objective problem');
end
ref_point = ones(1, num_obj) * 1.1;

% the kriging prediction and varince
u = zeros(num_x,num_obj);
mse = zeros(num_x,num_obj);

% expected f value from kriging model
for ii = 1:num_obj
    [u(:, ii), mse(:, ii)] = dace_predict(x, kriging_obj{ii});
end

% constraint problem handling
if nargin > 3
    num_con = length(kriging);
    uc = zeros(num_x, num_con);
    for ii = 1:num_con
        [uc(ii, :), ~] = dace_predict(x, kriging_con{ii});
    end
    uc(uc<=0) = 0; % prepare for feasiblity check
end

basehv = Hypervolume(f, ref_point);
fit = zeros(num_x, 1);                                   % contribution of each x-f to existing nd front

for ii = 1:num_x                                              %population size
    extendf = [f; u(ii, :)];
    % put constraint into consideration
    % skip computating if  infeasible
    if nargin > 3 && sum(uc(ii, :), 2) > 0 
        fit(ii) = 0;
        continue;
    end
    extendedhv = Hypervolume(extendf, ref_point);
    fit(ii) = -(extendedhv - basehv);
end

end