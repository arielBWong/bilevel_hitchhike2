function [fit] = krgbeliever_eval(x, f, kriging_obj, kriging_con)
% function of using expected hyper volume as fitness evaluation
% usage:
%
% input:
%        x - pop to evaluate
%        f - best f so far/feasible pareto front
%                                           in multi-objective probs
%        kriging_obj - kriging model for objectives
%        kriging_con - kriging model for constraints
% output:
%        fit - pop fitness
%--------------------------------------------------------------------------

% number of input designs
num_x = size(x, 1);
num_obj = size(f, 2);

% the kriging prediction and varince
u = zeros(num_x,num_obj);
for ii = 1:num_obj
    [u(:, ii), ~] = dace_predict(x, kriging_obj{ii});
end

% construct constraints
if isempty(f) && nargin > 3   % refer to no feasible solution
    % only calculate probability of feasibility
    num_con = length(kriging_con);
    % the kriging prediction and varince
    mu_g = zeros(size(x,1), num_con);
    for ii = 1: num_con
        [mu_g(:, ii), ~] = dace_predict(x, kriging_con{ii});
    end
end

% sort according to fitness and constraints


