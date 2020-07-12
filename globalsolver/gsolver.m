function [bestx, bestf, bestc] = gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param)
% This function is a wrapper on methods in global optimization/minimization
% folder. Main process is nsga, but reproducation is a DE operator
% Be aware, this method only handle inequality constraints 
% Don't use on equality constraints
% input:
%       funh_obj                                                : function handle to objective function
%       num_xvar                                              : number design variables
%       lb                                                             : upper bound of design varibles
%                                                                                      1d row array
%       up                                                           : lower bound of design variables
%       initmatrix                                               :  partial population to be embeded in
%                                                                                      initial population 
%       funh_con                                               : function handle to constraint functions
%       opts                                                        : structure specifying ea parameters
% output:
%       bestx                                                      : global search results of design variables  (best value or nd front)                      
%       bestf                                                       : global search results of objective values   (best value or nd front)      
%       bestc                                                      : global search results of constraints
%                                                                                       for constraint problems, if no feasible is found, return least infeasible one
%  ** under development of archive handling
%--------------------------------------------------------------------------

bestx =NaN;
bestf = NaN;
bestc = NaN;

% Initialization
[pop,archive] = initialize_pop(funh_obj, funh_con, num_xvar, lb, ub, initmatrix, param, []);

gen=1;
while gen<= param.gen
    % Recombination
    child.X=generate_child(pop,param,archive,prob,gen,objnum);      
    % Evaluate and Order
    [pop,archive]= evaluate_order(pop,prob,param,archive,evals,objnum,child.X,gen,track);     
    % Reduce 2N to N
     [pop]=reduce_pop(pop,param);
     
     gen = gen+1;
    % disp(gen);
end

num_obj = size(pop.F, 2);
num_con = size(pop.C, 2);

if num_obj == 1
    bestx = pop.X(1, :);
    bestf = pop.F(1, :);
    bestc = pop.C(1, :);
    return 
end

% return feasible nd front
if num_obj > 1
    if ~isempty(pop.C)                                    % constraint problem
        fy_ind = sum(pop.C, 2) ==num_con;
    else
        fy_ind = [1:para.popsize];                    % unconstraint problem
    end
end

% feasible exists for mo problem
if sum(fy_ind) > 0
    [fronts, ~, ~] = nd_sort(pop.F, find(fy_ind));
    bestf = pop.F(fronts(1).f, :);
    bestx = pop.X(fronts(1).f, :);
    bestc = pop.C(fronts(1).f, :);
else
    % no feasible solution in final population
    % return least infeasible on
    sumc= sum(pop.C, 2);
    [~, id] = sort(sumc);
    bestf = pop.F(id(1),  :);
    bestx = pop.X(id(1),  :);
    bestc = pop.C(id(1),  :);  
end 
end