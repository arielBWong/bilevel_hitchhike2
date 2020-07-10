function [bestx, bestf, bestc] = gsolver(obj, num_xvar, lb, ub, con, ops)
% This function is a wrapper on methods in global optimization/minimization
% folder. Main process is nsga, but reproducation is DE operators
% Be aware, this method only handle inequality constraints 
% Don't use on equality constraints
% input:
%       obj                               : function handle to objective function
%       num_xvar                          : number design variables
%       lb                                : upper bound of design varibles
%                                              1d row array
%       up                                : lower bound of design variables
%       con                               : function handle to constraint functions
%       ops                               : structure specifying ea parameters
% output:
%       bestx
%       bestf
%       bestc
%--------------------------------------------------------------------------


% Initialization
[pop,archive,evals,track]=Initialize_pop(prob,param,evals,objnum);

gen=1;
while gen<= param.gen
    % Recombination
    child.X=Generate_child(pop,param,archive,prob,gen,objnum);
    
    % Evaluate and Order
    [pop,archive,evals,track]= Evaluate_order(pop,prob,param,archive,evals,objnum,child.X,gen,track);
    
    % Reduce 2N to N
    [pop,archive]=Reduce(pop,param,archive,gen,prob);
    
    gen = gen+1;
    disp(gen);
end

end