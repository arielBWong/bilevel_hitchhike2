function obj = Infill_Standard_EIM_Maximin(x, kriging_obj, non_dominated_front)
%-----------------------------------------------------
% [1]  D. Zhan, Y. Cheng, J. Liu, Expected Improvement Matrix-based Infill
% Criteria for Expensive Multiobjective Optimization, IEEE Transactions
% on Evolutionary Computation, DOI: 10.1109/TEVC.2017.2697503
%-----------------------------------------------------
% the input parameters
f = non_dominated_front;
% number of non-dominated points
% number of objectives
[num_pareto,num_obj] = size(f);
% number of input designs
num_x = size(x,1);
y = zeros(num_x,1);
%-----------------------------------------------------
% the kriging prediction and varince
u = zeros(num_x,num_obj);
mse = zeros(num_x,num_obj);
for ii = 1:num_obj
    [u(:, ii),mse(:, ii)] = predictor(x, kriging_obj{ii});
end
s=sqrt(max(0,mse));
%-----------------------------------------------------
for ii = 1 : num_x
    u_matrix = repmat(u(ii,:),num_pareto,1);
    s_matrix = repmat(s(ii,:),num_pareto,1);
    EIM = (f - u_matrix).*Gaussian_CDF((f - u_matrix)./s_matrix) + s_matrix.*Gaussian_PDF((f - u_matrix)./s_matrix);
    y(ii) = min(max(EIM,[],2));
end
%-----------------------------------------------------
% the objective is maximized
obj = -y;
end




