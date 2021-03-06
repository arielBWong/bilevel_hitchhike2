%%
f = readtable('nd_front_scaled.txt')

[num_pareto,num_obj] = size(f);
% number of input designs
num_x = size(x,1);
r = 1.1*ones(1, num_obj);

r_matrix = repmat(r,num_pareto,1);
for ii = 1 : num_x
    u_matrix = repmat(u(ii,:),num_pareto,1);
    s_matrix = repmat(s(ii,:),num_pareto,1);    
    EIM = (f - u_matrix).*Gaussian_CDF((f - u_matrix)./s_matrix) + s_matrix.*Gaussian_PDF((f - u_matrix)./s_matrix);
    y(ii) =  min(prod(r_matrix - f + EIM,2) - prod(r_matrix - f,2));
end
%-----------------------------------------------------
% the objective is maximized
obj = -y;