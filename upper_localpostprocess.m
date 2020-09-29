function upper_localpostprocess(prob, xu, xl, fu, n_feval, seed, method)
% this method serves as an addon 
% for core algorithm including two more steps
% (1) bilevel local search
% (2) re-evaluation
% input
% prob
% xu: all xu training data
% xl: all xl training data
% fu: all fu objectives
% n_feval: lower level number of function evaluation
%----------------------------------------

%-bilevel local search

[xu_start,~, ~, ~, ~] = out_select(xu,  xl, prob);
penaltyf =  max(fu, [], 1); % for lower problem being constraint
[newxu, newxl, n_up, n_low] = blsovler(prob, xu_start, 20, 20, 20, 10, penaltyf);
n_up = n_up + size(xu, 1);
n_low = n_low + n_feval;
xu_all = [xu; newxu];

%-final hybrid ll search
%-- use newxu
hy_pop = 20;
hy_gen = 20;

[newxl, feval, flag] = hybrid_llsearch(newxu, newxl, prob, hy_pop, hy_gen);

n_low = n_low + feval;
xl_all = [xl; newxl];

%-performance record
%--constraints compatible
[fu, cu] = prob.evaluate_u(newxu, newxl);
 
% scatter(newxu, fu, 'r'); drawnow;

[fl, cl] = prob.evaluate_l(newxu, newxl);
num_conu = size(cu, 2);
num_conl = size(cl, 2);

% contraint tolerance adjust
cu(cu < 1e-6) = 0;
cl(cl < 1e-6) = 0;
% check feasibility
cu = sum(cu<=0, 2)==num_conu;
cl = sum(cl<=0, 2)==num_conl;

perf_record(prob, fu, cu, fl, cl, n_up, n_low, seed, method, true);
archive_record(xu_all, xl_all, prob, method, seed);

end



function archive_record(xu, xl, prob, method, seed)
num = prob.n_lvar;
savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num), '_', method, '_addon');

n = exist(savepath);
if n ~= 7
    mkdir(savepath)
end

savename = strcat(savepath, '\xu_', num2str(seed),'.csv');
csvwrite(savename, xu);
savename = strcat(savepath, '\xl_', num2str(seed),'.csv');
csvwrite(savename, xl);

[fu, ~] = prob.evaluate_u(xu, xl);
savename = strcat(savepath, '\fu_', num2str(seed),'.csv');
csvwrite(savename, fu);


end