function ulego_coreending(xu, fu, fc, xl, prob, seed, n_up, n_low, method)
% fu is converted to one from archive
[xu_best, fu, cu, ~, index] = out_select(xu,  xl, prob);

name = prob.name;
fprintf('returned index %d for problem %s seed %d\n', index, name, seed);

xl_best = xl(index, :);
[fl, cl] = prob.evaluate_l(xu_best, xl_best);

cu(cu < 1e-6) = 0;
cl(cl < 1e-6) = 0;
% check feasibility
num_conu = size(cu, 2);
num_conl = size(cl, 2);
cu = sum(cu<=0, 2)==num_conu;
cl = sum(cl<=0, 2)==num_conl;
perf_record(prob, fu, cu, fl, cl, n_up, n_low, seed, method);
archive_record(xu, xl, prob, method, seed);

end

function archive_record(xu, xl, prob, method, seed)
num = prob.n_lvar;
savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num), '_', method);

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