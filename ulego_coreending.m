function ulego_coreending(xu, fu, fc, xl, prob, seed, n_up, n_low, method)
% fu is converted to one from archive
[xu_best, fu, cu, ~, index] = localsolver_startselection(xu, fu, fc);
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

end