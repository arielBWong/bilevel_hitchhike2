%% reproduce upper level plots
% read xu and calculate best matching xl
% recreate plot
% read xu and fu
% recreate plot
% 
problems = {'smd5()'};
algs = {'hyb'}; 
seed = 1;

prob = eval(problems{1});

xu_savefile = strcat(pwd, '\result_folder\', prob.name,'_', num2str(prob.n_lvar), '_',  algs{1});
xu_savefile = strcat(xu_savefile, '\xu_', num2str(seed), '.csv' );

xu = csvread(xu_savefile);
xl = [];
for ii = 1: size(xu, 1)
    xl_prime = matching_prime(xu(i, :), prob);
    xl = [xl; xl_prime];
end

[fu, fc] = prob.evaluate_u(xu, xl);
[krg_obj, krg_con, ~] = update_surrogate(xu,fu, fc, 'normalization_z' );






function xl_prime =  matching_prime(xu, prob)
xl1 = ones(1, prob.q);
xu2 = xu(1, prob.p + 1:end);
xl2 = sqrt(abs(xu2));
xl_prime = [xl1, xl2];
end