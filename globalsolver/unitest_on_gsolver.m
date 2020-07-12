%% development test on gsolver

workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
problem_folder = strcat(upperfolder,'\problems\SMD');
addpath(problem_folder);

prob = smd1();
xu = [0, 0];

%-global search
obj = @(x)hyobj(x, xu, prob);
con = @(x)hycon(x, xu, prob);

opts.popsize = 20;
opts.gen = 50;

gsolver(obj,prob.n_lvar, prob.xl_bl, prob.xl_bu, [], con, opts);

rmpath(problem_folder)

function [f] = hyobj(x, xu, prob)
    [f, ~] = prob.evaluate_l(xu, x);
end

function [c, ceq] = hycon(x, xu, prob)
    [~, c] = prob.evaluate_l(xu, x);
    ceq = [];
end