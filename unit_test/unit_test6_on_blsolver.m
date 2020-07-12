%%
% this unit test is for blsolver

workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
problem_folder = strcat(upperfolder,'\problems\SMD');
addpath(problem_folder);
addpath(upperfolder);

addpath(problem_folder);

tic
prob = smd10();
xu_start = [-1.2136516906077681, -1.228916217819073];
inisize_l = 30;
numiter_l = 20;
[xu_end, xl_end, n_up, n_low] = blsovler(prob, xu_start, 100, 100, inisize_l, numiter_l);


toc

[fu,cu] = prob.evaluate_u(xu_end, xl_end);
[fl,cl] = prob.evaluate_l(xu_end, xl_end);
disp(fu);
disp(fl);

rmpath(problem_folder);
rmpath(upperfolder);