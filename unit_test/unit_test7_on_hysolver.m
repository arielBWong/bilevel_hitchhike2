%% Use this to check with SMD3
% one chance
workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
problem_folder = strcat(upperfolder,'\problems\SMD');
addpath(problem_folder);
addpath(upperfolder);

seed = 2;

hy_pop = 20;
hy_gen = 50;

prob = smd3();
xu = [0, 0];
%xl_start = [10, 10, pi/3];
xl_start = [];

[newxl, n_feval, flag] = hybrid_llsearch(xu, xl_start, prob, hy_pop, hy_gen);

[f, c] = prob.evaluate_u(xu, newxl);
disp(f);

rmpath(problem_folder);
rmpath(upperfolder);

