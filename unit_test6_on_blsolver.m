%%
% this unit test is for blsolver

problem_folder = strcat(pwd,'\problems\SMD');
addpath(problem_folder);

tic
prob = smd10();
xu_start = [-1.2136516906077681, -1.228916217819073];
inisize_l = 30;
numiter_l = 20;
[xu_end] = blsovler(prob, xu_start, 100, 100, inisize_l, numiter_l);
toc

rmpath(problem_folder);