%% This unittest is for testing llmatch
% test through all smd problems
% set xu to optimum, llmatch should find xl optimum
% ** Conclusion: pass, but smd3 needs more init size, and iteration number

clearvars;
close all;
%(1)create test problem
%(2)give xu optimum
%(3)print xl results

seed = 2;
rng(seed, 'twister');
problem_folder = strcat(pwd,'\problems\SMD');
addpath(problem_folder);

initsize = 20;
numiter = 30;
problem = smd3();
xu = [0, 0];  %smd1,2, 3, 4, 5, 6, 7, 8
%xu = [0, 0];  %smd9
%xu = [1, 1];  %smd10
%xu = [0, 0];  %smd11
%xu = [1, 1];  %smd12

rng(seed, 'twister');
[xl, n, flag] = llmatch(xu, problem, 100, 100, initsize, numiter);
disp(xl);
[f,c] = problem.evaluate_u(xu, xl);
disp(f);
acc = abs(f-problem.uopt);
disp(acc);


rmpath(problem_folder);
