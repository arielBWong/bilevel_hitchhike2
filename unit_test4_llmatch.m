%%
% This unittest is for testing llmatch
% no past yet

clearvars;
close all;
%(1)create test problem
%(2)give xu optimum
%(3)print xl results

seed = 2;
rng(seed, 'twister');
problem_folder = strcat(pwd,'\problems\SMD');
addpath(problem_folder);

problem = smd3();
xu = [0, 0];  %smd1,2, 3, 4, 5, 6, 7, 8
%xu = [0, 0];  %smd9
%xu = [1, 1];  %smd10
%xu = [0, 0];  %smd11
%xu = [1, 1];  %smd12

rng(seed, 'twister');
[xl, n, flag] = llmatch(xu, problem, 100, 100, 30, 20);
xl
[f,c] = problem.evaluate_u(xu, xl);
f
acc = abs(f-problem.uopt);
disp(acc)


rmpath(problem_folder);
