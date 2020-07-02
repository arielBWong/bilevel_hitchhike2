%%
% This unittest is for testing llmatch
% no past yet

clearvars;
close all;
%(1)create test problem
%(2)give xu optimum
%(3)print xl results

seed = 1;
problem_folder = strcat(pwd,'\problems\SMD');
addpath(problem_folder);

problem = smd9();
xu = [0, 0];
rng(seed, 'twister');
[xl, n, flag] = llmatch(xu, problem, 100, 100, 30, 20);



rmpath(problem_folder);
