%% This unittest is for testing llmatch_sao
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

workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
problem_folder = strcat(upperfolder,'\problems\SMD');
method_folder = strcat(upperfolder, '\globalsolver');
addpath(problem_folder);
addpath(upperfolder);
addpath(method_folder);


% xu, prob, num_pop, num_gen, iter_freq

%----smd test---
problem = smd8(); 
xu = [0, 0];  %smd1,2, 3, 4, 5, 6, 7, 8
% problem = smd7(); xu = [0, 0];
% problem = smd6(); xu = [0, 0];
% problem = smd5(); xu = [0, 0];
% problem = smd4(); xu = [0, 0];
% problem = smd3(); xu = [0, 0];
% problem = smd2(); xu = [0, 0];
% problem = smd1(); xu = [0, 0];
% problem = smd9(); xu = [0, 0];   %smd9
% problem = smd10(); xu = [1, 1];  %smd10
% problem = smd11(); xu = [0, 0];  %smd11
% problem = smd12(); xu = [1, 1];  %smd12

rng(seed, 'twister');
[xl, n, flag] = llmatch_sao_archiveinsert(xu, problem, 20, 100, 20);
disp(xl);
disp(n);
disp(flag);
[f,c] = problem.evaluate_u(xu, xl);
disp(f);
acc = abs(f-problem.uopt);
disp(acc);


rmpath(problem_folder);
rmpath(upperfolder);
rmpath(method_folder);