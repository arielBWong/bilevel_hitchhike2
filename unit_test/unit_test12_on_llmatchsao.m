%% This unittest is for testing llmatch_sao
% test through all smd problems
% set xu to optimum, llmatch should find xl optimum
% ** Conclusion: pass, but smd3 needs more init size, and iteration number
% smd3 parameter parameter 100(pop_size), 300(pop_gen), 100(feq), this
% problem is affected by init training data size 
clearvars;
close all;
%(1)create test problem
%(2)give xu optimum
%(3)print xl results

seed = 1;
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
% problem = smd8();  xu = [0, 0];  %smd1,2, 3, 4, 5, 6, 7, 8
%  problem = smd7(); xu = [0, 0];
% problem = smd6(); xu = [0, 0];
% problem = smd5(); xu = [0, 0];
 % problem = smd4(); xu = [0, 0];
problem = smd3(); xu = [0, 0]; % parameter 100(pop_size), 300(pop_gen), 100(feq)
% problem = smd2(); xu = [0, 0];
% problem = smd1(); xu = [0, 0];
% problem = smd9(); xu = [0, 0];   %smd9
% problem = smd10(); xu = [1, 1];  %smd10
% problem = smd11(); xu = [0, 0];  %smd11
% problem = smd12(); xu = [1, 1];  %smd12

rng(seed, 'twister');
[xl, n, flag] = llmatch_sao_archiveinsert(xu, problem, 	100, 200, 100);
% [xl, n, flag] = llmatch_sao_pop(xu, problem,100,200, 100);
% [xl, n, flag] = llmatch_sao_pop(xu, problem, 20, 400, 40);
disp(xl);
fprintf('number of function evaluation: %d\n', n);
fprintf('matching validation %d\n', flag);
[f,c] = problem.evaluate_u(xu, xl);
fprintf('upper level function value: %f\n', f);
acc = abs(f-problem.uopt);
fprintf('accuracy: %f\n', acc);


rmpath(problem_folder);
rmpath(upperfolder);
rmpath(method_folder);