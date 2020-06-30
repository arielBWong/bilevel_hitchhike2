%%
% This unittest is for testing llmatch
% no past yet

clearvars;
close all;
%(1)create test problem
%(2)give xu optimum
%(3)print xl results

problem = smd9();
xu = [0, 0];
[xl, n, flag] = llmatch(xu, problem, 100, 100, 30, 20);
xl
