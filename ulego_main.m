%%this script is to run ulego with multiple seeds
%
clearvars;
close all;
problem_folder = strcat(pwd,'\problems\SMD');
addpath(problem_folder);

problems = { 'smd9()', 'smd10()','smd11()','smd12()'};
seeds = linspace(1, 11, 11);

np = length(problems);
ns = length(seeds);

paras=cell(1, np*ns);
for i = 1:np
    for j = 1:ns
        paras{ (i-1)*ns +j } = {problems{i}, seeds(j)};
    end
end

nrun = length(paras);
parfor i = 1:nrun
    ulego(paras{i}{1}, paras{i}{2},'EIMnext_znorm' );
end