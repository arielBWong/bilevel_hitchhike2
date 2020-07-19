%% this script is to run ul_llea with multiple seeds
% 

clearvars;
close all;
problem_folder = strcat(pwd,'\problems\MOBP');
addpath(problem_folder);

problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
seeds = linspace(1, 15, 15);

np = length(problems);
ns = length(seeds);
paras=cell(1, np*ns);

nn = 1;
for i = 1:np
    for j = 1:ns
        paras{nn} = {problems{i}, seeds(j)};
    end
end

nrun = length(paras);
parfor i = 1:nrun
    ul_llea(paras{i}{1}, paras{i}{2});
end
