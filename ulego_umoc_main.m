%%this script is to run ulego with multiple seeds
%
clearvars;
close all;
problem_folder = strcat(pwd,'\problems\MOBP');
addpath(problem_folder);
problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
algs = {'EIM_eval', 'Ehv_eval'};

seeds = linspace(1, 1, 1);

np = length(problems);
ns = length(seeds);
na = length(algs);


paras=cell(1, np*ns * na);
nn = 1;
for i = 1:np
    for j = 1:ns
        for k = 1:na
            paras{ nn } = {problems{i}, seeds(j), algs{k}};
            nn = nn + 1;
        end
    end
end


nrun = length(paras);
parfor i = 1:nrun
    ulego_umoc(paras{i}{1}, paras{i}{2},'EIMnext_znorm' , paras{i}{3});
end

rmpath(problem_folder);