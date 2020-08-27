%% this script run to generate empirical pareto front
% for bilevel problems use ea-ea to speed up as surrogate takes time with
% number of training data increases


clearvars;
close all;
problem_folder = strcat(pwd,'\problems\MOBP');
solver_folder = strcat(pwd, '\globalsolver');
sort_folder = strcat(pwd,'\ND_Sort');
addpath(problem_folder);
addpath(solver_folder);
addpath(sort_folder);
problem_folder = strcat(pwd,'\problems\TP');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\DS');
addpath(problem_folder);


% problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
problems = { 'tp1()' ,'tp2(6)' ,'tp3()' ,'tp4()' , 'ds1(6)', 'ds2(6)', 'ds3(6)', 'ds4(3,2)', 'ds5(3, 2)'};
seeds = linspace(1, 5, 5);
np = length(problems);
ns = length(seeds);

paras = cell(1, np*ns);
nn = 1;
for jj = 1:ns
    for ii = 1:np
        problem = problems{ii};
        seed = seeds(jj);
        paras{nn} = {problem, seed};
        nn = nn + 1;
    end
end


for ii = 1: np*ns
    empirical_pfc(paras{ii}{1}, paras{ii}{2});
end


