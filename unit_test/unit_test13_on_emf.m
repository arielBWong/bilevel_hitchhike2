%% this script run to generate empirical pareto front
% for bilevel problems use ea-ea to speed up as surrogate takes time with
% number of training data increases


clearvars;
close all;



workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);


problem_folder = strcat(upperfolder,'\problems\MOBP');
solver_folder = strcat(upperfolder, '\globalsolver');
sort_folder = strcat(upperfolder,'\ND_Sort');
addpath(problem_folder);
addpath(solver_folder);
addpath(sort_folder);
problem_folder = strcat(upperfolder,'\problems\TP');
addpath(problem_folder);
problem_folder = strcat(upperfolder,'\problems\DS');
addpath(problem_folder);
addpath(upperfolder);


% problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
problems = { 'tp1()' ,'tp2(6)' ,'tp3()' ,'tp4()' , 'ds1(6)', 'ds2(6)', 'ds3(6)', 'ds4(3,2)', 'ds5(3, 2)'};
seeds = linspace(1, 1, 1);
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


