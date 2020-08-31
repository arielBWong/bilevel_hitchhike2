%% this test is single run of the main process,
% this test is to confirm the process and
% result of bilevel optimization on upper level moc lower level soc
% round 1: to be passed
% -----------------------
clearvars;
close all;
seed = 1;
rng(seed, 'twister');


workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
problem_folder = strcat(upperfolder,'\problems\MOBP');
sort_folder = strcat(upperfolder,'\ND_Sort');
addpath(problem_folder);
addpath(upperfolder);
addpath(sort_folder);
problem_folder = strcat(upperfolder,'\problems\DSM');
addpath(problem_folder);
solver_folder = strcat(upperfolder,'\globalsolver');
addpath(solver_folder);

% problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
prob ='dsm1(2)';
ul_llea(prob, seed);

prob = eval(prob);
savepath = strcat(pwd, '\result_folder\', prob.name, '_ea_ea' );
savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv');
nd_front = csvread(savename_fu);
scatter(nd_front(:,1), nd_front(:,2)); drawnow;

rmpath(problem_folder);
rmpath(upperfolder);
rmpath(sort_folder);