 %% this test is single run of the main process, 
% this test is to confirm the process and 
% result of bilevel optimization on upper level moc lower level soc
% round 1: to be passed
% -----------------------
seed = 1;
rng(seed, 'twister');


workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
problem_folder = strcat(upperfolder,'\problems\MOBP');
solver_folder = strcat(upperfolder,'\globalsolver');
sort_folder = strcat(upperfolder,'\ND_Sort');

addpath(problem_folder);
addpath(upperfolder);
addpath(solver_folder);
addpath(sort_folder);

tic;
prob = 'mobp7()';
% ulego_sao_pop(prob, seed, 'normalization_nd');
% ulego_sao_archiveinsert(prob, seed, 'normalization_nd');
ulego_sao(prob, seed, 'normalization_nd');

toc
rmpath(problem_folder); 
rmpath(upperfolder);
rmpath(solver_folder);
rmpath(sort_folder);