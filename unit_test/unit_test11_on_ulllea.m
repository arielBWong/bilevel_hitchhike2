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
addpath(problem_folder);
addpath(upperfolder);


prob = 'mobp5()';
ul_llea(prob, seed);

rmpath(problem_folder); 
rmpath(upperfolder);