%% this test is single run of the main process, prepare for pararun
% this test is to confirm the process and 
% result of bilevel optimization
% claim is using ulego, the optimization should be as good as 
% reported in paper, start with smd9/10/11/12
% -----------------------
seed = 1;
rng(seed, 'twister');

workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
problem_folder = strcat(upperfolder,'\problems\BLTP');
addpath(problem_folder);
addpath(upperfolder);



prob = 'bltp1()';
ulego(prob, seed, 'EIMnext_znorm', 'EIM_eval');

rmpath(problem_folder); 
rmpath(upperfolder);