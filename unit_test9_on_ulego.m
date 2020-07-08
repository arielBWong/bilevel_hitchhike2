%% this test is single run of the main process, prepare for pararun
% this test is to confirm the process and 
% result of bilevel optimization
% claim is using ulego, the optimization should be as good as 
% reported in paper, start with smd9/10/11/12
% -----------------------
seed = 1;
rng(seed, 'twister');
problem_folder = strcat(pwd,'\problems\BLTP');
addpath(problem_folder);

prob = 'bltp5()';
ulego(prob, seed, 'EIMnext_znorm');

rmpath(problem_folder); 