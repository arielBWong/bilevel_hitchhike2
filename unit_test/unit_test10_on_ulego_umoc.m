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

problem_folder = strcat(upperfolder,'\problems\DSM');
addpath(problem_folder);

addpath(problem_folder);
addpath(upperfolder);
addpath(solver_folder);
addpath(sort_folder);


tic;
prob = 'dsm1(2)';
ulego_umoc(prob, seed, 'EIMnext_znorm', 'EIM_eval', 'normalization_nd', 'EIMnext_znorm');
% ulego_umoc(prob, seed, 'EIMnext_znorm', 'Ehv_eval', 'normalization_nd', 'Believer_next');


prob = eval(prob);
savepath = strcat(pwd, '\result_folder\', prob.name, '_EIM_eval' );
savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv');
nd_front = csvread(savename_fu);
scatter(nd_front(:,1), nd_front(:,2)); drawnow;


toc
rmpath(problem_folder); 
rmpath(upperfolder);
rmpath(solver_folder);
rmpath(sort_folder);