clearvars;
close all;


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



prob = 'tp2(3)';
probc = eval(prob);
pf = probc.upper_pf(100);



% pop
ulego_sao_pop(prob, seed, 'normalization_nd');
savepath = strcat(pwd, '\result_folder\', probc.name, '_sao_popinsert' );
savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv');
nd_frontp = csvread(savename_fu);


% archive one by one
ulego_sao_archiveinsert(prob, seed, 'normalization_nd');
savepath = strcat(pwd, '\result_folder\', probc.name, '_sao_archiveinsert' );
savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv');
nd_fronta = csvread(savename_fu);




% eim
ulego_umoc(prob, seed, 'EIMnext_znorm', 'EIM_eval', 'normalization_nd', 'EIMnext_znorm');
savepath = strcat(pwd, '\result_folder\', probc.name, '_EIM_eval' );
savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv');
nd_front = csvread(savename_fu);

scatter(pf(:,1), pf(:,2),'bo');  hold on;
scatter(nd_frontp(:,1), nd_frontp(:,2), '+');
scatter(nd_front(:,1), nd_front(:,2), '^'); 
scatter(nd_fronta(:,1), nd_fronta(:,2), '*'); 




toc
rmpath(problem_folder); 
rmpath(upperfolder);
rmpath(solver_folder);
rmpath(sort_folder);