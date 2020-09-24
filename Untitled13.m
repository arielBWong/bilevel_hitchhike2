%
problem_folder = strcat(pwd,'\problems\SMD');
addpath(problem_folder);
problem_folder = strcat(pwd,'\globalsolver');
addpath(problem_folder);
problem_folder = strcat(pwd,'\ND_Sort');
addpath(problem_folder);


prob = smd10();
algs = {'eim', 'bel', 'gen'};

seed = 11;
rng(seed, 'twister');
% savepath = strcat(pwd, '\result_folder\', prob.name );
% filename = strcat(savepath, '_',  algs{2});
% xusave = strcat(filename, '\xu_', num2str(seed), '.csv')
% xu  = csvread(xusave);   
% 
% xlsave = strcat(filename, '\xl_', num2str(seed), '.csv')
% xl  = csvread(xlsave);    

xu = [1, 1];
[xl,  n_fev, flag] = llmatch(xu, prob, 20, 20, 'EIMnext', 40, 'EIM_eval', 9);

[fu, fc] = prob.evaluate_u(xu, xl)

[xl,  n_fev, flag] =  llmatch_sao_archiveinsert(xu, prob, 20, 800, 20);
[fu, fc] =prob.evaluate_u(xu, xl)
