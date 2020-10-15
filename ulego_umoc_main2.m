%%this script is to run ulego with multiple seeds
% need to update
clearvars;
close all;
problem_folder = strcat(pwd,'\problems\MOBP');
addpath(problem_folder);

dace_folder = strcat(pwd,'\dace');
addpath(dace_folder);

problem_folder = strcat(pwd,'\problems\DS');
addpath(problem_folder);

gsolver = strcat(pwd,'\globalsolver');
addpath(gsolver);

gsolver = strcat(pwd,'\ND_Sort');
addpath(gsolver);

problem_folder = strcat(pwd,'\problems\DSM');
addpath(problem_folder);



solver_folder = strcat(pwd,'\globalsolver');
addpath(solver_folder);

         
problems = {  'dsm1(2, 2)','dsm1d(2, 2)',...
            'dsm2(2, 2)','dsm2d(2, 2)',...
            'dsm3(2, 2)','dsm3d(2, 2)',...
            'dsm1dc1(2, 2)', 'dsm2dc1(2, 2)', 'dsm3dc1(2, 2)'};
        

seeds = linspace(1, 11, 11);

np = length(problems);
ns = length(seeds);


paras=cell(1, np*ns);
nn = 1;
for i = 1:np
    for j = 1:ns
            paras{ nn } = {problems{i}, seeds(j)};
            nn = nn + 1;
    end
end

nrun = length(paras);
parfor i = 1:nrun
    tic;
    ulego_umoc_hyb(paras{i}{1}, paras{i}{2},'EIMnext' , 'EIM_eval', 'normalization_nd',  'EIMnext');
    toc;
     % ulego_umoc(paras{i}{1}, paras{i}{2},'EIMnext' , 'EIM_eval', 'normalization_nd',  'EIMnext');
     % ulego_sao_archiveinsert(paras{i}{1}, paras{i}{2}, 'normalization_nd');
    % ulego_sao_pop(paras{i}{1}, paras{i}{2}, 'normalization_nd');
end

rmpath(problem_folder);