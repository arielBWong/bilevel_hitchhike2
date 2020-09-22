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

         
problems = { 'dsm1(3, 3)', 'dsm1d(3, 3)','dsm1dc1(3, 3)',...
                'dsm2(3, 3)', 'dsm2d(3, 3)','dsm21dc1(3, 3)', ...
                'dsm3(3, 3)', 'dsm3d(3, 3)','dsm3dc12(3, 3)'};


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
for i = 1:nrun
    ulego_umoc(paras{i}{1}, paras{i}{2},'EIMnext' , 'EIM_eval', 'normalization_nd',  'EIMnext');
    ulego_sao_archiveinsert(paras{i}{1}, paras{i}{2}, 'normalization_nd');
    ulego_sao_pop(paras{i}{1}, paras{i}{2}, 'normalization_nd');
end

rmpath(problem_folder);