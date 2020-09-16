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

problem_folder = strcat(pwd,'\problems\DSM');
addpath(problem_folder);



solver_folder = strcat(pwd,'\globalsolver');
addpath(solver_folder);

% problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
% problems = { 'mobp9(6)','mobp9(7)','mobp9(8)','mobp9(9)','mobp9(10)','mobp9(11)','mobp9(12)','mobp9(13)','mobp9(14)'};

% problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()', 'mobp11(6)' };
% problems = {'tp1()',  'tp2(3)' ,'tp3()' ,'tp4()' , 'ds1(6)', 'ds2(6)', 'ds3(6)', 'ds4(3,2)', 'ds5(3, 2)',  'mobp5()', 'mobp7()', 'mobp8()','mobp9(6)','mobp10()', 'mobp11(6)',  'dsm1(3)'  };


algs = {'EIM_eval'};  %  'Ehv_eval', 

problems = { 'dsm1(2, 2)', 'dsm1d(2, 2)','dsm1dc1(2, 2)','dsm1dc2(2,2)',...
             'dsm2(2, 2)', 'dsm2d(2, 2)','dsm2dc1(2, 2)','dsm2dc2(2, 2)',...
             'dsm3(2, 2)', 'dsm3d(2, 2)','dsm3dc1(2, 2)','dsm3dc2(2, 2)',...
             'dsm1(3, 3)', 'dsm1d(3, 3)','dsm1dc1(3, 3)','dsm1dc2(3, 3)',...
             'dsm2(3, 3)', 'dsm2d(3, 3)','dsm2dc1(3, 3)','dsm2dc2(3, 3)',...
             'dsm3(3, 3)', 'dsm3d(3, 3)','dsm3dc1(3, 3)','dsm3dc2(3, 3)'...
             'dsm1(4, 4)', 'dsm1d(4, 4)','dsm1dc1(4, 4)','dsm1dc2(4, 4)',...
             'dsm2(4, 4)', 'dsm2d(4, 4)','dsm2dc1(4, 4)','dsm2dc2(4, 4)',...
             'dsm3(4, 4)', 'dsm3d(4, 4)','dsm3dc1(4, 4)','dsm3dc2(4, 4)'};
seeds = linspace(1, 11, 11);

np = length(problems);
ns = length(seeds);
na = length(algs);

paras=cell(1, np*ns * na);
nn = 1;
for i = 1:np
    for j = 1:ns
        for k = 1:na
            paras{ nn } = {problems{i}, seeds(j), algs{k}};
            nn = nn + 1;
        end
    end
end

nrun = length(paras);
parfor i = 1:nrun
    ulego_umoc(paras{i}{1}, paras{i}{2},'EIMnext_znorm' , paras{i}{3}, 'normalization_nd',  'EIMnext_znorm');
    % ulego_umoc(prob, seed, 'EIMnext_znorm', 'Ehv_eval', 'normalization_nd', 'Believer_next');
    ulego_sao_archiveinsert(paras{i}{1}, paras{i}{2}, 'normalization_nd');
    ulego_sao_pop(paras{i}{1}, paras{i}{2}, 'normalization_nd');
end

rmpath(problem_folder);