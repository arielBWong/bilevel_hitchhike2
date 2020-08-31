%%this script is to run ulego with multiple seeds
%
clearvars;
close all;
problem_folder = strcat(pwd,'\problems\MOBP');
dace_folder = strcat(pwd,'\dace');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\TP');
addpath(problem_folder);

problem_folder = strcat(pwd,'\problems\DS');
addpath(problem_folder);
gsolver = strcat(pwd,'\globalsolver');
addpath(gsolver);

% problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
% problems = { 'mobp9(6)','mobp9(7)','mobp9(8)','mobp9(9)','mobp9(10)','mobp9(11)','mobp9(12)','mobp9(13)','mobp9(14)'};
addpath(dace_folder);
problems = { 'tp1()', 'tp2(6)' ,'tp3()' ,'tp4()' , 'ds1(6)', 'ds2(6)', 'ds3(6)', 'ds4(3,2)', 'ds5(3, 2)',  'mobp5()', 'mobp7()', 'mobp8()','mobp9(6)','mobp10()', 'mobp11(6)' };

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
    ulego_sao_archiveinsert(paras{i}{1}, paras{i}{2}, 'normalization_nd');
   
end

rmpath(problem_folder);