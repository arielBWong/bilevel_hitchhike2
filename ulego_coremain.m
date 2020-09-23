%%this script is to run ulego with multiple seeds
%
clearvars;
close all;
problem_folder = strcat(pwd,'\problems\SMD');
addpath(problem_folder);
problem_folder = strcat(pwd,'\globalsolver');
addpath(problem_folder);
problem_folder = strcat(pwd,'\ND_Sort');
addpath(problem_folder);

% problems = { 'smd1()','smd2()','smd3()','smd4()','smd5()','smd6()','smd7()',...
             % 'smd8()','smd9()', 'smd10()','smd11()','smd12()'};
% problems = { 'bltp1()','bltp2()','bltp3()','bltp4()','bltp5()','bltp6()','bltp7()',...
             % 'bltp8()','bltp9()', 'bltp10()','bltp11()'};
         
% problems = { 'smd1()','smd2()','smd3()','smd4()','smd5()','smd6()','smd7()', 'smd8()'};
problems = { 'smd12()'};
seeds = linspace(9, 9, 9);

np = length(problems);
ns = length(seeds);

paras=cell(1, np*ns);
for i = 1:np
    for j = 1:ns
        paras{ (i-1)*ns +j } = {problems{i}, seeds(j)};
    end
end

nrun = length(paras);

for i = 1:nrun
     % ulego_coreeim(paras{i}{1}, paras{i}{2}, 'EIMnext', 'EIM_eval',  'normalization_y');
     ulego_coreble(paras{i}{1}, paras{i}{2}, 'normalization_y');
     % ulego_coregen(paras{i}{1}, paras{i}{2}, 'normalization_y');
end

% parfor i = 1:nrun
%     ulego(paras{i}{1}, paras{i}{2},'EIMnext_znorm' , 'EIM_eval', 'normalization_y');
% end

rmpath(problem_folder);