%% partest


workdir = pwd;
problem_folder1 = strcat(pwd,'\problems\EGproblems');
problem_folder2 = strcat(pwd,'\problems\ZDT');
problem_folder3 = strcat(pwd,'\problems\DTLZ');

addpath(problem_folder1);
addpath(problem_folder2);
addpath(problem_folder);


%problems = { 'SRN()','BNH()', 'Welded_Beam()'}; % 'TNK()'};
eim_methods = {'EIMnext_znorm' ,'EIMnext'};
test_problems = {'ZDT1()','ZDT2()','ZDT3()','DTLZ2()','DTLZ5()','DTLZ7()'};
np = length(problems);
ne = length(eim_methods);


paras={};
for i=1:np
    for j = 1:ne
        paras{end+1}={problems{i}, eim_methods{j}};
    end
end

npara = length(paras);
parfor i = 1:npara
    moc_opt(paras{i}{2}, paras{i}{1});
end


rmpath(problem_folder1);
rmpath(problem_folder2);
rmpath(problem_folder3);