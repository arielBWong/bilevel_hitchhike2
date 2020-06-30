%% partest


workdir = pwd;
problem_folder = strcat(pwd,'\problems\EGproblems');
addpath(problem_folder);


problems = { 'SRN()','BNH()', 'Welded_Beam()'}; % 'TNK()'};
eim_methods = {'EIMnext_znorm'}; % ,'EIMnext_znorm'};
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


rmpath(problem_folder)