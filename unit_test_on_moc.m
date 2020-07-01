%%
% this unit test is for unit test on EG problems, which are MOC(multiple objective and constraint)
% similar to unit test 5 for MO
% 
clearvars;
close all;
workdir = pwd;
eim_methods = {'EIMnext'; 'EIMnext_znorm' };
test_problems = {'BNH()','SRN','Welded_Beam'}; 
np = length(test_problems);
ne = length(eim_methods);

problem_folder1 = strcat(pwd,'\problems\EGproblems');
addpath(problem_folder1);
% 
% paras={};
% for i=1:np
%     for j = 1:ne
%         paras{end+1}={test_problems{i}, eim_methods{j}};
%     end
% end


% moc_opt('EIMnext_znorm','BNH()', 50);

% npara = length(paras);
% parfor i = 1:npara
%     moc_opt(paras{i}{2}, paras{i}{1});
% end
rmpath(problem_folder1);


%------------where to start running paper demo
workdir = pwd;
problem_folder1 = strcat(pwd,'\Multiobjective_EGO_algorithms-master');
addpath(problem_folder1);
% Main_Constrained_Multiobjective_EGO('BNH');
% for i = 1:np
%     fun_name = test_problems{i}(1:end-2);
%     Main_Constrained_Multiobjective_EGO(fun_name);
%     
% end
rmpath(problem_folder1);
