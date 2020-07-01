%% partest

clearvars;
close all;
workdir = pwd;
eim_methods = {'EIMnext'; 'EIMnext_znorm' };
test_problems = {'ZDT1()','ZDT2()','ZDT3()'}; %,'DTLZ2()','DTLZ5()','DTLZ7()'};
np = length(test_problems);
ne = length(eim_methods);
%-------------------my code unit test--------------
% problem_folder1 = strcat(pwd,'\problems\EGproblems');
% problem_folder2 = strcat(pwd,'\problems\ZDT');
% problem_folder3 = strcat(pwd,'\problems\DLTZ');
% 
% addpath(problem_folder1);
% addpath(problem_folder2);
% addpath(problem_folder3);
% 
% paras={};
% for i=1:np
%     for j = 1:ne
%         paras{end+1}={test_problems{i}, eim_methods{j}};
%     end
% end
% %i = 1;
% %moc_opt('EIMnext_znorm','ZDT2()');
% 
% npara = length(paras);
% parfor i = 1:npara
%     moc_opt(paras{i}{2}, paras{i}{1});
% end
% rmpath(problem_folder1);
% rmpath(problem_folder2);
% rmpath(problem_folder3);
%-------------------my code unit test--------------

%------------where to start running paper demo
workdir = pwd;
problem_folder1 = strcat(pwd,'\Multiobjective_EGO_algorithms-master');
addpath(problem_folder1);
for i = 1:np
    fun_name = test_problems{i}(1:end-2);
    Main_Standard_Multiobjective_EGO(fun_name);
    
end
rmpath(problem_folder1);
