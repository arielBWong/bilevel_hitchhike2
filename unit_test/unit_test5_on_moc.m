%% this unit test is to call

clearvars;
close all;


workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
addpath(upperfolder);



eim_methods = {'EIMnext'; 'EIMnext_znorm' };
test_problems = {'ZDT1()','ZDT2()','ZDT3()'}; %,'DTLZ2()','DTLZ5()','DTLZ7()'};
np = length(test_problems);
ne = length(eim_methods);
%-------------------my code unit test--------------
% problem_folder1 = strcat(upperfolder,'\problems\EGproblems');
% problem_folder2 = strcat(upperfolder,'\problems\ZDT');
% problem_folder3 = strcat(upperfolder,'\problems\DLTZ');
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
problem_folder1 = strcat(upperfolder,'\Multiobjective_EGO_algorithms-master');
addpath(problem_folder1);
for i = 1:np
    fun_name = test_problems{i}(1:end-2);
    Main_Standard_Multiobjective_EGO(fun_name);
    
end
rmpath(problem_folder1);
rmpath(upperfolder);
