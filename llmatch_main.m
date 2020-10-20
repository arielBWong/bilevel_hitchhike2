%%
clearvars;
close all;
problem_folder = strcat(pwd,'\ND_Sort');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\TP2');
addpath(problem_folder);
 
problems = {'smd1()', 'smd2()','smd3()', 'smd4()',  'smd5()',   'smd6()', 'smd7()', 'smd8()',  'smd9()',   'smd10()', 'smd11()', 'smd12()',...
    'dsm1(2,2)','dsm1(3,3)', 'dsm1(4,4)','dsm1(5,5)','dsm1dc1(2,2)','dsm1dc1(3,3)', 'dsm1dc1(4,4)',  'dsm1dc1(5,5)'};

k = 2; 
init = 11 * k - 1;
problems ={
     'tp3(2, 2)', ...
    };


% 'tp3(5,5)','tp5(5,5)', 'tp6(5,5)', 'tp7(5,5)','tp8(5,5)','tp9(5,5)'
% 'tp3(4,4)','tp5(4,4)', 'tp6(4,4)', 'tp7(4,4)','tp8(4,4)','tp9(4,4)'
% 'tp3(3,3)','tp5(3,3)', 'tp6(3,3)','tp7(3,3)','tp8(3,3)','tp9(3,3)'
% 'tp3(2,2)','tp5(2,2)', 'tp6(2,2)','tp7(2,2)','tp8(2,2)','tp9(2,2)'
seeds = linspace(1, 29, 29);

match_methods = {'llmatch_believeradapt', 'llmatch', 'llmatch_sao_archiveinsert'}; % , 'llmatch_sao_archiveinsert', 'llmatch',  'llmatch_hyb', 'llmatch_believeradapt'


np = length(problems);
ns = length(seeds);
na = length(match_methods);

paras=cell(1, np * ns * na);
nn = 1;
for i = 1:np
    for j = 1:ns
        for k = 1:na
            paras{ nn } = {problems{i}, seeds(j), match_methods{k}};
            nn = nn + 1;
        end
    end
end

nrun = length(paras);

parfor i =1:nrun
% llmatch_globalmin_cmp(paras{i}{1}, paras{i}{3},  paras{i}{2});
llmatch_behaviourstudy(paras{i}{1}, paras{i}{3},  paras{i}{2}, init);
end
