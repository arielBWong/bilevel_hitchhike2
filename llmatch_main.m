%%
clearvars;
close all;
problem_folder = strcat(pwd,'\ND_Sort');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\TP2');
addpath(problem_folder);
 
problems = {'smd1()', 'smd2()','smd3()', 'smd4()',  'smd5()',   'smd6()', 'smd7()', 'smd8()',  'smd9()',   'smd10()', 'smd11()', 'smd12()',...
    'dsm1(2,2)','dsm1(3,3)', 'dsm1(4,4)','dsm1(5,5)','dsm1dc1(2,2)','dsm1dc1(3,3)', 'dsm1dc1(4,4)',  'dsm1dc1(5,5)'};


problems = {};
k = 2; 
init = 11 * k - 1;
problems ={
      'SHCBc()', 'newBranin5()','newBranin2()', 'Reverse_Mystery()' , 'Mystery()', 'Haupt_schewefel()', 'Gpc()', 'Gomez3()'...
    };

seeds = linspace(1, 29, 29);


match_methods = {'llmatch_sao_archiveinsert', 'llmatch', 'llmatch_believeradapt'}; % , 'llmatch_sao_archiveinsert', 'llmatch',  'llmatch_hyb', 'llmatch_believeradapt'


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

for i =1:nrun
% llmatch_globalmin_cmp(paras{i}{1}, paras{i}{3},  paras{i}{2});
llmatch_behaviourstudy(paras{i}{1}, paras{i}{3},  paras{i}{2}, init);
end
