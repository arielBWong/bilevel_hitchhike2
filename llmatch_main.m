%%
clearvars;
close all;
problem_folder = strcat(pwd,'\ND_Sort');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\EGproblems');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\TP2');
addpath(problem_folder);
 
% unimodal 'cec2010(18,2)', 'cec2010(17,2)', cec2010(14,2), cec2010(9,2),
% 'cec2010(7,2)', 'cec2010(5,2)', 'cec2010(4,2)', 'cec2010(3,2)', 'cec2010(2,2)';

% multimodal 'cec2010(16,2)', 'cec2010(13,2)', 'cec2010(12,2)',
% 'cec2010(1,2)';
k = 4; 
init = 11 * k - 1;
problems ={ 'cec2010(1,2)', 'tp3(2, 2)', 'tp5(2, 2)', 'tp7(2, 2)','Shekel(2, 2)','rastrigin(2, 2)', ... 
     'SHCBc()', 'newBranin5()','newBranin2()', 'Reverse_Mystery()' , 'Mystery()', 'Haupt_schewefel()', 'Gpc()', 'Gomez3()',...
      'cec2010(1,3)', 'tp3(3, 3)', 'tp5(3, 3)', 'tp7(3, 3)','Shekel(3, 3)','rastrigin(3, 3)', ... 
};

% 
seeds = linspace(1, 29, 29);

match_methods = {'llmatch_believeradapt' }; % , 'llmatch_sao_archiveinsert', 'llmatch',  'llmatch_believeradapt'


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
