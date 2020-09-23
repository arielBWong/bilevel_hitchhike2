%%
clearvars;
close all;
% problems = { 'dsm1(3, 3)', 'dsm1d(3, 3)','dsm1dc1(3, 3)','dsm1dc2(3, 3)',...
%          'dsm2(3, 3)', 'dsm2d(3, 3)','dsm2dc1(3, 3)','dsm2dc2(3, 3)',...
%               'dsm3(3, 3)', 'dsm3d(3, 3)','dsm3dc1(3, 3)','dsm3dc2(3, 3)' };
%           
 
%           problems = { 'dsm1(2,2)', 'dsm1d(2,2)','dsm1dc1(2,2)','dsm1dc2(2,2)',...
%          'dsm2(2,2)', 'dsm2d(2,2)','dsm2dc1(2,2)','dsm2dc2(2,2)',...
%               'dsm3(2,2)', 'dsm3d(2,2)','dsm3dc1(2,2)','dsm3dc2(2,2)' }; % change p3 term back to scale 10
%           
          
 % problems = { 'dsm1(4,4)', 'dsm1d(4,4)','dsm1dc1(4,4)','dsm1dc2(4,4)' }; % change p3 term back to scale 10
  % problems = { 'dsm1(5, 5)', 'dsm1d(5, 5)','dsm1dc1(5, 5)','dsm1dc2(5, 5)' };     
  % problems = { 'dsm1(3,3)', 'dsm1d(3,3)','dsm1dc1(3,3)','dsm1dc2(3,3)' };   
 % problems = { 'dsm1(2,2)', 'dsm1d(2,2)','dsm1dc1(2,2)','dsm1dc2(2,2)' }; 
 
 problems = {'smd12()'};
seeds = linspace(1, 11, 11);

match_methods = {'llmatch',  'llmatch_sao_archiveinsert'}; % 'llmatch_sao_pop', ,  'llmatch_sao_archiveinsert'

% llmatch_globalmin_cmp(problems{1}, match_methods{1}, 1);

np = length(problems);
ns = length(seeds);
na = length(match_methods);

paras=cell(1, np*ns * na);
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
llmatch_globalmin_cmp(paras{i}{1}, paras{i}{3},  paras{i}{2});
end