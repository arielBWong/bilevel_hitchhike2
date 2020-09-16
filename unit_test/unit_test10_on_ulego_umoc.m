%% this test is single run of the main process, 
% this test is to confirm the process and 
% result of bilevel optimization on upper level moc lower level soc
% round 1: to be passed
% -----------------------
seed = 1;



workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
problem_folder = strcat(upperfolder,'\problems\MOBP');
solver_folder = strcat(upperfolder,'\globalsolver');
sort_folder = strcat(upperfolder,'\ND_Sort');

problem_folder = strcat(upperfolder,'\problems\DSM');
addpath(problem_folder);

addpath(problem_folder);
addpath(upperfolder);
addpath(solver_folder);
addpath(sort_folder);


tic;
% prob = 'dsm1dc1(3, 3)';
% ulego_umoc(prob, seed, 'EIMnext_znorm', 'EIM_eval', 'normalization_nd', 'EIMnext_znorm');
% ulego_umoc(prob, seed, 'EIMnext_znorm', 'Ehv_eval', 'normalization_nd', 'Believer_next');


% prob = eval(prob);
% savepath = strcat(pwd, '\result_folder\', prob.name, '_EIM_eval' );
% savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv');
% nd_front = csvread(savename_fu);
% scatter(nd_front(:,1), nd_front(:,2)); drawnow;

problems = {'dsm1(3, 3)', 'dsm1d(3, 3)','dsm1dc1(3, 3)','dsm1dc2(3,3)',...
             'dsm2(3, 3)', 'dsm2d(3, 3)','dsm2dc1(3, 3)','dsm2dc2(3, 3)',...
             'dsm3(3, 3)', 'dsm3d(3, 3)','dsm3dc1(3, 3)','dsm3dc2(3, 3)'};
seeds=[1];
algs = {'EIM_eval'}; 
np = length(problems);
ns = length(seeds);
na = length(algs);

paras=cell(1, np*ns * na);
nn = 1;
for i = 1:np
    for j = 1:ns
        for k = 1:na
            paras{ nn } = {problems{i}, seeds(j), algs{k}};
            nn = nn + 1;
        end
    end
end

nrun = length(paras);
parfor i = 1:nrun
    ulego_umoc(paras{i}{1}, paras{i}{2},'EIMnext_znorm' , paras{i}{3}, 'normalization_nd',  'EIMnext_znorm');
    % ulego_umoc(prob, seed, 'EIMnext_znorm', 'Ehv_eval', 'normalization_nd', 'Believer_next');
    ulego_sao_archiveinsert(paras{i}{1}, paras{i}{2}, 'normalization_nd');
    ulego_sao_pop(paras{i}{1}, paras{i}{2}, 'normalization_nd');
end

toc
rmpath(problem_folder); 
rmpath(upperfolder);
rmpath(solver_folder);
rmpath(sort_folder);