%% this script extract empirical_pf from exghausive search results
% read result folder from seeds
% combine them together and extract nd front
% save them into a pf file named 'empf_{problemname}'
clearvars;
close all;

problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
problem_folder = strcat(pwd,'\problems\MOBP');
addpath(problem_folder);

np = length(problems);

for ii = 1:np
    prob = eval(problems{ii});
    fu = [];
    for jj = 1:seedmax
        % extract nd front from each seed
        savename =  strcat(pwd, '\result_folder\', prob.name, '_emp_pf');
        empf_filename = strcat(savename, '\fu_', str(jj), '.csv');
        nd_front = csvread(empf_filename);
        fu = [nd_front; fu];
    end
    
    % reconduct pareto front on this collected results
    nd_index = Paretoset(fu);
    nd_front = fu(nd_index, :);
    
    % save empf to nd_front
    outname =  strcat(pwd, '\result_folder\', prob.name, '_emp_pf.csv');
    csvread(outname, nd_front);
    
end