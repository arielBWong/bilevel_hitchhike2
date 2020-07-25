%% this script is used to extract empirical pf from generated results
% similar to resconvert script
% (1) read across seed
% (2) stack nd and re-extract nd
% (3) save to certain file

clearvars;
close all;


seedmax = 4;
% read a seed and
problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
problem_folder = strcat(pwd,'\problems\MOBP');
addpath(problem_folder);

np = length(problems);

for ii = 1: np
    prob = eval(problems{ii});
    allnds = [];
    % (1) read across seed
    for jj = 1: seedmax
        savepath = strcat(pwd, '\result_folder\', prob.name, '_emp_pf');
        savename_fu = strcat(savepath, '\fu_', num2str(jj),'.csv');
        
        nd_front = csvread(savename_fu);
        % (2) stack nd
        allnds = [allnds; nd_front];
    end
    % (2) re-extract nd
    nd_index = Paretoset(allnds);
    nd_front = allnds(nd_index, :);
    num_front = size(nd_front, 1);
    disp(num_front);
    
    % (3) save to certain file
    savename_pf = strcat(pwd, '\result_folder\', prob.name, '_empf_', num2str(num_front), '.csv');
    csvwrite(savename_pf, nd_front);
    
    % (4) plot to have a look
    fig1 = gcf;
    scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
    t = [prob.name,' empirical pf'];
    title(t);
    savename = strcat(pwd, '\result_folder\', prob.name, '_empf.fig');
    savefig(savename);
    savename = strcat(pwd, '\result_folder\', prob.name, '_empf.png');
    
    saveas(fig1, savename);
    close(fig1);

end