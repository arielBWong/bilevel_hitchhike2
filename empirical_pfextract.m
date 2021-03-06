%% this script is used to extract empirical pf from generated results
% similar to resconvert script
% (1) read across seed
% (2) stack nd and re-extract nd
% (3) save to certain file

clearvars;
close all;


seedmax =3;
% read a seed and
% problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
problems = { 'tp1()' ,'tp2(6)' ,'tp3()' ,'tp4()' , 'ds1(6)', 'ds2(6)', 'ds3(6)', 'ds4(3,2)', 'ds5(3, 2)'};
% problems = { 'mobp7()'};
problems = { 'dsm1(3)', 'dsm2(3)', 'dsm3(3)'};

problem_folder = strcat(pwd,'\problems\MOBP');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\DSM');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\TP');
addpath(problem_folder);
problem_folder = strcat(pwd,'\evenpf');
addpath(problem_folder);


np = length(problems);

for ii = 1:np
    prob = eval(problems{ii});
    
    if ~strcmp(prob.name(1:end-1), 'dsm')
        fu = [];
        if  strcmp(prob.name, 'mobp11')
            seedmax =3;
        end
        for jj = 1:seedmax
            % extract nd front from each seed
            savename =  strcat(pwd, '\result_folder\', prob.name, '_emp_pf');
            empf_filename = strcat(savename, '\fu_', num2str(jj), '.csv')
            nd_front = csvread(empf_filename);
            fu = [nd_front; fu];
        end
        
        % reconduct pareto front on this collected results
        nd_index = Paretoset(fu);
        nd_front = fu(nd_index, :);
        
        % select evenly 100 points
        if size(nd_front,1) > 200
            guarantee = 100;
            [id_fronts,f_fronts,NumComp,NumFronts] = E_NDSort_c(nd_front);
            updated_order=Sparse_selection(id_fronts,f_fronts,guarantee);
            nd_front = nd_front(updated_order(1:100), :);
        end
    else
        nd_front = prob.upper_pf(100);
    end
    
    % save empf to nd_front
    outname =  strcat(pwd, '\result_folder\', prob.name, '_emp_pf.csv');
    csvwrite(outname, nd_front);
end

% plot nd and count number of nd save it in the file name
for ii = 1: np
    prob = eval(problems{ii});
    if  ~strcmp(prob.name(1:end-1), 'dsm')
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
    else
        
        % save empf to nd_front
        outname =  strcat(pwd, '\result_folder\', prob.name, '_emp_pf.csv');
        nd_front = csvread(outname);
        
        fig1 = gcf;
        scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
        t = [prob.name,' mathmetical pf'];
        title(t);
        savename = strcat(pwd, '\result_folder\', prob.name, '_empf.fig');
        savefig(savename);
        savename = strcat(pwd, '\result_folder\', prob.name, '_empf.png');
        
        saveas(fig1, savename);
        close(fig1);
    end
    
end