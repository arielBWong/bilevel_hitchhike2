%% this method process experiment results
%
clearvars;
close all;


seedmax = 1;
% read a seed and
problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
methods = {'Ehv_eval', 'EIM_eval'};
problem_folder = strcat(pwd,'\problems\MOBP');
addpath(problem_folder);

np = length(problems);
nm = length(methods);

ndmatrix_problems = zeros(seedmax, nm * np);
for ii = 1: np
    prob = eval(problems{ii});
    %compare results by problem
    nd_fu = cell(seedmax, nm);
    
    
    for seed = 1: seedmax
        for jj = 1:nm
            
            method = methods{jj};
            savepath = strcat(pwd, '\result_folder\', prob.name, '_', method);
            savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv');
            
            % --- save plots----
            nd_front = csvread(savename_fu);
            
            fig1 = gcf;
            scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
            t = [prob.name,' ', method, ' seed ', num2str(seed) ];
            title(t);
            savename = strcat(savepath, '\nd_', num2str(seed),'.fig');
            savefig(savename);
            savename = strcat(savepath, '\nd_', num2str(seed),'.png');
            
            saveas(fig1, savename);
            close(fig1);
            %----
            
            %--- calculate hyper volume results---
            nd_fu{seed, jj} = nd_front;
        end
    end
    % after record the raw nd front
    % calculate hyper volume
    % (1) find nadir and ideal (2) normalize all nd with nadir and ideal
    % (3) calculate hyper volume
    allnd = [];
    for seed = 1:seedmax
        for jj = 1:nm
            allnd= [allnd; nd_fu{seed, jj} ];
        end
    end
    num_obj = size(allnd, 2);
    % (1) find nadir and ideal
    nadir = max(allnd) ;
    ideal = min(allnd);
    ref = ones(1, num_obj) * 1.1;
    for seed = 1:seedmax
        for jj = 1:nm
            %  (2) normalize all nd with nadir and ideal
            nd_fu{seed, jj} = (nd_fu{seed, jj} - ideal) ./ (nadir - ideal);
            % (3) calculate hyper volume
            ndmatrix_problems(seed, (ii -1) *nm +jj)= Hypervolume(nd_fu{seed, jj} , ref);
        end
    end
    
    %----for ii -- problem routine
end

%---process ndmatrix_problems----
statistic_matrix = zeros(4, nm * np); %(mean, std, median, median_id)
for ii = 1:np
    for jj = 1: nm
        one_column = ndmatrix_problems(:, (ii-1) * nm + jj) ;
        % mean
        statistic_matrix(1, (ii-1) * nm + jj) =  mean(one_column);
        % std
        statistic_matrix(2, (ii-1) * nm + jj) =  std(one_column);
        % median
        statistic_matrix(3, (ii-1) * nm + jj)  = median(one_column);
        % median id
        middle = round(seedmax/2);
        [~, ids] = sort(one_column);
        statistic_matrix(4, (ii-1) * nm + jj)  =  ids(middle);
    end
end


% save into csv
savepath = strcat(pwd, '\result_folder\mobp_res.csv');
fp=fopen(savepath,'w');
fprintf(fp, 'seed,');
% format header
for ii = 1:np
    for jj = 1:nm
        single_header = strcat(problems{ii}, '_', methods{jj});
        fprintf(fp, '%s,', single_header);
    end
end
fprintf(fp, '\n');
% format raw record
for seed = 1:seedmax
    fprintf(fp, '%d,', seed);
    for ii = 1:np
        for jj = 1:nm
            ele = ndmatrix_problems(seed, (ii -1) *nm +jj);
            fprintf(fp, '%f,', ele);
        end
    end
    fprintf(fp, '\n');
end
% format statistics
st = {'mean', 'std', 'median', 'median id'};
for s = 1:4
    fprintf(fp, st{s});
    for ii = 1:nm * np
        fprintf(fp, '%f,', statistic_matrix(s, ii) );
    end
    fprintf('\n');
end
fclose(fp);
rmpath(problem_folder);
