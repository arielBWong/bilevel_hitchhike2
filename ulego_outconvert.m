% this script generate median igd and its corresponding
%

clearvars;
close all;

seedmax = 11;
problems = { 'dsm1(2, 2)','dsm1d(2, 2)','dsm1dc1(2, 2)', ...
             'dsm2(2, 2)', 'dsm2d(2, 2)','dsm2dc1(2, 2)', ...
             'dsm3(2, 2)', 'dsm3d(2, 2)', 'dsm3dc1(2, 2)'};

methods = { 'EIM_eval', 'sao_archiveinsert', 'hyb'};
np = length(problems);
nm = length(methods);
ndmatrix_problems = zeros(seedmax, nm * np);
for ii = 1: np
    
    prob = eval(problems{ii});
    % compare results by problem
    % (1) for each problem, load empirical pf
    if strcmp(prob.name(1:3), 'dsm')
        empf = prob.upper_pf(100);
    else
        empfsave = strcat(pwd, '\result_folder\', prob.name, '_emp_pf.csv');
        empf = csvread(empfsave);
    end
    
    
    igd_fu = cell(seedmax, nm);
    
    for seed = 1: seedmax
        for jj = 1:nm
            
            method = methods{jj};
            savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(prob.n_lvar), '_', method);
            % savepath = strcat(pwd, '\result_folder\', prob.name,  '_', method);
            savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv')
            nd_front = csvread(savename_fu);
            size(empf);
            size(nd_front);
            % (3) calculate igd for this seed
            igd_fu{seed, jj} = mean(min(pdist2(empf,nd_front),[],2));
            ndmatrix_problems(seed, (ii -1) *nm +jj)= igd_fu{seed, jj};
        end
    end
    %----for ii -- problem routine
end
%---process ndmatrix_problems----
statistic_matrix = zeros(6, nm * np); %(mean, std, median, median_id, nn)
for ii = 1 : np
    % plot across problems
    % (1) read in empf
    prob = eval(problems{ii});
    if strcmp(prob.name(1:3), 'dsm')
        empf = prob.upper_pf(100);
    else
        empfsave = strcat(pwd, '\result_folder\', prob.name, '_emp_pf.csv');
        empf = csvread(empfsave);
    end
    % (2) read median of three different methods
    medianlist = [];
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
        medianlist = [medianlist, ids(middle)];
        
        % use this median id to create lower level match accuracy
        % median function evaluation
        prob = eval(problems{ii});
        method = methods{jj};
        lower_median = lower_performance(ids(middle), prob, method);
        statistic_matrix(6, (ii-1) * nm + jj)  = lower_median;
        
        savepath = strcat(pwd, '\result_folder\', prob.name, '_',  num2str(prob.n_lvar), '_',method);
        % savepath = strcat(pwd, '\result_folder\', prob.name, '_', method);
        savename_nn = strcat(savepath, '\nn_', num2str(ids(middle)),'.csv')
        nn = csvread(savename_nn);
        statistic_matrix(5, (ii-1) * nm + jj)  =  sum(nn);
    end
end


%-to csv
prob = eval(problems{1});

savepath = strcat(pwd, '\result_folder\',prob.name, '_K', num2str(prob.n_lvar), '_median.csv')
fp=fopen(savepath,'w');
fprintf(fp, 'problem,');
for jj = 1:nm
    fprintf(fp, '%s,', methods{jj}  );
    fprintf(fp, 'lower_performance,');
end
fprintf(fp, '\n');


for ii = 1:np
    prob = eval(problems{ii});
    name = strcat(prob.name, '_', num2str(prob.n_lvar));
     fprintf(fp, '%s,', name);
     for jj = 1:length(methods)
        fprintf(fp, '%f,',  statistic_matrix(3, (ii-1) * length(methods)+jj));
        fprintf(fp, '%f,',  statistic_matrix(6, (ii-1) * length(methods)+jj));
     end
     fprintf(fp, '\n');
end
fclose(fp);

function median_accuracy = lower_performance(seed, prob, method)

savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(prob.n_lvar), '_', method);
% savepath = strcat(pwd, '\result_folder\', prob.name,  '_', method);
savename_xu = strcat(savepath, '\xu_raw_', num2str(seed),'.csv');
savename_xl = strcat(savepath, '\xl_raw_', num2str(seed),'.csv');

xu = csvread(savename_xu);
xl = csvread(savename_xl);

fl = prob.evaluate_l(xu, xl);
median_accuracy = median(fl);


end