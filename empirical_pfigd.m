%% this function is used to generate igd
% from empirical pf
% steps
% (1) for each problem, load empirical pf
% (2) for each seed, load nd front
% (3) calculate idg for this seed
% (4) similar to converter, change write to csv file
clearvars;
close all;


seedmax = 11;
% read a seed and
problems = { 'dsm1(3, 3)', 'dsm1d(3, 3)','dsm1dc1(3, 3)',...
             'dsm2(3, 3)', 'dsm2d(3, 3)','dsm2dc1(3, 3)', ...
             'dsm3(3, 3)', 'dsm3d(3, 3)','dsm3dc1(3, 3)'};

met = {'EIM', 'BEL', 'GEN'};
met = {'HYB', 'EIM', 'BEL'};
problems ={ 'dsm1(2, 2)',  'dsm1(3, 3)'};
% methods = {'Ehv_eval', 'EIM_eval', 'ea_ea'};
methods = { 'hyb', 'EIM_eval', 'sao_archiveinsert'}; %, 'ea_ea','sao_onerand', 'sao_popinsert'

problem_folder = strcat(pwd,'\problems\MOBP');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\DS');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\DSM');
addpath(problem_folder);
problem_folder = strcat(pwd,'\problems\TP');
addpath(problem_folder);
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
            size(empf)
            size(nd_front)
            % (3) calculate igd for this seed
            igd_fu{seed, jj} = mean(min(pdist2(empf,nd_front),[],2));
            ndmatrix_problems(seed, (ii -1) *nm +jj)= igd_fu{seed, jj};
        end
    end
    %----for ii -- problem routine
end

%---process ndmatrix_problems----
statistic_matrix = zeros(5, nm * np); %(mean, std, median, median_id, nn)
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
        
        % median function evaluation
        prob = eval(problems{ii});
        method = methods{jj};
        savepath = strcat(pwd, '\result_folder\', prob.name, '_',  num2str(prob.n_lvar), '_',method);
        % savepath = strcat(pwd, '\result_folder\', prob.name, '_', method);
        savename_nn = strcat(savepath, '\nn_', num2str(ids(middle)),'.csv')
        nn = csvread(savename_nn);
        statistic_matrix(5, (ii-1) * nm + jj)  =  sum(nn);
    end
    
    % plot three of them
    fig1 = gcf;
    plot(empf(:,1), empf(:,2),'Color',[0.4660 0.6740 0.1880], 'LineWidth', 2, 'LineStyle', ':'); hold on;
    
    pattern = cell(1, nm);
    pattern{1} = '^';
    pattern{2} = 'o';
    pattern{3} = 'd';
    pattern{4} = 'o';
    

    color{1} = 'b';
    color{2} = 'k';
    color{3} = 'r';
    color{4} = 'g';
    
    for kk = 1:nm
        seed = medianlist(kk);
        method = methods{kk};
        savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(prob.n_lvar), '_',method);
        % savepath = strcat(pwd, '\result_folder\', prob.name, '_', method);
        savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv')
        nd_front = csvread(savename_fu);
        scatter(nd_front(:,1), nd_front(:,2), 80, pattern{kk}, color{kk},'filled'); drawnow;
    end
    numk = prob.n_lvar;
    t = [prob.name,' k', num2str(numk),' IGD median compare'];
    title(t,  'FontSize', 16);
    legend('PF', met{1}, met{2},  met{3},  'FontSize', 14); % methods{4}
    
    xlabel('F1', 'FontSize', 14);
    ylabel('F2', 'FontSize', 14);
    
    a = get(gca,'XTickLabel');
    set(gca, 'XTickLabel',a, 'FontSize', 12);
    a = get(gca,'YTickLabel');
    set(gca, 'YTickLabel',a, 'FontSize', 12);
    
    
    savename = strcat(pwd, '\result_folder\', prob.name,'_', num2str(numk), '_igdcompare.fig');
    savefig(savename);
    savename = strcat(pwd, '\result_folder\', prob.name,'_', num2str(numk), '_igdcompare.png');
    saveas(fig1, savename);
    close(fig1);
    
end


% save into csv
savepath = strcat(pwd, '\result_folder\dsm_igdres.csv');
fp=fopen(savepath,'w');
fprintf(fp, 'seed,');
% format header
for ii = 1:np
    for jj = 1:nm
        prob = eval(problems{ii});
        single_header = strcat(prob.name, '_', methods{jj},'_', num2str(prob.n_lvar));
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
st = {'mean', 'std', 'median', 'median id', 'nfev_median'};
for s = 1:length(st)
    fprintf(fp, '%s,', st{s});
    for ii = 1:nm * np
        fprintf(fp, '%f,', statistic_matrix(s, ii) );
    end
    fprintf(fp, '\n');
end
fclose(fp);



% save to csv with problem row wise
savepath = strcat(pwd, '\result_folder\mo_sodsm_igdmedian.csv');
fp=fopen(savepath,'w');
fprintf(fp, 'problem,');
for jj = 1:nm
    fprintf(fp, '%s,', methods{jj}  );
end
fprintf(fp, '\n');


for ii = 1:np
    prob = eval(problems{ii});
     fprintf(fp, '%s,', prob.name);
     for jj = 1:length(methods)
        fprintf(fp, '%f,',  statistic_matrix(3, (ii-1) * length(methods)+jj));
     end
     fprintf(fp, '\n');
end
fclose(fp);
