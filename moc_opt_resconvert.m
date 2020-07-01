% function moc_opt_resconvert()
% this function create collective results and mean std median
% for 'ZDT1-3', 'DTLZ2/5/7'
% for EIMnext and EIMnext_znorm comparison
% run unit_test5_on_parfor to generate results first
% this function output serve as a support for eim working properly on
% mo problems
%
% -------------------------------------------------------------------------
%%
clearvars;
close all;
eim_methods = {'EIMnext', 'EIMnext_znorm', 'paperdemo'};
test_problems = {'ZDT1','ZDT2','ZDT3'}; %,'DTLZ2','DTLZ5','DTLZ7'};

s = 29;
np = length(test_problems);
ne = length(eim_methods);
out_matrix = zeros(s+3, np*ne);


for j = 1:np % for each problem list two columns
    for i = 1:ne % method output
        k = (j-1) * ne + i;
        % read in each result file
        fn = strcat(pwd, '\result_folder\', eim_methods{i},'_', test_problems{j}, '_hv.csv');
        m = csvread(fn);
        
        out_matrix(1:s, k) = m(1:s);
    end
end

for i=1: np*ne
    
    out_matrix(s+1, i) = mean(out_matrix(1:s, i));
    out_matrix(s+2, i) = std(out_matrix(1:s, i));
    out_matrix(s+3, i) = median(out_matrix(1:s, i));
% no more drawing
%     % draw a graph?
%     col = out_matrix(1:s, i);
%     [~, I] = sort(col);
%     plot_index = I(round(s/2));
%     plot_finalnd(plot_index, i, test_problems, eim_methods)
% 
%     plot_seed = I(round(s/2));
%     plot_finalnd(plot_seed, i, test_problems, eim_methods)

end
output_path = strcat(pwd, '\result_folder\moc_opt_resconvert.csv');

% put header and index when write in file
fp=fopen(output_path,'w');
fprintf(fp, 'seed,');
for i = 1:np
    header = strcat(test_problems{i},'_', eim_methods{1});
    fprintf(fp, '%s,', header);
    % header = strcat(test_problems{i},'_', eim_methods{2});
    % fprintf(fp, '%s,', header);
end
fprintf(fp, '\n');

% format matrix
for i = 1:s
    fprintf(fp, '%d,', i);
    for j = 1: ne*np
        fprintf(fp, '%.4f,', out_matrix(i,j));
    end
    fprintf(fp, '\n');
end

% format statistics
st={'mean', 'std', 'median'};
for i = s+1:s+3
    fprintf(fp, '%s,', st{i-s});
    for j = 1: ne*np
        fprintf(fp, '%.4f,', out_matrix(i,j));
    end
    fprintf(fp, '\n');
end
fclose(fp);


%--------------------------------------
% not using the following function
function plot_finalnd(seed, problem_i, test_problems, methods)
np = length(test_problems);
ne = length(methods);

% p =  round(problem_i/ne); % problem index
% e = mod(problem_i, ne); % method index

problem = test_problems{problem_i};
method = methods{1};


fn = strcat(pwd, '\result_folder\', method,'_', problem, '_',num2str(seed),'_trainy.csv');
train_y = csvread(fn);
% fn = strcat(pwd, '\result_folder\', method,'_', problem, '_',num2str(seed), '_trainc.csv');
% train_c = csvread(fn);
num_con = size(train_y, 2);


% index_c = sum(train_c <= 0, 2) == num_con;
% if sum(index_c) ~=0
%     feasible_y = train_y(index_c, :);
%     nd_index = Paretoset(feasible_y);
%     nd_front = feasible_y(nd_index, :);
% else
% nd_front = [];
% end

nd_index = Paretoset(train_y);
nd_front = train_y(nd_index, :);
    

fn = strcat(pwd, '\result_folder\', method,'_', problem,'_', num2str(seed), '_trainy.csv');
train_y = csvread(fn);
fn = strcat(pwd, '\result_folder\', method,'_', problem,'_', num2str(seed), '_trainc.csv');
train_c = csvread(fn);
num_con = size(train_y, 2);

index_c = sum(train_c <= 0, 2) == num_con;
if sum(index_c) ~=0
    feasible_y = train_y(index_c, :);
    nd_index = Paretoset(feasible_y);
    nd_front = feasible_y(nd_index, :);
else
    nd_front = [];
end

fig = figure(1);
scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled');
title(sprintf('%s problem with %s, seed: %d',fun_name,ite, problem, seed));
savename = strcat(pwd, '\result_folder\', method,'_', problem, '_median.jpeg');
saveas(fig, savename);
close(fig);

end
