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
eim_methods = {'EIMnext_znorm' ,'EIMnext'};
test_problems = {'ZDT1()','ZDT2()','ZDT3()','DTLZ2()','DTLZ5()','DTLZ7()'};

np = length(test_problems);
ne = length(eim_methods);
out_matrix = zeros(30+3, np*ne);


for j = 1:np % for each problem list two columns
    for i = 1:ne % method output
        k = (j-1) * ne + i;
        % read in each result file
        fn = strcat(pwd, '\result_folder', eim_methods{i},'_', test_problems{j}, '_hv.csv');
        m = csvread(fn);
        
        out_matrix(1:30, k) = m;
    end
end

for i=1: np*ne
    
    out_matrix(31, i) = mean(out_matrix(1:30, i));
    out_matrix(32, i) = std(out_matrix(1:30, i));
    out_matrix(33, i) = median(out_matrix(1:30, i));
    
end
output_path = strcat(pwd, '\result_folder\moc_opt_resconvert.csv');

% put header and index when write in file
fp=fopen(output_path,'w');
fprintf(fp, 'seed,');
for i = 1:np  
    header = strcat(test_problems{i},'_', eim_methods{1});
    fprintf(fp, '%s,', header);
    header = strcat(test_problems{i},'_', eim_methods{2});
    fprintf(fp, '%s,', header);
end
fprintf(fp, '\n');

% format matrix
for i = 1:30
    fprintf(fp, '%d,', i);
    for j = 1: ne*np
        fprintf(fp, '%.4f,', out_matrix(i,j));
    end
    fprintf(fp, '\n');
end

% format statistics
st={'mean', 'std', 'median'};
for i = 31:33
    fprintf(fp, '%s,', st{i-30});
    for j = 1: ne*np
        fprintf(fp, '%.4f,', out_matrix(i,j));
    end
    fprintf(fp, '\n');
end
fclose(fp);
