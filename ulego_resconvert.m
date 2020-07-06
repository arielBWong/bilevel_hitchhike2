%% this method process experiment results 
%
clearvars;
close all;
problem_folder = strcat(pwd,'\problems\SMD');
addpath(problem_folder);

problems = { 'smd9()', 'smd10()', 'smd11()','smd12()'};
s = 11;
seeds = linspace(1, s, s);


np = length(problems);
ns = length(seeds);
collectmatrix = zeros(ns + 3, 4 * 2);

for i = 1:np
    savepath = strcat(pwd, '\result_folder\', problems{i}(1:end-2) );
    for j = 1:ns
        seed = seeds(j);
        singlerun_file = strcat(savepath, '\acc_con_fea_', num2str(seed),'.csv');
        smatrix = csvread(singlerun_file);
        collectmatrix(j, (i-1)*2 + 1) = smatrix(1,1);
        collectmatrix(j, (i-1)*2 + 2) = smatrix(1,2);       
    end
end
for i =1: np*2
   collectmatrix(ns+1, i) = mean(collectmatrix(1:ns, i));
   collectmatrix(ns+2, i) = std(collectmatrix(1:ns, i));
   collectmatrix(ns+3, i) = median(collectmatrix(1:ns, i));
end

output_path = strcat(pwd, '\result_folder\SMD_resconvert.csv');
fp=fopen(output_path,'w');
fprintf(fp, 'seed,');
for i = 1:np
    header = strcat(problems{i}(1:end-2),'_ul');
    fprintf(fp, '%s,', header);
    header = strcat(problems{i}(1:end-2),'_ll');
    fprintf(fp, '%s,', header);
end
fprintf(fp, '\n');

% format matrix
for i = 1:ns
    fprintf(fp, '%d,', i);
    for j = 1: 4 * 2
        fprintf(fp, '%.4f,', collectmatrix(i,j));
    end
    fprintf(fp, '\n');
end
%--format statistics
st={'mean', 'std', 'median'};
for i = ns+1:ns+3
    fprintf(fp, '%s,', st{i-ns});
    for j = 1: 4*2
        fprintf(fp, '%.4f,', collectmatrix(i,j));
    end
    fprintf(fp, '\n');
end
fclose(fp);


%-----
rmpath(problem_folder);