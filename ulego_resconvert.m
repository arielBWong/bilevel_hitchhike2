%% this method process experiment results
%
clearvars;
close all;
%problem_folder = strcat(pwd,'\problems\SMD');
%addpath(problem_folder);

%problems = { 'smd9()', 'smd10()', 'smd11()','smd12()'};

problems = { 'smd1()','smd2()','smd3()','smd4()','smd5()','smd6()','smd7()',...
    'smd8()','smd9()', 'smd10()','smd11()','smd12()',...
    'bltp1()','bltp2()','bltp3()', 'bltp5()','bltp6()','bltp7()',...
    'bltp8()','bltp9()', 'bltp10()','bltp11()'};


s = 11;
seeds = linspace(1, s, s);
algs = {'eim', 'gen', 'bel'};

np = length(problems);
ns = length(seeds);
na = length(algs);
collectmatrix = zeros(ns + 5, np * 2);
feasimatrix = zeros(ns, np * 2);
nfevmatrix = zeros(ns, np*2);

accuracyupper = zeros(ns, np * na);
accuracylower = zeros(ns, np * na);
feasiupper = zeros(ns, np * na);
feasilower = zeros(ns, np * na);

for i = 1:np
    prob = eval(problems{i});
    savepath = strcat(pwd, '\result_folder\', prob.name );
    for k = 1:na
        filename = strcat(savepath, '_',  algs{k});
        for j = 1:ns
            seed = seeds(j);
            singlerun_file = strcat(savepath, '\out_', num2str(seed),'.csv');
            smatrix = csvread(singlerun_file);                           % performance records
            accuracyupper(j, (i-1)*np + k) = smatrix(1,1);     % upper accuracy
            accuracylower(j, (i-1)*np + k) = smatrix(1,2);     % lower accuracy
            feasiupper(j, (i-1)*np + k)    = smatrix(3, 1);        % upper feasi
            feasilower(j, (i-1)*np + k)   = smatrix(3, 2);             % lower feasi
        end
    end
end


output_path = strcat(pwd, '\result_folder\smd_convertmean.csv');
fp=fopen(output_path,'w');
fprintf(fp, 'seed,');
for i = 1:na
    fprintf(fp, '%s,', 'problem');
    method = algs{i};
    header = strcat(method,'_ul');
    fprintf(fp, '%s,', header);
    header = strcat(method,'_ll');
    fprintf(fp, '%s,', header);
end
fprintf(fp, '\n');

% upper level accuracy
for i = 1:np
    prob =  eval(problems{i});
    fprintf(fp, '%s, ',  prob.name);
    
    for jj = 1: na
        index = (i-1) * np + jj;
        
        % calculate upper level mean accuracy
        mean =  mean( accuracyupper(:, index ));
        std = std( accuracyupper(:, index ));
        item = strcat(num2str(mean), char(177), num2str(std));
        fprintf(fp, '%s,', item);
        
        % calculate lower level mean accuracy
        mean =  mean( accuracylower(:, index ));
        std = std( accuracylower(:, index ));
        item = strcat(num2str(mean), char(177), num2str(std));
        fprintf(fp, '%s,', item);
    end
    fprintf(fp, '\n');
end
fclose(fp);


output_path = strcat(pwd, '\result_folder\smd_convertmedian.csv');
fp=fopen(output_path,'w');
fprintf(fp, 'seed,');
for i = 1:na
    fprintf(fp, '%s,', 'problem');
    
    method = algs{i};
    header = strcat(method,'_ul');
    fprintf(fp, '%s,', header);
    header = strcat(method,'_ll');
    fprintf(fp, '%s,', header);
    
    header = strcat(method,'_ulfeasibility');
    fprintf(fp, '%s,', header);
    
    header = strcat(method,'_llfeasibility');
    fprintf(fp, '%s,', header);
    
    header = strcat(method,'_medianseed');
    fprintf(fp, '%s,', header);
end
fprintf(fp, '\n');

%
for ii = 1: np
    prob =  eval(problems{i});
    fprintf(fp, '%s, ',  prob.name);
    
    middle = round(seedmax/2);
    for jj = 1:na
        index =  (i-1) * np + jj;
        
        % calculate upper level median accuracy
        one_column = accuracyupper(:, index);
        [~, ids] = sort(one_column);
        med = one_column(ids(middle));
        fprintf(fp, '%s,', num2str(med));
        
        % extract corresponding lower level accuracy
        one_column = accuracylower(:, index);
        med = one_column(ids(middle));
        fprintf(fp, '%s,', num2str(med));
        
        %  feasibility ul
        one_column = feasiupper(:, index);
        feasi = one_column(ids(middle));
        fprintf(fp, '%s,', num2str(feasi));
        
        %  feasibility ll
        one_column = feasilower(:, index);
        feasi = one_column(ids(middle));
        fprintf(fp, '%s,', num2str(feasi));
    end
    fprintf(fp, '\n');
end
fclose(fp);

%-----
%rmpath(problem_folder);