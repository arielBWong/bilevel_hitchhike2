%% this method process experiment results
%
clearvars;
close all;
%problem_folder = strcat(pwd,'\problems\SMD');
%addpath(problem_folder);

% problems = {'smd1()','smd2()','smd3()','smd4()','smd5()','smd6()','smd7()',...
%    'smd8()',};

problems = { 'smd1()','smd2()','smd3()','smd4()','smd5()','smd6()','smd7()',...
    'smd8()','smd9()', 'smd10()','smd11()','smd12()'};


s = 11;
seeds = linspace(1, s, s);
algs = {'eim', 'bel', 'gen'};

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
            singlerun_file = strcat(filename, '\out_', num2str(seed),'.csv')
            smatrix = csvread(singlerun_file);                           % performance records
            accuracyupper(j, (i-1)*na + k) = smatrix(1,1);     % upper accuracy
            accuracylower(j, (i-1)*na + k) = smatrix(1,2);     % lower accuracy
            feasiupper(j, (i-1)*na + k)    = smatrix(3, 1);        % upper feasi
            feasilower(j, (i-1)*na + k)   = smatrix(3, 2);             % lower feasi
        end
    end
end


output_path = strcat(pwd, '\result_folder\smd_convertmean.csv');
fp=fopen(output_path,'w');
fprintf(fp, 'problem,');
for i = 1:na
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
        index = (i-1) * na + jj;
        
        % calculate upper level mean accuracy
        mea =  mean( accuracyupper(:, index ));
        st= std( accuracyupper(:, index ));
        item = strcat(num2str(mea), char(177), num2str(st));
        fprintf(fp, '%s,', item);
        
        % calculate lower level mean accuracy
        mea =  mean( accuracylower(:, index ));
        st = std( accuracylower(:, index ));
        item = strcat(num2str(mea), char(177), num2str(st));
        fprintf(fp, '%s,', item);
    end
    fprintf(fp, '\n');
end
fclose(fp);


output_path = strcat(pwd, '\result_folder\smd_convertmedian.csv');
fp=fopen(output_path,'w');
fprintf(fp, 'problem,');
for i = 1:na

    
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
    prob =  eval(problems{ii});
    fprintf(fp, '%s, ',  prob.name);
    
    middle = round(s/2);
    for jj = 1:na
        index =  (ii-1) * na + jj;
        
        % calculate upper level median accuracy
        one_column = accuracyupper(:, index);
        [~, ids] = sort(one_column);
        med = one_column(ids(middle));
        medseed = ids(middle);
        fprintf(fp, '%s,', num2str(med));
        
        % extract corresponding lower level accuracy
        one_column = accuracylower(:, index);
        med = one_column(medseed);
        fprintf(fp, '%s,', num2str(med));
        
        %  feasibility ul
        one_column = feasiupper(:, index);
        feasi = one_column(medseed);
        fprintf(fp, '%s,', num2str(feasi));
        
        %  feasibility ll
        one_column = feasilower(:, index);
        feasi = one_column(medseed);
        fprintf(fp, '%s,', num2str(feasi));
        
        % seed
        
        fprintf(fp, '%s,', num2str(medseed));
    end
    fprintf(fp, '\n');
end
fclose(fp);

%-----
%rmpath(problem_folder);