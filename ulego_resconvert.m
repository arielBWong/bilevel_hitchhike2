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

        
s = 29;
seeds = linspace(1, s, s);


np = length(problems);
ns = length(seeds);
collectmatrix = zeros(ns + 5, np * 2);
feasimatrix = zeros(ns, np * 2);
nfevmatrix = zeros(ns, np*2);

for i = 1:np
    savepath = strcat(pwd, '\result_folder\', problems{i}(1:end-2) );
    for j = 1:ns
        seed = seeds(j);
        singlerun_file = strcat(savepath, '\out_', num2str(seed),'.csv');
        smatrix = csvread(singlerun_file);                % performance records
        collectmatrix(j, (i-1)*2 + 1) = smatrix(1,1);     % upper accuracy
        collectmatrix(j, (i-1)*2 + 2) = smatrix(1,2);     % lower accuracy
        feasimatrix(j, (i-1)*2 + 1)  = smatrix(3, 1);     % upper feasi
        feasimatrix(j, (i-1)*2 + 2) = smatrix(3, 2);      % lower feasi
        nfevmatrix(j, (i-1)*2 + 1) = smatrix(2, 1);       % number of function evaluation
        nfevmatrix(j, (i-1)*2 + 2) = smatrix(2, 2);       % number of function evaluation
    end
end
for i =1: np
   %--upper level statistics
   collectmatrix(ns+1, (i-1)*2 + 1) = mean(collectmatrix(1:ns,   (i-1)*2 + 1));
   collectmatrix(ns+2, (i-1)*2 + 1) = std(collectmatrix(1:ns,    (i-1)*2 + 1));
   collectmatrix(ns+3, (i-1)*2 + 1) = median(collectmatrix(1:ns, (i-1)*2 + 1));
   
   %--lower level statistics
   collectmatrix(ns+1, (i-1)*2 + 2) = mean(collectmatrix(1:ns,  (i-1)*2 + 2));
   collectmatrix(ns+2, (i-1)*2 + 2) = std(collectmatrix(1:ns,   (i-1)*2 + 2));
   
   %--lower median, use corresponding upper median
   upper = collectmatrix(1:ns,   (i-1)*2 + 1);
   lower = collectmatrix(1:ns,   (i-1)*2 + 2);
   
   [sortedupper, id] = sort(upper);                     % upper median seed
   medianind = id(round(s/2));          
   
   collectmatrix(ns+3, (i-1)*2 + 2) = lower(medianind); % extract lower
   
   %-- number of function valuations
   collectmatrix(ns+4, (i-1)*2 + 1) = nfevmatrix(medianind, (i-1)*2 + 1);
   collectmatrix(ns+4, (i-1)*2 + 2) = nfevmatrix(medianind, (i-1)*2 + 2);
   
   %--feasibility
   collectmatrix(ns+5, (i-1)*2 + 1) = feasimatrix(medianind, (i-1)*2 + 1);
   collectmatrix(ns+5, (i-1)*2 + 2) = feasimatrix(medianind, (i-1)*2 + 2);
   
end


output_path = strcat(pwd, '\result_folder\smd_resconvert.csv');
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
    for j = 1: np * 2
        fprintf(fp, '%.4f,', collectmatrix(i,j));
    end
    fprintf(fp, '\n');
end
%--format statistics mean/std/median
st={'mean', 'std', 'median', 'nf', 'feasibility'};
for i = ns+1:ns+5
    fprintf(fp, '%s,', st{i-ns});
    for j = 1: np*2
        fprintf(fp, '%.4f,', collectmatrix(i,j));
    end
    fprintf(fp, '\n');
end


fclose(fp);


%-----
%rmpath(problem_folder);