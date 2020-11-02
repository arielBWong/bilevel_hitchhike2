%% profiling
clearvars;
close all;


savepath = strcat(pwd, '\result_folder\matchxl_invest\');
outputname = 'singleObjCon_median.csv';

result_file = strcat(savepath, outputname);
leg = { 'HYB','EIM', 'BLE'};

mat = csvread(result_file, 1, 1);
median_out = mat(:,[1,3,5]);
min_s = min(median_out, [], 2);


np  = size(median_out, 1);
nm  = 3;

%
n = 10;
tao_list = linspace(1, 2, n);

% plot rows
plotrows = zeros(nm, n);

for ii = 1:n
    for jj = 1:nm
        count = 0;
        for kk = 1:np
            % ---------
            performance_p = median_out(kk, jj);
            r = performance_p/min_s(kk);
            if r < tao_list(ii)
                count = count + 1;
            end
            
        end
        
        plotrows(jj, ii) = count/np;
    end
end

for ii = 1:nm
    plot(tao_list, plotrows(ii, :)); hold on;
end

 legend( leg{1}, leg{2}, leg{3},'FontSize', 14);
 title('profile');
 xlim([1, 2]);

a = 0;



