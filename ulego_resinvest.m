clearvars;
close all;

seedmax = 11;
          
problems ={'smd9()'};     
prob = eval(problems{1});
                       
methods = {'eim',  'bel'};  % 'llmatchpop',
leg = {'EIM', 'BLE'};
color = {'r', 'b'};
seeds = [1, 1];
np= length(methods);
nm = length(methods);
    fig1 = gcf;
for i = 1:np
    savepath = strcat(pwd, '\result_folder\', prob.name ,'_', methods{i}, '\xu_', num2str(seeds(i)), '.csv')
    xu = csvread(savepath);
    
    scatter(xu(:, 1), xu(:,2), 40, color{i}, 'filled'); hold on;
    
    
end