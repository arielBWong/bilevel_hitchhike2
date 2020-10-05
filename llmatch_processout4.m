%% plot median mse in every 10 new training data
%
% for a seed per problem
% (1) load xl, generate xu
% (2) iterate from  index 20 every 10 more until 60,
%       read xl, match xu, train dace, predict best solution
% (3) save mse to another matrix
% (4)

clearvars;
close all;

seedmax = 12;
problems ={'smd12()'};     

prob = eval(problems{1});                    
methods = {'eim', 'bel'};  % 'llmatchpop',
seed = [9,3];
np= length(problems);
nm = length(methods);

filename = strcat(pwd, '\result_folder\', prob.name ,'_',methods{1});
filename = strcat(filename, '\xu_', num2str(seed(1)), '.csv');
xu = csvread(filename);

filename = strcat(pwd, '\result_folder\', prob.name ,'_',methods{1});
filename = strcat(filename, '\xl_', num2str(seed(1)), '.csv');
xl = csvread(filename);
%-----------------------
[fu, cu] = prob.evaluate_u(xu, xl);
cu(cu < 1e-6) = 0;
numcon = size(cu, 2);
index_c = sum(cu <= 0, 2) == numcon;
xu_f = xu(index_c, :);
xu_c = xu(~index_c, :);

scatter(xu_f(:, 1), xu_f(:, 2), 60, 'ro', 'filled'); hold on;
scatter(xu_c(:, 1), xu_c(:, 2), 60, 'ro'); hold on;

%========================================================

filename = strcat(pwd, '\result_folder\', prob.name ,'_',methods{2});
filename = strcat(filename, '\xu_', num2str(seed(1)), '.csv');
xu = csvread(filename);

filename = strcat(pwd, '\result_folder\', prob.name ,'_',methods{2});
filename = strcat(filename, '\xl_', num2str(seed(1)), '.csv');
xl = csvread(filename);

%-----------------------
[fu, cu] = prob.evaluate_u(xu, xl);
cu(cu < 1e-6) = 0;
numcon = size(cu, 2);
index_c = sum(cu <= 0, 2) == numcon;
xu_f = xu(index_c, :);
xu_c = xu(~index_c, :);

scatter(xu_f(:, 1), xu_f(:, 2), 60, 'bo', 'filled'); hold on;
scatter(xu_c(:, 1), xu_c(:, 2), 60, 'bo'); hold on;
 
xu = [1, 1];
scatter(xu(1), xu(2), 80, 'yd','filled'); hold on;


legend('EIM feasible','EIM infeasible','BEL feasible','BEL infeasible','xu prime','FontSize', 16);
xlabel('xu1', 'FontSize', 16);
ylabel('xu2',  'FontSize', 16);
t = prob.name;
title(t,  'FontSize', 16);




% filename = strcat(pwd, '\result_folder\', prob.name ,'_',methods{2});
% filename = strcat(filename, '\xu_', num2str(seed(2)), '.csv');
% xu = csvread(filename);
% scatter(xu(:, 1), xu(:, 2), 80, 'bo'); 
% legend('EIM','BEL','FontSize', 16);
% xlabel('xu1', 'FontSize', 16);
% ylabel('xu2',  'FontSize', 16);
% t = prob.name;
% title(t,  'FontSize', 16);

