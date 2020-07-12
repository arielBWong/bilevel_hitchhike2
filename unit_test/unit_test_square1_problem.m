%%
% This unittest is for single objective one variable
clearvars;
close all;

% (1) intialize training data
% (2) start loop
%   (2-1) call EIMnext to generate next x  point
%   (2-2) combine new x into training data

workdir = pwd;
idcs = strfind(workdir, '\');
upperfolder = workdir(1: idcs(end)-1);
problem_folder = strcat(upperfolder,'\problems');
addpath(problem_folder);
addpath(upperfolder);


num_samples = 3;
num_vari = 1;
train_x = lhsdesign(num_samples,num_vari); %,'criterion','maximin','iterations',1000);
train_y = testproblem1(train_x);
xu_bound = [1];
xl_bound = [0];
test_x = linspace(0, 1, 100);
test_y = testproblem1(test_x');

figure(1);
for iter=1:10
    
    [newx, info] = EIMnext_znorm(train_x, train_y, xu_bound, xl_bound, 20, 50, []);
    newy = testproblem1(newx);
    train_x = [train_x; newx];
    train_y = [train_y; newy];
    krg = info.krg;
    
    norm_x = (test_x' - info.train_xmean)/info.train_xstd;
    prey = dace_predict(norm_x, krg{1});
    prey = prey * info.train_ystd + info.train_ymean;
    
   
    testplot = plot(test_x, test_y, 'b-');hold on;
    predplot = plot(test_x,prey,'r-');
    trainplot = plot(train_x(1:end-1), train_y(1:end-1), 'bo');
    newpont = plot(train_x(end),  train_y(end), 'go');
    legend('real','predict');
    pause(2)
    
    set(predplot, 'Visible', 'off');
    t = 0;
    
    
    
end


rmpath(problem_folder)