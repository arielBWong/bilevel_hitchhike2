%% this method process experiment results 
%
clearvars;
close all;


seed = 1;
% read a seed and 
problems = { 'mobp5()', 'mobp7()','mobp8()','mobp9(6)','mobp10()','mobp11(6)' };
methods = {'Ehv_eval', 'EIM_eval'};

ii = 1;
jj = 1;

prob = eval(problems{ii});
method = methods{jj};

savepath = strcat(pwd, '\result_folder\', prob.name, '_', method);
savename_fu = strcat(savepath, '\fu_', num2str(seed),'.csv');

nd_front = csvread(savename_fu);
fig1 = gcf;
scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
t = [prob.name,' ', method, ' seed ', num2str(seed) ];
title(t);
savename = strcat(savepath, '\nd_', num2str(seed),'.fig');
savefig(savename);
savename = strcat(savepath, '\nd_', num2str(seed),'.png');

saveas(fig1, savename);
close(fig1);
