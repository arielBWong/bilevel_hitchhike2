clear all; close all;clc;
% addpath(genpath(strcat(pwd,strcat(filesep,'dace'))));

%%
gaussianpdf(inf)
function y =  gaussianpdf(x)
  y = 1/sqrt(2*pi) * exp(-x.^2/2);
end



% 
% 
% 
% tic;
% 
% xtrg = readtable('train_x.csv' );
% ytrg = readtable('train_y.csv' );
% ytrg = ytrg(:,1);
% 
% xtrg = xtrg{:,:};
% ytrg = ytrg{:,:};
% x = readtable('test_x.csv' );
% x = x{:,:};
% 
% model = dace_train(xtrg, ytrg);
% y = dace_predict(x, model)