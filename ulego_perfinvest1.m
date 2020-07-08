%%
% function  ulego_perfinvest1(xu, fu, prob, info)
% this function create a plot about the real problem function and
% the kriging recreated function, info comes from EIMnext_znorm
% so this is not an independent function, follow the main process
% of ulego after step main ego
% sprinkle with xu and fu
% assume eim method is eim_znorm
% usage
%   input
%   output
%-------------------------------------------------------------------
problem_folder = strcat(pwd,'\problems\BLTP');
addpath(problem_folder);
prob = bltp5();

% load saved data
variables = load('vari.mat');
info = load('info.mat');
xu = variables.xu;
xl = variables.xl;
fu = variables.fu;

% generate test data from prob upper range
testx = linspace(prob.xu_bl, prob.xu_bu, 50);
testx = testx';

%-exghausive search recreate upper level
% xl = [];
% llfeasi_flag = [];
% for  i = 1:50
%     xu = testx(i, :);
%     [newxl, n_feval, flag] = hybrid_llsearch(xu, [], prob, 100, 100);
%     llfeasi_flag = [llfeasi_flag, flag];
%     xl = [xl; newxl];
% end
% 
% fu = prob.evaluate_u(testx, xl);
% for i=1:50
%     fu = llfeasi_modify(fu, llfeasi_flag, i);
% end
% plot(testx, fu);
%--------------------------------------------------------


%----------------------------------------------

% znorm with info
testx_norm = (testx - info.train_xmean)./info.train_xstd;

% prediction test data with info.ego
testy_norm = dace_predict(testx_norm, info.krg{1});

% denormalization
testy = testy_norm .* info.train_ystd + info.train_ymean;

% plot in figure
predplot = scatter(testx,testy,'r'); hold on;
pointplot = plot(xu, fu, 'bo');

%end
