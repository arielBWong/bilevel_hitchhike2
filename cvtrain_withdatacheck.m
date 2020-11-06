function krg =  cvtrain_withdatacheck(trainx, trainy, nd)
% this function trains one kriging model from trainx and trainy
% trainy only has one column data
% usage  -- to be filled
% input
% ouput
%----------------------------------------------------------------
% trainy is one column
if size(trainy, 2)>1
    error('response variable should have only one variable');
end

mse_crossnd = zeros(1, nd);
krg_crossnd = cell(1, nd);
k=5;
num_vari =  size(trainx,2);

% for i = 1:nd % nd number of digits
% --elimination scheme
% --r potentially reduced data size
% [trainx_r,  trainy_r]  = close_elimination(trainx, trainy, i);
trainx_r = trainx;
trainy_r = trainy;
% --create krg with cross validation
cv_par = cvpartition(trainy_r, 'k', k);     %cv class use respond to partition
mse_perfolder = zeros(cv_par.NumTestSets,1);
krg_perfolder  =  cell(1, k);
for j = 1:k
    %---train
    trIdx = cv_par.training(j);
    teIdx = cv_par.test(j);
    krg_perfolder {j}= dacefit(trainx_r(trIdx , :), trainy_r(trIdx, 1),...
        'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    %---test
    [test_y, ~]  = dace_predict(trainx_r(teIdx, :), krg_perfolder {j});
    mse_perfolder(j)  =  mse(test_y, trainy_r(teIdx, 1));
end
%-- pick krg
[mse_min, ind] = min(mse_perfolder);
mse_crossnd = mse_min;
krg = krg_perfolder{ind};
% end
% -pick across nd
[~, ind] = min(mse_crossnd);

% -- retrain krg use all data of 'ind' (decimal accuracy) eliminated data
% krg = krg_crossnd(ind);  % use (k-1) folds data trained model, accuracy
% might be comprmised when training size is small
% [trainx_r,  trainy_r]  = close_elimination(trainx, trainy, ind);
% trainx_r = trainx;
% trainy_r = trainy;

% krg = dacefit(trainx_r, trainy_r,...
%    'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
end


function [x_after, y_after] = close_elimination(x, y, ndigi)
x = round(x, ndigi);
[x_after, ia] = unique(x, 'row');
y_after = y(ia, :);
end