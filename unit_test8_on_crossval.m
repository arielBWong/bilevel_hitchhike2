%% this test is to test dace training on crossvalidation
% this cross validation should have better accuracy
% cv method also has close point elimination, so it should have better
% modelling ability than non point elimination version
% tests should  support above  two claims
%------
% (1) test on normal and cross validation this first should be single level
% problem, use SMD problems upper level (single objective), and fix xl to optimum
% ** Conclusion: crossvalidation improves performance when training data
% ** exceed certain amount.
% ** if training data is not enough, then crossvalidation is not as good as
% ** using all data for training
%------
% (2) test on point elimination, with point elimination the 
% algorithm is supposed to overcome situations when 
%  ulego has too close point in training data
%  kriging approximation is not able to approximate function
%  to validate above claim, I need to reproduce the data 
%  where prediction on test data are better than only use whole data
% can also be used to compare with crossvalidation train
%-------------------------------------------------------



%----verbose
seed = 2;
rng(seed, 'twister');
problem_folder = strcat(pwd,'\problems\SMD');
addpath(problem_folder);

test1();

%----verbose
rmpath(problem_folder);


function test1()
%--test1
xl  =  [0,  0 , 0];
num_train = 80;
num_test = 30;
prob = smd3();

upper_bound = prob.xu_bu;
lower_bound = prob.xu_bl;
nvar = prob.n_uvar;


train_xu = lhsdesign(num_train, nvar,'criterion','maximin','iterations',1000);
train_xu = repmat(lower_bound, num_train, 1) ...
    + repmat((upper_bound - lower_bound), num_train, 1) .* train_xu;
train_xl = repmat( xl, num_train, 1);
[train_yu, ~] = prob.evaluate_u(train_xu, train_xl );
nobj = size(train_yu, 2);


test_xu = lhsdesign(num_test, nvar,'criterion','maximin','iterations',1000);
test_xu = repmat(lower_bound, num_test, 1) ...
    + repmat((upper_bound - lower_bound), num_test, 1) .* test_xu;
test_xl = repmat( xl, num_test, 1);
[test_yu, ~] = prob.evaluate_u(test_xu, test_xl );

%----normal dace process, test on 1st obj
num_vari = nvar;
for ii = 1:1
    kriging_obj = dacefit(train_xu,train_yu(:,ii),...
        'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
end
%----cv dace, test on 1st obj
for ii = 1:1
    cvkriging_obj = cvtrain(train_xu, train_yu(:,ii));
end

% test part
test_y = dace_predict(test_xu, kriging_obj);
test_cvy = dace_predict(test_xu, cvkriging_obj);

m1 = mse(test_y,test_yu );
disp(m1);
mcv = mse(test_cvy, test_yu);
disp(mcv);
end



function krg = cvtrain(trainx, trainy)
k=5;
num_vari =  size(trainx,2);

% too lazy to change
trainx_r = trainx;
trainy_r = trainy;

% --create krg with cross validation
cv_par = cvpartition(trainy_r, 'KFold',  5);
mse_perfolder = zeros(cv_par.NumTestSets,1);
krg_perfolder  =  cell(1, k);
for j = 1:k
    %---train
    trIdx = cv_par.training(j);
    teIdx = cv_par.test(j);
    krg_perfolder {j}= dacefit(trainx_r(trIdx , :), trainy_r(trIdx, 1),...
        'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    %--test
    [test_y, mse_test]  = dace_predict(trainx_r(teIdx, :), krg_perfolder {j});
    mse_perfolder(j) =  mse(test_y, trainy_r(teIdx, 1));
end

% --pick
[msemin, mseind] = min(mse_perfolder);
krg = krg_perfolder{mseind};


end





