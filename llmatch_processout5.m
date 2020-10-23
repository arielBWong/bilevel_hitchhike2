%% evaluate 100 median prediction and real function
% multiple xu-xl

clearvars;
close all;

seedmax = 11;
problems = { 'dsm1(5,5)'};


% methods = {'llmatcheim',  'llmatchble',  'llmatchpop'};  % 'llmatchpop',
% leg = {'EIM', 'BEL', 'GEN'};

leg = {'EIM','BEL', 'HYB'};
np  = length(problems);
nm  = length(methods);

k = 5;
init_size = 11 * k - 1;
methods = { 'llmatcheim',  'llmatchble'};

for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
     for jj = 1:nm
        method = methods{jj};
        num = length(prob.xl_bl);
        %savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_llstudy_',method);
        savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_',method);
        savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_',method, '_init_', num2str(init_size));
     end
end
    
%-----auxiliary function ---
function [krg_obj, krg_con, info] = update_surrogate(trainx, trainy, trainc)
% this function updates train kriging model from train x and train y
%
%
if size(trainy,2)>1
    error(' following zscore norm only applies to single objective problems');
end
[train_y_norm, y_mean, y_std] = zscore(trainy, 1, 1);
num_obj = size(trainy, 2);
krg_obj = cell(1, num_obj);
num_vari = size(trainx, 2);
for ii = 1:num_obj
    % kriging_obj{ii} = dace_train(train_x_norm,train_y_norm(:,ii));
    krg_obj{ii} = dacefit(trainx,train_y_norm(:,ii),...
        'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  % for test
end

info = struct();
info.ymean = y_mean;
info.ystd = y_std;

% deal with constraints
if ~isempty(trainc)
    num_con = size(trainc, 2);
    krg_con = cell(1, num_con);
    [train_c_norm, c_mean, c_std] = zscore(trainc, 1, 1);
    
    for ii = 1:num_con
        krg_con{ii} = dacefit(trainx, train_c_norm(:,ii),...
            'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  % for test
    end
    info.cmean = c_mean;
    info.cstd = c_std;
else
    krg_con = [];
    info.cmean = [];
    info.cstd = [];
end

end


% believer objectives 
function  f = llobj(x, kriging_obj)
num_obj = length(kriging_obj);   % krg cell array?
num_x = size(x, 1);
f = zeros(num_x, num_obj);
for ii =1:num_obj
    [f(:, ii), ~] = dace_predict(x, kriging_obj{ii});
end
end


% believer constraints
function c = llcon(x, krging_con)
num_con = length(krging_con);
num_x = size(x, 1);
c = zeros(num_x, num_con);
for ii =1:num_con 
    [c(:, ii), ~] = dace_predict(x, krging_con{ii});
end
end




