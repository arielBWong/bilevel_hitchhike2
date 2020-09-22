
%-----auxiliary function ---
function [krg_obj, krg_con, info] = update_surrogate(trainx, trainy, trainc, normhn)
% this function updates train kriging model from train x and train y
%
%
train_y_norm = normhn(trainy);
num_obj = size(trainy, 2);
krg_obj = cell(1, num_obj);
num_vari = size(trainx, 2);
for ii = 1:num_obj
    % kriging_obj{ii} = dace_train(train_x_norm,train_y_norm(:,ii));
    krg_obj{ii} = dacefit(trainx,train_y_norm(:,ii),...
        'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  % for test
end

info = struct();

% deal with constraints
if ~isempty(trainc)
    num_con = size(trainc, 2);
    krg_con = cell(1, num_con);
    
    % constraints should not be normalized
    for ii = 1:num_con
        krg_con{ii} = dacefit(trainx, trainc(:,ii),...
            'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  % for test
    end
else
    krg_con = [];
end

end