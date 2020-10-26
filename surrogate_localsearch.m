function[newx] = surrogate_localsearch(xu, newx, prob, trainx, trainf, trainc, normalization)
%----------------------------
% normalization: str
%
%-------------------------

if size(trainf, 2) > 1
    % no local search for MO
     return
end


normhn = str2func(normalization);
[krg_obj, krg_con, info] = update_surrogate(trainx, trainf, trainc, normhn);

%  [f_predict, ~]= predictor(newx, krg_obj{1})
 funh_obj = @(x)llobj(x, krg_obj);
 funh_con = @(x)llcon(x, krg_con);

opts = optimset('fmincon');
opts.Algorithm = 'sqp';
opts.Display = 'off';
opts.MaxFunctionEvaluations = 100;
[newx, newf, ~, output] = fmincon(funh_obj, newx, [], [],[], [],  ...
    prob.xl_bl, prob.xl_bu, [],opts);

% [f_predict, ~] = predictor(newx, krg_obj{1})


 %prob.evaluate_l(xu, newx)
% newx
% a = 0


end


function  f = llobj(x, kriging_obj)
num_obj = length(kriging_obj);   % krg cell array
num_x = size(x, 1);
f = zeros(num_x, num_obj);
for ii =1:num_obj
    [f(:, ii), ~] = dace_predict(x, kriging_obj{ii});
end
end


function c = llcon(x, krging_con)
num_con = length(krging_con);
num_x = size(x, 1);
c = [];
for ii =1:num_con 
    [c(:, ii), ~] = dace_predict(x, krging_con{ii});
end
end