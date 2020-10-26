function llmatch_behaviourstudy(prob, match_method, seed, lower_init)
% this function takes xu and problem as input
% go through different lower level match method
% return global match accuracy
% save to certain files
% input

prob = eval(prob);
% xu = select_xu(prob, 1);
xu = [1, 1];
match_check(match_method, xu, lower_init, prob, seed);

end

function xu = select_xu(prob, samplesize)
rng(999, 'twister');
u_nvar = prob.n_uvar;

upper_bound = prob.xu_bu;
lower_bound = prob.xu_bl;



xu = lhsdesign(samplesize,u_nvar,'criterion','maximin','iterations',1000);
xu = repmat(lower_bound, samplesize, 1) ...
    + repmat((upper_bound - lower_bound), samplesize, 1) .* xu;


end

function match_check(method, xu, lower_init, prob, seed)
samplesize = size(xu, 1);
xl = [];

lower_iter =300; 
normalization = 'normalization_nd';
% normalization = 'normalization_z';

if ~contains(method, '_')
    rng(seed, 'twister');
    match_method =  str2func(method);
    tic;
    for i = 1:samplesize
        xu_i = xu(i, :);
        [match_xl, n_fev, flag] = match_method(xu_i, prob, 20, 20, 'EIMnext',lower_init, lower_iter, 'EIM_eval', seed, normalization);
        
        xl = [xl; match_xl];
    end
    toc;
elseif contains(method, 'hyb')
    rng(seed, 'twister');
    match_method =  str2func( method);
    for i = 1:samplesize
        xu_i = xu(i, :);
        [match_xl, n_fev, flag] =match_method(xu_i, prob, 20, 20, 'EIMnext', 20, 'EIM_eval', seed);
        xl = [xl; match_xl];
    end
elseif contains(method, 'archive')
    rng(seed, 'twister');
    match_method =  str2func( method);
    freq =20;
    tic;
    for i = 1:samplesize
        xu_i = xu(i, :);
        [match_xl, n_fev, flag] = match_method (xu_i, prob, 20, lower_init, lower_iter, freq, seed, normalization);      
        xl = [xl; match_xl];
    end
    toc;
elseif contains(method, 'adapt')
   rng(seed, 'twister');
   match_method =  str2func( method);
    for i = 1:samplesize
        xu_i = xu(i, :);
        [match_xl, n_fev, flag] =  llmatch_believeradapt(xu, prob, 20, 20, 'EIMnext', lower_init, lower_iter, 'EIM_eval', seed, normalization);
        xl = [xl; match_xl];
    end
   
    
end

% save_study(xu, xl, prob,method, seed);
end

function save_study(xu, xl, prob, method, seed)
num = length(prob.xl_bl);
savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_llstudy_',method);
% savepath = strcat(pwd, '\result_folder\', prob.name, '_',method);
n = exist(savepath);
if n ~= 7
    mkdir(savepath)
end
[fl, cl] = prob.evaluate_l(xu, xl);

savename_xu = strcat(savepath, '\xu_', num2str(seed),'.csv');
savename_xl = strcat(savepath, '\xl_', num2str(seed),'.csv');
savename_fl = strcat(savepath, '\fl_', num2str(seed),'.csv');
savename_cl = strcat(savepath, '\cl_', num2str(seed),'.csv');

csvwrite(savename_xu, xu);
csvwrite(savename_xl, xl);
csvwrite(savename_fl, fl);
csvwrite(savename_cl, cl);

end