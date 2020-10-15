function ulego_corehyb(prob, seed, alg_next, fitnesshandle, normhn, coresteps)
% method of main optimization process of upper level ego
% adapt to upper level problems of "single objective"
% usage:
%     input
%               prob                          : problem instance
%               seed                          : random process seed
%               alg_next                            : string, the name of eim function
%               fitnesshandle          : string, the function that eim use to
%                                                                evaluate fitess in its ea process of proposing next point
%               normhn                    : string, normalization function used in EIMnext_znorm
%
%     output
%               csv files saved in result folder
%               performance statistics include 3*3 matrix
%                                                       [  ul  accuracy,                       ll accuracy;
%                                                          upper number of feval,   lower number of feval;
%                                                          upper feasibility,               lower feasibility]
%
%
%--------------------------------------------------------------------------
tic;
rng(seed, 'twister');
% performance record variable
n_feval = 0;




% algo parameter
inisize_l = 10;
numiter_l = 30;  % lower dual iteration
numiter_u = 40;
inisize_u = 10;
num_pop   = 20;
num_gen   = 20;


% parallel compatible setting
prob = eval(prob);
numl = prob.n_lvar;
% 
% savepath = strcat(pwd,  '\result_folder\', prob.name, '_', num2str(numl),'_hyb');
% file = strcat(savepath, '\out_', num2str(seed),'.csv');
% if exist(file,'file') == 2  % ignore existing runs
%     fprintf('skip');
%     return;
% end

eim    = str2func(alg_next);
fithn  = str2func(fitnesshandle);
normhn = str2func(normhn);

%--upper problem variable
u_nvar      = prob.n_uvar;
upper_bound = prob.xu_bu;
lower_bound = prob.xu_bl;


%-xu initialization
xu = lhsdesign(inisize_u,u_nvar,'criterion','maximin','iterations',1000);
xu = repmat(lower_bound, inisize_u, 1) + repmat((upper_bound - lower_bound), inisize_u, 1) .* xu;

xl = [];
llfeasi_flag = [];
%-xu match its xl and evaluate fu
sname = prob.name;
fprintf('init problem %s, seed %d', sname, seed);
for i=1:inisize_u
    disp(i);
    [xl_single, n, flag] = llmatch_hyb(xu(i, :), prob,num_pop, num_gen,alg_next, numiter_l, fitnesshandle, seed);
    xl = [xl; xl_single];
    llfeasi_flag = [llfeasi_flag, flag];
    n_feval = n_feval + n; %record lowerlevel nfeval
end

%--xu evaluation
[fu, fc] = prob.evaluate_u(xu, xl);
[lu, lc] = prob.evaluate_l(xu, xl);
num_con  = size(fc, 2);
% scatter(xu, fu, 'b'); drawnow;

%--fu adjust
for i=1:inisize_u
    fu = llfeasi_modify(fu, llfeasi_flag, i);
end
%-disp('main ego')
%-main ulego routine
for i = 1:numiter_u
    %--search next xu
    
    [newxu, info] = eim(xu, fu, upper_bound, lower_bound,num_pop, num_gen, fc, fithn, normhn);
    %--get its xl
    [newxl, n, flag] = llmatch_hyb(newxu, prob,num_pop, num_gen,alg_next, numiter_l, fitnesshandle, seed);
    llfeasi_flag = [llfeasi_flag, flag];
    n_feval = n_feval + n;  
    [xu, xl, fu, lu, fc, lc] = postnew_process(prob, newxu, newxl, xu, xl, fu, lu, fc, lc, llfeasi_flag);
   
    %-----------------------------------
    %--dual process of believer
    [krg_obj, krg_con, ~] = update_surrogate(xu, fu, fc, str2func('normalization_z'));
    funh_obj = @(x)ulobj(x, krg_obj);
    funh_con = @(x)ulcon(x, krg_con);
    param.gen     = num_gen;
    param.popsize = num_pop;
    [~,~,~, archive] = gsolver(funh_obj, prob.n_uvar,  prob.xu_bl, prob.xu_bu, [], funh_con, param);
    
    [newxu, ~]     = believer_select(archive.pop_last.X, xu, prob, true);
    [newxl, n, flag] = llmatch_hyb(newxu, prob, num_pop, num_gen, alg_next, numiter_l, fitnesshandle, seed);
    llfeasi_flag = [llfeasi_flag, flag];
    n_feval = n_feval + n;
    
    [xu, xl, fu, lu, fc, lc] = postnew_process(prob, newxu, newxl, xu, xl, fu, lu, fc, lc, llfeasi_flag);
    
end
if coresteps % no local search
    n_up =  size(xu, 1);
    n_low = n_feval;
    method = strcat('hyb_init', num2str(inisize_u));
    ulego_coreending(xu, fu, fc, xl, prob, seed, n_up, n_low, method);
else
    upper_localpostprocess(prob, xu, xl, fu, n_feval, seed, 'hyb');
end
end


function [xu, xl, fu, lu, fc, lc] = postnew_process(prob, newxu, newxl, xu, xl, fu, lu, fc, lc, llfeasi_flag)
% one by one post process with new xu and xl

%--evaluate xu
[newfu, newfc] = prob.evaluate_u(newxu, newxl);
[newlu, newlc] = prob.evaluate_l(newxu, newxl);

%--assemble xu fu fc
xu = [xu; newxu];
xl = [xl; newxl];
fu = [fu; newfu];    lu = [lu; newlu];
fc = [fc; newfc];    lc = [lc; newlc];

%--adjust fu by lower feasibility
checkindex = size(xu, 1);
fu = llfeasi_modify(fu, llfeasi_flag, checkindex);  % --?
end


%----------------------------------------------------
%-- surrogate objective --
function  f = ulobj(x, kriging_obj)
num_obj = length(kriging_obj);   % krg cell array?
num_x = size(x, 1);
f = zeros(num_x, num_obj);
for ii =1:num_obj
    [f(:, ii), ~] = dace_predict(x, kriging_obj{ii});
end
end

function c = ulcon(x, krging_con)
num_con = length(krging_con);
num_x = size(x, 1);
c = zeros(num_x, num_con);
for ii =1:num_con
    [c(:, ii), ~] = dace_predict(x, krging_con{ii});
end
end




