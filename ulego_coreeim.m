function ulego_coreeim(prob, seed, alg_next, fitnesshandle, normhn, coresteps)
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
numiter_l = 60;
numiter_u = 80;
inisize_u = 10;
num_pop   = 20;
num_gen   = 20;


% parallel compatible setting
prob = eval(prob);
numl = prob.n_lvar;

% savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(numl),'_eim_addon');
% file = strcat(savepath, '\out_', num2str(seed),'.csv');
% if exist(file,'file') == 2  % ignore existing runs
%     fprintf('skip');
%     return;
% end




eim = str2func(alg_next);
fithn = str2func(fitnesshandle);
normhn = str2func(normhn);

%--upper problem variable
u_nvar = prob.n_uvar;
upper_bound = prob.xu_bu;
lower_bound = prob.xu_bl;


%-xu initialization
xu = lhsdesign(inisize_u,u_nvar,'criterion','maximin','iterations',1000);
xu = repmat(lower_bound, inisize_u, 1) + repmat((upper_bound - lower_bound), inisize_u, 1) .* xu;

xl = [];
llfeasi_flag = [];
% -xu match its xl and evaluate fu
sname = prob.name;
fprintf('init problem %s, seed %d',sname, seed);
for i=1:inisize_u
    disp(i);
    [xl_single, n, flag] = llmatch(xu(i, :), prob,num_pop, num_gen,alg_next, numiter_l, fitnesshandle, seed);
    
    xl = [xl; xl_single];
    llfeasi_flag = [llfeasi_flag, flag];
    n_feval = n_feval + n; %record lowerlevel nfeval
end
%--xu evaluation
[fu, fc] = prob.evaluate_u(xu, xl);
[lu, lc] = prob.evaluate_l(xu, xl);
num_con = size(fc, 2);
% scatter(xu, fu, 'b'); drawnow;

%--fu adjust
for i=1:inisize_u
    fu = llfeasi_modify(fu, llfeasi_flag, i);
end
% disp('main ego')
%-main ulego routine
for i = 1:numiter_u
    %--search next xu
    [newxu, info] = eim(xu, fu, upper_bound, lower_bound,num_pop, num_gen, fc, fithn, normhn);
    %--get its xl
    [newxl, n, flag] = llmatch(newxu, prob,num_pop, num_gen,alg_next, numiter_l, fitnesshandle, seed);
    
    n_feval = n_feval + n;
    %--evaluate xu
    [newfu, newfc] = prob.evaluate_u(newxu, newxl);
    [newlu, newlc] =  prob.evaluate_l(newxu, newxl);
    
    %--assemble xu fu fc
    xu = [xu; newxu];
    xl = [xl; newxl];
    fu = [fu; newfu];    lu = [lu; newlu];
    fc = [fc; newfc];    lc = [lc; newlc];
    
    llfeasi_flag = [llfeasi_flag, flag];
    %--adjust fu by lower feasibility
    fu = llfeasi_modify(fu, llfeasi_flag, inisize_u+i);  % --?
    % scatter(xu, fu, 'r'); drawnow;
    disp(i);
end
if coresteps % no local search
    n_up =  size(xu, 1);
    n_low = n_feval;
    method = strcat('eim_init', num2str(inisize_u));
    ulego_coreending(xu, fu, fc, xl, prob, seed, n_up, n_low, method);
else
    upper_localpostprocess(prob, xu, xl, fu, n_feval, seed, 'eim');
end
end
