function ulego_umoc(prob, seed, eim)
% method of main optimization process of upper level ego
% adapt to upper level problems of "multiple objectives"
% usage: 
%     input
%               prob                          : problem instance                  
%               seed                          : random process seed
%     output  
%               csv files saved in result folder
%               performance statistics include 3*3 matrix
%                                                       [  ul  accuracy, ll accuracy;
%                                                          upper number of feval, lower number of feval;
%                                                          upper feasibility, lower feasibility]
%                                                                      
%                                                                       
%--------------------------------------------------------------------------
tic;
rng(seed, 'twister');
% performance record variable
n_feval = 0;

% algo parameter
inisize_l                  = 30;
numiter_l               = 20;
numiter_u             = 50;
num_pop              = 100;
num_gen              = 100;
hy_pop                  = 20;
hy_gen                  = 50;

% parallel compatible setting
prob = eval(prob);
eim = str2func(eim);

%--upper problem variable
u_nvar = prob.n_uvar;
upper_bound = prob.xu_bu;
lower_bound = prob.xu_bl;
inisize_u = 11 * u_nvar - 1;

%-xu initialization
xu = lhsdesign(inisize_u,u_nvar,'criterion','maximin','iterations',1000);
xu = repmat(lower_bound, inisize_u, 1) + repmat((upper_bound - lower_bound), inisize_u, 1) .* xu;

xl = [];
llfeasi_flag = [];
% -xu match its xl and evaluate fu
for i=1:inisize_u
    [xl_single, n, flag] = llmatch(xu(i, :), prob,num_pop, num_gen,inisize_l, numiter_l);
    xl = [xl; xl_single];
    llfeasi_flag = [llfeasi_flag, flag];
    n_feval = n_feval + n; %record lowerlevel nfeval
end
%--xu evaluation
[fu, fc] = prob.evaluate_u(xu, xl);
num_con = size(fc, 2);

%--fu adjust
for i=1:inisize_u
    fu = llfeasi_modify(fu, llfeasi_flag, i);
end

%-main ulego routine
for i = 1:numiter_u
    %--search next xu
    [newxu, ~] = eim(xu, fu, upper_bound, lower_bound,num_pop, num_gen, fc);
    %--get its xl
    [newxl, n, flag] = llmatch(newxu, prob,num_pop, num_gen,inisize_l, numiter_l);
    n_feval = n_feval + n;
    %--evaluate xu
    [newfu, newfc] = prob.evaluate_u(newxu, newxl);
    %--assemble xu fu fc
    xu = [xu; newxu];
    fu = [fu; newfu];
    fc = [fc; newfc];
    llfeasi_flag = [llfeasi_flag, flag];
    %--adjust fu by lower feasibility
    disp(i);
    fu = llfeasi_modify(fu, llfeasi_flag, inisize_u+i); 
end

%-bilevel local search
[xu_start, ~, ~, ~] = localsolver_startselection(xu, fu, fc); % --?
[newxu, newxl, n_up, n_low] = blsovler(prob, xu_start, num_pop, num_gen, inisize_l, numiter_l);
n_up = n_up + size(xu, 1);
n_low = n_low + n_feval;


%-final hybrid ll search
%-- use newxu
[newxl, feval, flag] = hybrid_llsearch(newxu, newxl, prob, hy_pop, hy_gen);
n_low = n_low + feval;


%-performance record
%--constraints compatible
[fu, cu] = prob.evaluate_u(newxu, newxl);
[fl, cl] = prob.evaluate_l(newxu, newxl);
num_conu = size(cu, 2);
num_conl = size(cl, 2);
% contraint tolerance adjust
cu(cu < 1e-6) = 0;
cl(cl<1e-6) = 0;
% check feasibility
cu = sum(cu<=0, 2)==num_conu;
cl = sum(cl<=0, 2)==num_conl;
perf_record(prob, fu, cu, fl, cl, n_up, n_low, seed);

toc
end
