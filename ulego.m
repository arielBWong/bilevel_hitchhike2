function ulego(prob, seed, num_pop,num_gen,iter_size)
% method of main optimization process of upper level ego
% usage:
%
% input: prob    -problem instance
%--------------------------------------------------------------------------

% performance record
n_feval = 0;

% algo parameter
inisize_l = 30;
numiter_l = 20;
numiter_u = 50;

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
[fu, fc] = pro.evaluate_u(xu, xl);
num_con = size(fc, 2);

%--fu adjust
for i=1:inisize_u
    fu = llfeasi_modify(fu, ll_feasi_flag, i);
end

%-main ulego routine
for i = 1:numiter_u
    %--search next xu
    [newxu, ~] = EIMnext_znorm(xu, fu, upper_bound, lower_bound,num_pop, num_gen, fc);
    %--get its xl
    [newxl, n, flag] = llmatch(newxu, prob,num_pop, num_gen,low_init_size, numiter_l);
    %--evaluate xu
    [newfu, newfc] = prob.evaluate_u(newxu, newxl);
    %--assemble xu fu fc
    xu = [xu; newxu];
    fu = [fu; newfu];
    fc = [fc; newfc];
    llfeasi_flag = [llfeasi_flag, flag];
    %--adjust fu by lower feasibility
    fu = llfeasi_modify(fu, llfeasi_flag, inisize_u+i);  
end

%-bilevel local search
[xu_start, ~, ~] = localsolver_startselection(xu, fc, fu);
 newxu, n_up, n_low = blsovler(prob, xu_start, num_pop, num_gen, inisize_l, numiter_l);
 n_up = n_up + size(xu, 1);
 n_low = n_low + feval;

%-final hybrid ll search
%-- use newxu 
[newxl, feval] = hybrid_llsearch(newxu, prob, hy_pop, hy_gen);
n_low = n_low + feval;

%-performance record
perf_record(prob, newxu, newxl, n_up, n_low);
end

function fu = llfeasi_modify(fu, feasi_list, ind)
% this function is to varify the feasibility of xl
% if xl is infeasible, fu is modified to a higher value
% so that this point is not preferred in later search
% consider range(1:ind) to deal with both one instance and
% a list of instances
if ~feasi_list(ind)
    f = max(fu(1:ind));
    f =  f + 1;
    modify_ind = feasi_list(1:ind) == 0;
    fu(modify_ind) = f;
end
end