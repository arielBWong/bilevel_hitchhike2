function ulego_coregen(prob_str, seed, normhn, coresteps)
% this function implements the sabla like ego method
% instead of one by one propose next x
% a whole population is evaluated periodically
% use population based krg to accumulate xu and return nd front from xu
% input
%   prob_str: problem instance
%   seed: seed
% output
%   saved nd front (xu, fu, fc) in csv file and number of function
%   evaluation
% main steps:
% (1) initialize a population of xu
% (2) find ll xl match
% (3) train krg for ul problem
% (4) use krg and ea to search xu that maximize predicted fu
% (5) after iter_frequ generations, evaluate population with real problem
% evaluation. update krg and continue evolving on krg
%--------------------------------

rng(seed, 'twister');
% algorithm parameter

prob = eval(prob_str);
numl = prob.n_lvar;
% save some runs
savepath = strcat(pwd, '\result_folder\', prob.name, prob.name, '_', num2str(numl), '_gen_addon');
file = strcat(savepath, '\out_', num2str(seed),'.csv');
if exist(file,'file') == 2  % ignore existing runs
    return;
end


num_popu   =   20;   % 80 in total
num_genu   =   60;
iter_frequ =   20;

num_popl   =   20;   % 60 in total
num_genl   =   40;
iter_freql =   20;
evaln      =   num_popu;
max_nl     =   20000;   % control on max
%--------
%
normhn= str2func(normhn);
n_feval = 0;
%--upper problem variable
u_nvar = prob.n_uvar;
upper_bound = prob.xu_bu;
lower_bound = prob.xu_bl;

%--xu initialization
xu = lhsdesign(num_popu,u_nvar,'criterion','maximin','iterations',1000);
xu = repmat(lower_bound, num_popu, 1) + repmat((upper_bound - lower_bound), num_popu, 1) .* xu;

xl = [];
llfeasi_flag = [];
% - xu match its xl and evaluate fu
for i=1:num_popu
    [xl_single, nl, flag]  = llmatch_sao_pop(xu(i, :), prob, num_popl, num_genl, iter_freql);
    xl = [xl; xl_single];
    llfeasi_flag = [llfeasi_flag, flag];
    n_feval = n_feval + nl;           %record lowerlevel nfeval
end

%--xu evaluation
[fu, fc] = prob.evaluate_u(xu, xl);
num_con = size(fc, 2);

%--fu adjust
for i=1:num_popu
    fu = llfeasi_modify(fu, llfeasi_flag, i);
end

% -- create surrogate for first round evolution
[krg_obj, krg_con, ~] = update_surrogate(xu, fu, fc, normhn);

%--main population based optimization  routine
n = round(num_genu/iter_frequ);       % how many  evolution  subroutine
initmatrix = xu;                      % no normalization on x
for g=1:n
    disp(g);
    % first subroutine has first population evaluated, and last generation
    % evaluated
    % (krg prediction is exactly corresponding f value)
    if g==1
        param.gen=iter_frequ;
    else
        % for other subroutines, its previous subroutine's last population is
        % evaluated, so for every freq(e.g. 10) generations, evaluate once, it means
        % gsolver needs evolve freq + 1 times
        param.gen=iter_frequ + 1;
    end
    param.popsize= num_popu;
    
    funh_obj = @(x)ulobj(x, krg_obj);
    funh_con = @(x)ulcon(x, krg_con);
    
    [~,~,~, archive] = gsolver(funh_obj, u_nvar,  prob.xu_bl, prob.xu_bu, initmatrix, funh_con, param);
    
    new_xu = archive.pop_last.X(1:evaln, :);
    % replace evaluate whole population with evaluate unseen data
    repeat_index = ismember(new_xu, xu, 'row');
    new_index = ~repeat_index;
    num_new =  sum(new_index);
    new_xl  = [];
    new_xu = new_xu(new_index, :);
    
    % if new point does not exist
    % continue evolution with re-start evolution with random initialization
    num_add = param.popsize - num_new;
    if num_add > 0
        fprintf('compensate unseen data');
        add_xu = add_randompoints(prob.xu_bl,prob.xu_bu, prob.n_uvar, xu, num_add);
        new_xu = [new_xu; add_xu];
    end
 
    % --find xl for new_xu
    num_new = size(new_xu, 1);
    if num_new ~= param.popsize
        error('true evaluation size is not population size');
    end
    for i=1:num_new
        % match new_xl for new_xu
        [xl_single, nl, flag]  = llmatch_sao_pop(new_xu(i, :), prob, num_popl, num_genl, iter_freql);
        new_xl = [new_xl; xl_single];
        llfeasi_flag = [llfeasi_flag, flag];
        n_feval = n_feval + nl;                         %record lowerlevel nfeval
    end
    [newfu, newfc] =  prob.evaluate_u(new_xu, new_xl);
    
    % add to training
    xu = [xu; new_xu];  xl = [xl; new_xl];
    fu = [fu; newfu];   fc = [fc; newfc];
    
    % adjust fu according to lower flag
    tr_size = size(xu, 1);
    for i = tr_size - num_new + 1: tr_size
        fu = llfeasi_modify(fu, llfeasi_flag, i);
    end
    
    [krg_obj, krg_con, ~] = update_surrogate(xu, fu, fc, normhn);
    % initmatrix = new_xu;
    % initmatrix = archive.pop_last.X;
    initmatrix = [];
end

if coresteps
    n_up =  size(xu, 1);
    n_low = n_feval;
    ulego_coreending(xu, fu, fc, xl, prob, seed, n_up, n_low, 'eim');
else
    upper_localpostprocess(prob, xu, xl, fu, n_feval , seed, 'gen');
   
end
end



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
