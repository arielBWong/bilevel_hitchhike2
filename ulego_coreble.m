function ulego_coreble(prob_str, seed, normhn, coresteps)
% this function implements the sabla like ego method
% instead of one by one propose next x
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

rng(seed);
prob = eval(prob_str);
numl = prob.n_lvar;

% % save some runs
% savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(numl), '_ble_addon');
% filename = strcat(savepath, '\out_', num2str(seed),'.csv');
% if exist(filename,'file') == 2  % ignore existing runs 
%     disp(filename);
%     fprintf('exist');
%     return;
% end

% algorithm parameter
evaln = 1;
n = 80;
initsize_u  = 10;
num_popu    = 20;   % 80 in total
num_genu    = 1200;
iter_frequ  = 20;

num_popl    = 20;   % 60 in total
num_genl    = 800;
iter_freql  = 20;
iter_sizel = 60;   % 70 intotal


max_nl = 20000;
%--------

normhn= str2func(normhn);
n_feval = 0;
%--upper problem variable
u_nvar = prob.n_uvar;
upper_bound = prob.xu_bu;
lower_bound = prob.xu_bl;

%--xu initialization
xu = lhsdesign(initsize_u,u_nvar,'criterion','maximin','iterations',1000);
xu = repmat(lower_bound, initsize_u, 1) + repmat((upper_bound - lower_bound), initsize_u, 1) .* xu;

xl = [];
llfeasi_flag = [];
% -xu match its xl and evaluate fu
for i=1:initsize_u
    [xl_single, nl, flag] = llmatch_sao_archiveinsert(xu(i, :), prob, num_popl, iter_sizel, iter_freql);
     xl = [xl; xl_single];
    llfeasi_flag = [llfeasi_flag, flag];
    n_feval = n_feval + nl;           %record lowerlevel nfeval
end

%--xu evaluation
[fu, fc] = prob.evaluate_u(xu, xl);
num_con = size(fc, 2);

%--fu adjust
for i=1:initsize_u
    fu = llfeasi_modify(fu, llfeasi_flag, i);
end

% -- create surrogate for first round evolution
[krg_obj, krg_con, ~] = update_surrogate(xu, fu, fc, normhn);

%--main population based optimization  routine
                            % how many  evolution  subroutine
% initmatrix = xu; 
initmatrix = [];
% no normalization on x
for g=1:n
    disp(g);
    param.gen = iter_frequ;
    param.popsize = num_popu;
        
    funh_obj = @(x)ulobj(x, krg_obj);
    funh_con = @(x)ulcon(x, krg_con);
    
    [~,~,~, archive] = gsolver(funh_obj, u_nvar,  prob.xu_bl, prob.xu_bu, initmatrix, funh_con, param);
    new_xu = archive.pop_last.X;
    
    % replace with checking whether new xu exists in archive
    repeat_index = ismember(new_xu, xu, 'row');
    new_index = ~repeat_index;
    num_new =  sum(new_index);
    new_xl  = [];
    new_xu = new_xu(new_index, :);
    % if new point does not exist
    % continue evolution with re-start evolution with random initialization
    if num_new ==0
        fprintf('evolution converge and no new point is found, select a random new one \n');
        new_xu = add_randompoints(prob.xu_bl,prob.xu_bu, prob.n_uvar, xu, evaln);
    end
    
    % --add new_xu to xu
    for i=1:evaln
        % match new_xl for new_xu
        % tic;
        [xl_single, nl, flag]  = llmatch_sao_archiveinsert(new_xu(i, :), prob, num_popl, iter_sizel, iter_freql);
        % toc;
        new_xl = [new_xl; xl_single];
        llfeasi_flag = [llfeasi_flag, flag];
        n_feval = n_feval + nl;           %record lowerlevel nfeval
    end
    [newfu, newfc] =  prob.evaluate_u(new_xu(1:evaln, :), new_xl);
    
    % add to training
    xu = [xu; new_xu(1:evaln, :)];  xl = [xl; new_xl];
    fu = [fu; newfu];               fc = [fc; newfc];
    
    % adjust fu according to lower flag
    tr_size = size(xu, 1);
    for i = tr_size - evaln + 1: tr_size
        fu = llfeasi_modify(fu, llfeasi_flag, i);
    end
    
    % check under level number of evaluation
%     if n_feval > max_nl
%         break;
%     end
    % update krg  and initmatrix, continue to evolve
    [krg_obj, krg_con, ~] = update_surrogate(xu, fu, fc, normhn);
    initmatrix = [];
end

if coresteps
    n_up =  size(xu, 1);
    n_low = n_feval;
    method = strcat('ble_init', num2str(initsize_u));
    ulego_coreending(xu, fu, fc, xl, prob, seed, n_up, n_low, method);
else
    upper_localpostprocess(prob, xu, xl, fu, n_feval, seed, 'ble');
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
