function ulego_sao_archiveinsert(prob_str, seed, normhn)
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
% algorithm parameter
evaln = 1;

num_popu    = 20;   % 80 in total
num_genu    = 1200;
iter_frequ  = 20;

num_popl    = 20;   % 60 in total
num_genl    = 800;
iter_freql  = 20;


max_nl = 20000;
%--------
prob = eval(prob_str);
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
% -xu match its xl and evaluate fu
for i=1:num_popu
    [xl_single, nl, flag] = llmatch_sao_archiveinsert(xu(i, :), prob, num_popl, num_genl, iter_freql);
    xl = [xl; xl_single];
    llfeasi_flag = [llfeasi_flag, flag];
    n_feval = n_feval + nl;           %record lowerlevel nfeval
end

%--xu evaluation
[fu, fc] = prob.evaluate_u(xu, xl);
num_con = size(fc, 2);
scatter(fu(:, 1), fu(:, 2), 'ro', 'filled');


%--fu adjust
for i=1:num_popu
    fu = llfeasi_modify(fu, llfeasi_flag, i);
end

% -- create surrogate for first round evolution
[krg_obj, krg_con, ~] = update_surrogate(xu, fu, fc, normhn);

%--main population based optimization  routine
n = round(num_genu/iter_frequ);                             % how many  evolution  subroutine
initmatrix = xu;                                            % no normalization on x
for g=1:n
    disp(g);
    param.gen = iter_frequ;
    param.popsize = num_popu;
    
    ndnorm = get_ndfront(fu, fc, normhn);
    
    
    funh_obj = @(x)ulobj(x, krg_obj, ndnorm);
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
        fprintf('evolution converge and no new point is found\n');
        initmatrix = [];
        continue;
    end
    
    % --add new_xu to xu
    for i=1:evaln
        % match new_xl for new_xu
        % tic;
        [xl_single, nl, flag]  = llmatch_sao_archiveinsert(new_xu(i, :), prob, num_popl, num_genl, iter_freql);
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
    if n_feval > max_nl
        break;
    end
    %
    %-plot ----
    num_obj = size(fu, 2);
    ref_point = ones(1, num_obj) * 1.1;
    if ~isempty(fc) % constraint problems
        index_c = sum(fc <= 0, 2) == num_con;
        if sum(index_c) ~=0
            feasible_y = fu(index_c, :);
            nd_index = Paretoset(feasible_y);
            nd_front = feasible_y(nd_index, :);
            % f1 = scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
            % f2 = scatter(newfu(1), newfu(2), 'go', 'filled');
            num_nd = size(nd_front, 1);
            if num_nd > 1
                nd_front = (nd_front - min(nd_front))./(max(nd_front) - min(nd_front));
                h = Hypervolume(nd_front,ref_point);
                % fprintf(' iteration: %d, nd normalised hypervolume: %f\n',  i,  h);
            end
        end
    else  % unconstraint problems
        nd_index = Paretoset(fu);
        nd_front = fu(nd_index, :);
        clf('reset');
        % f1 = scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); hold on ;
        % f2 =scatter(newfu(:, 1), newfu(:,2), 'go', 'filled');drawnow;
        % f3 = scatter(expfu(1), expfu(2), 'bo', 'filled'); drawnow;
        num_nd = size(nd_front, 1);
        if num_nd >1
            nd_front = (nd_front - min(nd_front)) ./ (max(nd_front) - min(nd_front));
            h = Hypervolume(nd_front,ref_point);
            % fprintf(' iteration: %d, hypervolume: %f\n',  i,  h);
        end
    end
    %-plot ----
    
    % update krg  and initmatrix, continue to evolve
    [krg_obj, krg_con, ~] = update_surrogate(xu, fu, fc, normhn);
    % initmatrix = new_xu;
    % initmatrix = archive.pop_last.X;
    % [sf, sx, sc] = initmatrix_pick(xu,fu, fc);
    % initmatrix = sx(1:num_popu, :);
    % [initmatrix, ~] = unique(initmatrix,'rows','stable');
    initmatrix = [];
    
end
% plot nd fronts

% save data
nxu = n + num_popu;  % first generation and then every freq generations
nxl = n_feval;

disp(nxu);
disp(nxl);
method = 'sao_archiveinsert';
perfrecord_sao(xu, fu, fc, prob, seed, method, nxu, nxl);
end

function nd_frontnorm = get_ndfront(fu, fc, normhn)
% find nd front from fu, compatible with constraint problem
fu_norm = normhn(fu);
num_con = size(fc, 2);
if ~isempty(fc)     % constraint problems
    index_c = sum(fc <= 0, 2) == num_con;
    if sum(index_c) ~=0
        feasible_y = fu_norm(index_c, :);
        nd_index = Paretoset(feasible_y);
        nd_frontnorm = feasible_y(nd_index, :);
    else
        nd_frontnorm = [];
    end
else                % unconstraint problems
    nd_index = Paretoset(fu_norm);
    nd_frontnorm = fu_norm(nd_index, :);
end
end

function [sf, sx, sc] = initmatrix_pick(x, f, c)
[sf, sx, sc] = pop_sort(f, x, c);
end


function  f = ulobj(x, kriging_obj, nd_front)
num_obj = length(kriging_obj);   % krg cell array?
num_x = size(x, 1);
obj = zeros(num_x, num_obj);
for ii =1:num_obj
    [obj(:, ii), ~] = dace_predict(x, kriging_obj{ii});
end

% convert to hv contribution
ref_point = [1.1, 1.1];
basehv = Hypervolume(nd_front, ref_point);
f = zeros(num_x, 1);                                        % contribution of each x-f to existing nd front

for ii = 1:num_x                                              % population size
    extendf = [nd_front; obj(ii, :)];
    extendedhv = Hypervolume(extendf, ref_point);
    f(ii) = -(extendedhv - basehv);
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