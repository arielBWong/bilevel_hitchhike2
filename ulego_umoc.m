function ulego_umoc(prob, seed, str_nextxhn, fitnesshandle, normhn, llmatch_nextx)
% method of main optimization process of upper level ego
% adapt to upper level problems of "multiple objectives"
% usage:
%     input
%               prob: problem instance
%               seed: random process seed
%               str_nextxhn: string, upper level method to propose next x
%               fitnesshandle: string, the function that upper level next x method uses for its evolution process
%               normhn: string, normalization function used in EIMnext_znorm
%               llmatch_nextx: string, lower level match next x propose method
%             
%
%
%     output
%               csv files saved in result folder
%               performance statistics include 3*3 matrix
%                                                       [  ul  accuracy,                        ll accuracy;
%                                                          upper number of feval,    lower number of feval;
%                                                          upper feasibility,                lower feasibility]
%
%
%--------------------------------------------------------------------------

rng(seed, 'twister');
fprintf(prob);
prob = eval(prob);

% save some runs
% savepath = strcat('C:\Users\z3276872\matlab_scripts\bilevel_hitchhike2\result_folder\', prob.name, '_EIM_eval');
% file = strcat(savepath, '\fu_', num2str(seed),'.csv');
% if exist(file,'file') == 2  % ignore existing runs 
%     return;
% end

% algo parameter
numiter_l               = 40; %  100 intotal
initsize_l              = 20;
numiter_u               = 60;
inisize_u               = 20;
num_pop                 = 20;
num_gen                 = 20;
max_nl                  = 20000;

% parallel compatible setting

nextxhn = str2func(str_nextxhn);
fithn = str2func(fitnesshandle);
normhn = str2func(normhn);

n_feval = 0;
%--upper problem variable
u_nvar = prob.n_uvar;
upper_bound = prob.xu_bu;
lower_bound = prob.xu_bl;
% inisize_u = 11 * u_nvar - 1;

%-xu initialization
xu = lhsdesign(inisize_u,u_nvar,'criterion','maximin','iterations',1000);
xu = repmat(lower_bound, inisize_u, 1) + repmat((upper_bound - lower_bound), inisize_u, 1) .* xu;

xl = [];
llfeasi_flag = [];
% -xu match its xl and evaluate fu
for i=1:inisize_u
    [xl_single, n, flag] = llmatch(xu(i, :), prob,num_pop, num_gen,llmatch_nextx, numiter_l, fitnesshandle);
    xl = [xl; xl_single];
    llfeasi_flag = [llfeasi_flag, flag];
    n_feval = n_feval + n; %record lowerlevel nfeval
end
%--xu evaluation
[fu, fc] = prob.evaluate_u(xu, xl);
num_con = size(fc, 2);
% scatter(fu(:, 1), fu(:, 2), 'ro', 'filled');

%--fu adjust
for i=1:inisize_u
    fu = llfeasi_modify(fu, llfeasi_flag, i);
end

%-main ulego routine
for i = 1:numiter_u
    %--search next xu
    [newxu, info] = nextxhn(xu, fu, upper_bound, lower_bound,num_pop, num_gen, fc, fithn, normhn);
    
    %---test on recreating expected fu from kriging
    % expfu = expectedfu_fromkrg(newxu, info);
    
    %--get its xl
    [newxl, n, flag] = llmatch(newxu, prob,num_pop, num_gen, llmatch_nextx, numiter_l, fitnesshandle);
    % fprintf('xl matching feasibility is %d \n', flag);
    n_feval = n_feval + n;
    %--evaluate xu
    [newfu, newfc] = prob.evaluate_u(newxu, newxl);
    %--assemble xu fu fc
    xu = [xu; newxu];
    xl = [xl; newxl];
    fu = [fu; newfu];
    fc = [fc; newfc];
    llfeasi_flag = [llfeasi_flag, flag];
    %--adjust fu by lower feasibility
    fu = llfeasi_modify(fu, llfeasi_flag, inisize_u+i);                    % upper mo compatible
    
    if n_feval > max_nl
        fprintf(num2str(n_feval));
        break;
    end
    
    %--plot ----
    num_obj = size(fu, 2);
    ref_point = ones(1, num_obj) * 1.1;
    if ~isempty(fc) % constraint problems
        index_c = sum(fc <= 0, 2) == num_con;
        if sum(index_c) ~=0
            feasible_y = fu(index_c, :);
            nd_index = Paretoset(feasible_y);
            nd_front = feasible_y(nd_index, :);
             f1 = scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
             f2 = scatter(newfu(1), newfu(2), 'go', 'filled');
           
        end
    else  % unconstraint problems
        nd_index = Paretoset(fu);
        nd_front = fu(nd_index, :);
        clf('reset');
        f1 = scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); hold on ;
        f2 = scatter(newfu(1), newfu(2), 'go', 'filled');drawnow;
        % f3 = scatter(expfu(1), expfu(2), 'bo', 'filled'); drawnow;
        num_nd = size(nd_front, 1);
        if num_nd >1
            nd_front = (nd_front - min(nd_front)) ./ (max(nd_front) - min(nd_front));
            h = Hypervolume(nd_front,ref_point);
           %  fprintf(' iteration: %d, hypervolume: %f\n',  i,  h);
        end
    end
    %--plot ----
end
nxu = size(xu, 1);
nxl = n_feval;
perfrecord_umoc(xu, xl, fu, fc, prob, seed, fitnesshandle,nxu, nxl);

end

% auxiliary function for investigation
function [expfu] = expectedfu_fromkrg(xu, info)
% xu: returned newxu from eim
% info: struct returned from eim
% expfu: fu from  krigin
xu_norm = (xu - info.train_xmean) ./ info.train_xstd;
num_obj = length(info.krg);
yu_norm = zeros(1, num_obj);

for ii = 1:num_obj
    [yu_norm(ii), ~ ]= dace_predict(xu_norm, info.krg{ii});
end

% denormalize
yu_znorm = yu_norm .* (info.train_ynormmax - info.train_ynormmin) + info.train_ynormmin;
yu = yu_znorm .* info.train_ystd + info.train_ymean;

expfu = yu;
end
