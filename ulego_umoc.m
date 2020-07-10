function ulego_umoc(prob, seed, eim)
% method of main optimization process of upper level ego
% adapt to upper level problems of "multiple objectives"
% usage:
%     input
%               prob                          : problem instance
%               seed                          : random process seed
%               eim                            : string, the name of eim function 
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
hy_gen                   = 50;

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
scatter(fu(:, 1), fu(:, 2), 'ro', 'filled');


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
    fprintf('xl is %d \n', flag);
    n_feval = n_feval + n;
    %--evaluate xu
    [newfu, newfc] = prob.evaluate_u(newxu, newxl);
    %--assemble xu fu fc
    xu = [xu; newxu];
    fu = [fu; newfu];
    fc = [fc; newfc];
    llfeasi_flag = [llfeasi_flag, flag];
    %--adjust fu by lower feasibility
    fu = llfeasi_modify(fu, llfeasi_flag, inisize_u+i);                    % upper mo compatible
    
    %-plot ----
    num_obj = size(fu, 2);
    ref_point = ones(1, num_obj) * 1.1;
    if ~isempty(fc) % constraint problems 
        index_c = sum(fc <= 0, 2) == num_con;
        if sum(index_c) ~=0
            feasible_y = fu(index_c, :);
            nd_index = Paretoset(feasible_y);
            nd_front = feasible_y(nd_index, :);
            f1 = scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
            num_nd = size(nd_front, 1);
            if num_nd > 1
                nd_front = (nd_front - min(nd_front))./(max(nd_front) - min(nd_front));
                h = Hypervolume(nd_front,ref_point);
                fprintf(' iteration: %d, nd normalised hypervolume: %f\n',  i,  h);
            end
        end
    else  % unconstraint problems
        nd_index = Paretoset(fu);
        nd_front = fu(nd_index, :);
        clf('reset');
        scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
        num_nd = size(nd_front, 1);
        if num_nd >1
            nd_front = (nd_front - min(nd_front)) ./ (max(nd_front) - min(nd_front));
            h = Hypervolume(nd_front,ref_point);
            fprintf(' iteration: %d, hypervolume: %f\n',  i,  h);
        end
    end
    %-plot ----
end
toc
end