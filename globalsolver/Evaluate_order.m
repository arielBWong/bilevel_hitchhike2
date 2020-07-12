function [pop, archive]= evaluate_order(pop, archive, funh_obj, funh_con, cx, gen, param)
% This function evaluate given x population, add to pop (N to 2N)
% and sort pop with nd and feasibility consideration
% input 
%           pop : population of previous generation
%           archive : record of all evolutionary process
%           funh_obj : function handle of objective function
%           funh_con: function handle of constraints
%           cx: child population generated from last step
%           gen: current generation
%           param: evolution parameter (gen popsize)
% output
%         pop: extended and sorted population
%         archive: extended archive
%-------------------------
child.X=cx;
% This is NSGA-II
child.F = funh_obj(child.X);
child.C = funh_con(child.X);
numcon = size(child.C, 2);

archive.sols=[archive.sols;[repmat(gen,param.popsize,1), child.X, child.F, child.C]];

% Appending X F C to pop
pop.X=[pop.X;child.X];
pop.F=[pop.F;child.F];
pop.C = [pop.C; pop.C];

[pop.F, pop.X, pop.C] = pop_sort(pop.F, pop.X, pop.C);
end

function [sf, sx, sc] = pop_sort(f, x, c)
% auxiliary function
% this function sort evoluation popluation w.r.t. 
% number of objectives
% constraints
% input
%               f                                                       %  objective population 
%               x                                                      %  design variable population 
%               c                                                      %  constraints population 
% output
%               sf                                                     % sorted objective population 
%               sx                                                     % sorted design variable population 
%               sc                                                     % sorted constraints population 
%-----------------------------------------------------

numcon = size(c, 2);

if numcon == 0                                             % unconstraint problem
    sc = [];
    [~, ids, ~] = nd_sort(f, [1:size(f, 1)]);
    sf = f(ids,:);
    sx = c(ids, :); 
end


if numcon>0
% Ordering considers constraints
c(c<=0) = 0;
fy_ind = sum(c, 2) == numcon;                 % feasibility index 
cv_ind =~fy_ind;                                            % constraint violation index


% seperate feasible and infeasible
% sort two subset seperately
fy_F = f(fy_ind, :);  cv_F = f(cv_ind, :);         % operation should be valid when no feasible solutions
fy_C = c(fy_ind, :); cv_C = c(cv_ind, :);
fy_X = x(fy_ind, :); cv_X = x(cv_ind, :);

% sort feasible
[~, ids, ~] = nd_sort(fy_F, [1: size(fy_F, 1)]);  fy_F = fy_F(ids, :); fy_C = fy_C(ids, :); fy_X = fy_X(ids, :);

% sort infeasible
[~, idc] = sort(cv_C); cv_F = cv_F(idc, :); cv_C = cv_C(idc, :); cv_X = cv_X(idc, :);

% replace unsorted each fields of pop
sf = [fy_F; cv_F]; sc= [fy_C; cv_C]; sx = [fy_X; cv_X];
end
end




%----------------------------------------------------------------------------------------------------
% if(param.strategy==1)
%     flag=0;tmp_evals=zeros(1,prob.nf);
%     % Appending
%     pop.X=[pop.X;child.X];
%     pop.F=[pop.F;NaN*ones(param.popsize,prob.nf)];
%     pop.ischild=[pop.ischild; ones(param.popsize,1)];
%     pop.Gen=[pop.Gen; gen*ones(param.popsize,1)];
%     child.F=NaN*ones(param.popsize,prob.nf);
%
%     % Updating the archive
%     archive.sols=[archive.sols;[repmat(gen,param.popsize,1) child.X child.F]];
%
%     while flag==0
%
%         % Computing mu and sigma of the offsprings
%         [mu,sigma]=Compute_mu_sigma(pop,prob,param,archive);
%
%         % Selecting which solution to be evaluated in which objective
%         [who,which_obj]=Select(mu,sigma,prob);
%         tmp=unique([who' which_obj'],'rows');
%         who=tmp(:,1);which_obj=tmp(:,2);
%
%         % Stepping solutions in who list which_obj order
%         [F_tmp,~] = funh(objnum,pop.X(who,:));
%         for i=1:length(who)
%             % Updating pop.F
%             pop.F(who(i),which_obj(i))=F_tmp(i,which_obj(i));
%             tmp_evals(which_obj(i))=tmp_evals(which_obj(i))+1;
%             % Updating Archive
%             [~,ipx]= ismember(pop.X(who(i),:),archive.sols(:,2:prob.nx+1),'rows');
%             archive.sols(ipx,1+prob.nx+which_obj(i)) = pop.F(who(i),which_obj(i));
%         end
%
%         % Number of fully evaluated solutions
%         if length(find(sum(isnan(pop.F),2)==0))>param.popsize
%             flag=1;
%         end
%     end
% end
%------------------------------------------------------------------

% % Updating evaluations
% evals=[evals;tmp_evals];
%
% % Among the offsprings did the scheme choose a ND and fully evaluated it
% [child.F,~]=funh(objnum,child.X);
% tmp=[pop.F(1:param.popsize,:) ;child.F];
% [fronts,~,~] = nd_sort(tmp,[1:size(tmp,1)]');
% % These are ND solutions if the childpop was evaluated
% id1=fronts(1).f;
%
% % The scheme fully evaluated
% tmp=find(sum(isnan(pop.F),2)==0);
% id_evaluated=tmp(find(tmp>param.popsize));
% track.perc_corr=[track.perc_corr 100*length(intersect(id1,id_evaluated))/length(id_evaluated)];
% return
