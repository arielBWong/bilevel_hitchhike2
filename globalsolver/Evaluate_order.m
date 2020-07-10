function [pop,archive,evals,track]= Evaluate_order(pop,prob,param,archive,evals,objnum,cx,gen,track)
child.X=cx;
funh = str2func(param.prob_name);

% This is NSGA-II
if(param.strategy==2)
    [child.F,~]=funh(objnum,child.X);
    
    % Updating the archive
    archive.sols=[archive.sols;[repmat(gen,param.popsize,1) child.X child.F]];
    
    % Appending
    pop.X=[pop.X;child.X];
    pop.F=[pop.F;child.F];
    pop.ischild=[pop.ischild; ones(param.popsize,1)];
    pop.Gen=[pop.Gen; gen*ones(param.popsize,1)];
    
    % Ordering 2N
    [~,ids,~] = nd_sort(pop.F,[1:size(pop.F,1)]');
    
    
    % Reassign pop based on order
    pop.X=pop.X(ids,:);
    pop.F=pop.F(ids,:);
    pop.Gen=pop.Gen(ids);
    pop.ischild=pop.ischild(ids);
    tmp_evals=[repmat(param.popsize,1,prob.nf)];
    evals=[evals;tmp_evals];
end

if(param.strategy==1)
    flag=0;tmp_evals=zeros(1,prob.nf);
    % Appending
    pop.X=[pop.X;child.X];
    pop.F=[pop.F;NaN*ones(param.popsize,prob.nf)];
    pop.ischild=[pop.ischild; ones(param.popsize,1)];
    pop.Gen=[pop.Gen; gen*ones(param.popsize,1)];
    child.F=NaN*ones(param.popsize,prob.nf);
    
    % Updating the archive
    archive.sols=[archive.sols;[repmat(gen,param.popsize,1) child.X child.F]];
    
    while flag==0
        
        % Computing mu and sigma of the offsprings
        [mu,sigma]=Compute_mu_sigma(pop,prob,param,archive);
        
        % Selecting which solution to be evaluated in which objective
        [who,which_obj]=Select(mu,sigma,prob);
        tmp=unique([who' which_obj'],'rows');
        who=tmp(:,1);which_obj=tmp(:,2);
        
        % Stepping solutions in who list which_obj order
        [F_tmp,~] = funh(objnum,pop.X(who,:));
        for i=1:length(who)
            % Updating pop.F
            pop.F(who(i),which_obj(i))=F_tmp(i,which_obj(i));
            tmp_evals(which_obj(i))=tmp_evals(which_obj(i))+1;
            % Updating Archive
            [~,ipx]= ismember(pop.X(who(i),:),archive.sols(:,2:prob.nx+1),'rows');
            archive.sols(ipx,1+prob.nx+which_obj(i)) = pop.F(who(i),which_obj(i));
        end
        
        % Number of fully evaluated solutions
        if length(find(sum(isnan(pop.F),2)==0))>param.popsize
            flag=1;
        end
    end
end

% Updating evaluations
evals=[evals;tmp_evals];

% Among the offsprings did the scheme choose a ND and fully evaluated it
[child.F,~]=funh(objnum,child.X);
tmp=[pop.F(1:param.popsize,:) ;child.F];
[fronts,~,~] = nd_sort(tmp,[1:size(tmp,1)]');
% These are ND solutions if the childpop was evaluated
id1=fronts(1).f;

% The scheme fully evaluated
tmp=find(sum(isnan(pop.F),2)==0);
id_evaluated=tmp(find(tmp>param.popsize));
track.perc_corr=[track.perc_corr 100*length(intersect(id1,id_evaluated))/length(id_evaluated)];
return
