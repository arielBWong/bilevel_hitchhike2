function [pop,archive,evals,track]=Initialize_pop(prob,param,evals,objnum)
N = param.popsize;
pop.G=[];
infeas_ids=[];
G_pop=[];
id2=[];
funh = str2func(param.prob_name);
X_pop = repmat(prob.range(:,1)',N,1) + repmat((prob.range(:,2)'-prob.range(:,1)'),N,1).*lhsdesign(N,prob.nx);
[F_pop,~]=funh(objnum,X_pop);
[~,ids,~] = nd_sort(F_pop,[1:size(F_pop,1)]');

% Storing relevant information pop and archive
pop.X=X_pop(ids,:);
pop.F=F_pop(ids,:);
pop.Gen(1:param.popsize,:)=0;
pop.ischild=zeros(param.popsize,1);
archive.sols=[repmat(0,N,1) pop.X pop.F];
archive.parents=[repmat(0,N,1) pop.X pop.F pop.ischild];
evals=[evals;repmat(param.popsize,1,prob.nf)];

% Storing information in track
track.perc_corr=[100];

return

