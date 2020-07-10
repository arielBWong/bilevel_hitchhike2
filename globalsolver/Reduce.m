function [pop,archive] = Reduce(pop,param,archive,gen,prob)
% Reduction is based on information of solutions that have been completely evaluated 
ff=archive.sols(:,prob.nx+2:prob.nf+prob.nx+1);
id=find(sum(isnan(ff),2)==0);
fnd=ff(id,:);
[~,ids,~] = nd_sort(fnd,[1:size(fnd,1)]');
pop.X = archive.sols(id(ids(1:param.popsize)),2:prob.nx+1);
pop.F = archive.sols(id(ids(1:param.popsize)),2+prob.nx:prob.nf+prob.nx+1);
pop.Gen = pop.Gen(1:param.popsize);
pop.ischild = zeros(param.popsize,1);
archive.parents=[archive.parents ;[repmat(gen+1,param.popsize,1) pop.X pop.F pop.ischild]];
return
