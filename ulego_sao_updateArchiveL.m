function [trainx, trainf, trainc, growflag] = ulego_sao_updateArchiveL(xu, popx, prob,  trainx, trainf, trainc, evalOne)
% this method add new member(s) to achive picking
% from current population, new member has to be unseen data 
% w.r.t. archive
% input
% xu: matching xu
% popx: last population of evolutionary process
% prob:  problem instance
% trainx/trainf/train: archive
% evalOne: flag to signal whether only evaluate one or all
% ouput
% trainx/trainf/trainc: archive
%  growflag: whether archive has increased
%------------------------------------------------------------------------
repeatindex = ismember(popx, trainx, 'row');
n_new = sum(~repeatindex);
if n_new == 0
    fprintf('Evolutional process unable to find new  matching xl');
    growflag = false;
    return;
end
newindex = ~repeatindex;
new_xl = popx(newindex, :);
xu_expand =  repmat(xu, n_new, 1);

if evalOne  % acommodate one by one add
    [new_fl, new_fc] = prob.evaluate_l(xu_expand(1, :), new_xl(1,:));
    trainx = [trainx; new_xl(1, :)];
else
    [new_fl, new_fc] = prob.evaluate_l(xu_expand, new_xl);
    trainx = [trainx; new_xl];
end

trainf = [trainf; new_fl];
trainc = [trainc; new_fc];
growflag = true;


end
