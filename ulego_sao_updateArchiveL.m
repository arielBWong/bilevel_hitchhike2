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
    %------if new_xl is less than population size, add more random unseen  ones --
    num_new = sum(newindex);
    num_pop = size(popx, 1);
    k = num_new + 1;
    lower_bound = prob.xl_bl;
    upper_bound = prob.xl_bu;
    n_var = prob.n_lvar;
    while k <= num_pop
        xl_add = lhsdesign(1, n_var,'criterion','maximin','iterations',1000);
        xl_add = lower_bound +(upper_bound - lower_bound) .* xl_add;
        if ~ismember(xl_add, trainx, 'row')
            new_xl = [new_xl; xl_add];
            k = k + 1;
        else
            fprintf('added lower level to evaluation is seen one');
        end
    end
    %-------------------------------------------------------   
    % --add new_xl to xu
    num_new = size(new_xl, 1);
    xu_expand =  repmat(xu, num_new, 1);
    if num_new ~= num_pop
        error('something wrong with adding random xl to evaluation, should be the same size as population');
    end
    
    [new_fl, new_fc] = prob.evaluate_l(xu_expand, new_xl);
    trainx = [trainx; new_xl];
end

trainf = [trainf; new_fl];
trainc = [trainc; new_fc];
growflag = true;


end
