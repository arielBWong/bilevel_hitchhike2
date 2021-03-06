function [trainx, trainf, trainc, growflag] = ulego_sao_updateArchiveL(xu, popx, prob,  trainx, trainf, trainc, evalOne, distancecontrol)
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
    fprintf('Evolutional process converge to seen  xl\n');
    growflag = false;
    return;
end
newindex = ~repeatindex;
new_xl = popx(newindex, :);
xu_expand =  repmat(xu, n_new, 1);

% ---------------
% use distance control
if distancecontrol
    count = 0;
    for i =1:sum(newindex)
        newx = new_xl(i, :);
        tooclose = archive_check(newx, trainx, prob);
        if ~tooclose
            [new_fl, new_fc] = prob.evaluate_l(xu, newx);
            trainx = [trainx; newx];
            trainf = [trainf; new_fl];
            trainc = [trainc; new_fc];
            growflag = true;
            return;
        else
            count = count + 1;
        end
        
    end
    
    if count == sum(newindex)
        fprintf('Evolutional process unable to find new  matching xl\n');
        growflag = false;
        return;
    end
end


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


function tooclose = archive_check(newx, trainx, prob)
% ---check newx whether it is
tooclose = false;
eps_dist = sqrt(prob.n_lvar) * 0.01;

upper_bound = prob.xl_bu;
lower_bound = prob.xl_bl;

trainx_norm = (trainx - lower_bound) ./ (upper_bound - lower_bound);
newx_norm = (newx - lower_bound) ./ (upper_bound - lower_bound);

%---
mindistance = min(pdist2(newx_norm,trainx_norm));

if mindistance < eps_dist
    tooclose =  true;
end
end
