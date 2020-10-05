
%-----------------------------------------------------
%-- repeat point --
function [newx, flag] = newpoint_randselection(prob, trainx, upper)
% function finds a new point that is unseen from trainx archive
% 50 tries, if all fails, return a random one


if upper  % Upper level point
    lower_bound = prob.xu_bl;
    upper_bound = prob.xu_bu;
    n_var = prob.n_uvar;
else
    lower_bound = prob.xl_bl;
    upper_bound = prob.xl_bu;
    n_var = prob.n_lvar;
end

n = 50;  % Control max number of search attempts for unseen data
m = 1;
while m<n
    x_add = lhsdesign(1, n_var, 'criterion', 'maximin', 'iterations', 1000);
    x_add = lower_bound +(upper_bound - lower_bound) .* x_add;
    if ~ismember(x_add, trainx, 'row')
        new_x = x_add;
        flag = true;
        break;
    else
        fprintf('new point is a seen one, re-search');
        m = m + 1;
        if m == n
            fpritnf('no new point in 50 tries, returns a random one');
            newx = lhsdesign(1, n_var, 'criterion', 'maximin', 'iterations', 1000);
            newx = lower_bound +(upper_bound - lower_bound) .* newx;
            flag = false;
        end
    end
end

end

