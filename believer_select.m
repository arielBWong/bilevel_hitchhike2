%-----------------------------------------------------
%-- select one point --
function [newx, growflag] = believer_select(popx, trainx, prob, upper)
% select one point from believer's last population
% and avoid repeated point

repeatindex = ismember(popx, trainx, 'row');
n_new = sum(~repeatindex);
if n_new == 0
    fprintf('Evolutional process unable to find new  matching xl');
    growflag = false;
    [newx, growflag] = newpoint_randselection(prob, trainx, upper);
    return;
end
newindex = ~repeatindex;
new_xl = popx(newindex, :);

newx = new_xl(1,:);
growflag = true;

end