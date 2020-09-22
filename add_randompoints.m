function newx = add_randompoints(lower_bound,upper_bound, n_var, trainx, num_eval)
% this method is called when unseen data is less than required number
% 
k = 1;
newx = [];
while k <= num_eval
    x_add = lhsdesign(1, n_var,'criterion','maximin','iterations',1000);
    x_add = lower_bound +(upper_bound - lower_bound) .* x_add;
    if ~ismember(x_add, trainx, 'row')
        newx = [newx; x_add];
        k = k + 1;
    else
        fprintf('added lower level to evaluation is seen one, continue initialize');
    end
end