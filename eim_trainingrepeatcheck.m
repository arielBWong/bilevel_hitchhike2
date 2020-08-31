function flag =  eim_trainingrepeatcheck(newx, trainx)
%  this function checks whether there is repeated point in training data
% flag: true - repeated 
%          false- no repeated
%-----------------------
repeat_index = ismember(newx, trainx, 'row');
if sum(repeat_index) > 0
    flag = true;
else
    flag = false;
end
end