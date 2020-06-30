%%
% This unittest is for 2 objectives 2 variables and 2 constraints
clearvars;
close all;

% (1) intialize training data
% (2) start loop
%   (2-1) call EIMnext to generate next x  point
%   (2-2) combine new x into training data

num_vari = 2;
num_samples = 11 * num_vari - 1;
num_con = 2;
fun_name = 'BNH';
rng(2, 'twister');
max_evaluation = 50;

% problem x range
xu_bound = [ 5, 3];
xl_bound = [0, 0];
ref_point = [140, 50];

%train data preparation
train_x = lhsdesign(num_samples,num_vari,'criterion','maximin','iterations',1000);
train_x = repmat(xl_bound, num_samples, 1) + repmat((xu_bound - xl_bound), num_samples, 1) .* train_x;
[train_y, train_c] = feval(fun_name, train_x);
figure(1)
for iter=1:29
    
    index_c = sum(train_c <= 0, 2) == num_con;
    if sum(index_c) ~=0
        feasible_y = train_y(index_c, :);
        nd_index = Paretoset(feasible_y);
        nd_front = feasible_y(nd_index, :);
        f1 = scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
        h = Hypervolume(nd_front,ref_point);
        fprintf(' iteration: %d, hypervolume: %f\n',  iter,  h);

    end
    
    [newx, info] = EIMnext(train_x, train_y, xu_bound, xl_bound, 20, 50,train_c);
     
    [newy, newc] =  feval(fun_name, newx);
    train_x = [train_x; newx];
    train_y = [train_y; newy];
    train_c = [train_c; newc];
    
    % f2 = scatter(newy(:,1), newy(:,2),'go', 'filled');    
end