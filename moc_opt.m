function moc_opt(eim_process_name, prob)
% (1) intialize training data
% (2) start loop
%   (2-1) call EIMnext to generate next x  point
%   (2-2) combine new x into training data

workdir = pwd;
problem_folder = strcat(pwd,'\problems\EGproblems');
addpath(problem_folder);

prob = eval(prob);

hv_record = zeros(1, 30);
eim_function = str2func(eim_process_name);

for seed = 1:30
    fprintf(' seed: %d\n', seed);
    num_vari = prob.n_var;
    num_samples = 11 * num_vari - 1;
    num_con = prob.n_con;
    rng(seed, 'twister');

    % problem x range
    xu_bound = prob.xu;
    xl_bound = prob.xl;
    ref_point = prob.ref;
    
    %train data preparation
    train_x = lhsdesign(num_samples,num_vari,'criterion','maximin','iterations',1000);
    train_x = repmat(xl_bound, num_samples, 1) + repmat((xu_bound - xl_bound), num_samples, 1) .* train_x;
    % [train_y, train_c] = feval(fun_name, train_x);
    [train_y, train_c] = prob.evaluate(train_x);
    %figure(1)
    for iter=1:29
        
        index_c = sum(train_c <= 0, 2) == num_con;
        if sum(index_c) ~=0
            feasible_y = train_y(index_c, :);
            nd_index = Paretoset(feasible_y);
            nd_front = feasible_y(nd_index, :);
            % f1 = scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
            h = Hypervolume(nd_front,ref_point);
            %fprintf(' iteration: %d, hypervolume: %f\n',  iter,  h);
            
        end
        
        [newx, info] = eim_function(train_x, train_y, xu_bound, xl_bound, 20, 50, train_c);
        
        [newy, newc] =  prob.evaluate(newx);
        train_x = [train_x; newx];
        train_y = [train_y; newy];
        train_c = [train_c; newc];
        
        % f2 = scatter(newy(:,1), newy(:,2),'go', 'filled');
    end
    
    hv_record(seed) = h;
end

%record hv
filename=strcat(pwd, '\result_folder\', prob.name, '_cstill_hv.csv' );
csvwrite(filename, hv_record');
rmpath(problem_folder)
end