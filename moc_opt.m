function moc_opt(eim_process_name, prob)
% this test is a unit test function for compare EIMnext, EIMnext_znorm
% and paper demo
% (1) intialize training data
% (2) start loop
%   (2-1) call EIMnext to generate next x  point
%   (2-2) combine new x into training data




prob = eval(prob);
fprintf('%s', eim_process_name);
fprintf('%s', prob.name);

% number of iteration setting
if prob.n_obj>2
    maxeval = 200;
else
    maxeval = 100;
end

if nargin>2
    maxeval = 50;
end

hv_record = zeros(1, 10);
eim_function = str2func(eim_process_name);

for seed = 1:11
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
    % figure(1)
    
    maxiter = maxeval - num_samples;
    for iter=1:maxiter
        
        [newx, info] = eim_function(train_x, train_y, xu_bound, xl_bound, 50, 200, train_c);
        
        [newy, newc] =  prob.evaluate(newx);
        
        train_x = [train_x; newx];
        train_y = [train_y; newy];
        train_c = [train_c; newc];
        
        % for plot
        if ~isempty(train_c) % constraint problems
            index_c = sum(train_c <= 0, 2) == num_con;
            if sum(index_c) ~=0
                feasible_y = train_y(index_c, :);
                nd_index = Paretoset(feasible_y);
                nd_front = feasible_y(nd_index, :);
                % f1 = scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); drawnow;
                h = Hypervolume(nd_front,ref_point);
                fprintf(' iteration: %d, hypervolume: %f\n',  iter,  h);
            end
        else  % unconstraint problems
            nd_index = Paretoset(train_y);
            nd_front = train_y(nd_index, :);
            %             clf('reset');
            %             scatter(nd_front(:,1), nd_front(:,2),'ro', 'filled'); hold on;
            %             scatter(pareto_front(:,1), pareto_front(:,2),'bo', 'filled');
            %             drawnow;
            h = Hypervolume(nd_front,ref_point);
            fprintf(' iteration: %d, hypervolume: %f\n',  iter,  h);
        end
    end
    nd_index = Paretoset(train_y);
    nd_front = train_y(nd_index, :);
    h = Hypervolume(nd_front,ref_point);
    hv_record(seed) = h;
    
    fprintf('hv redord %.6f', h);
    
    
    %     filename2=strcat(pwd, '\result_folder\',eim_process_name,'_', prob.name, '_',num2str(seed), '_trainy.csv' );
    %     filename3=strcat(pwd, '\result_folder\',eim_process_name,'_', prob.name, '_',num2str(seed), '_trainc.csv' );
    %     csvwrite(filename2, train_y); % for plot
    %     csvwrite(filename3, train_c); % for plot
    %
    
end
%record hv
filename1=strcat(pwd, '\result_folder\',eim_process_name,'_', prob.name, '_hv.csv' );
csvwrite(filename1, hv_record'); % make sure column

fprintf('hv mean: %.6f', mean(hv_record));
fprintf('hv median: %.6f', median(hv_record));

end