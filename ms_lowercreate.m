%% recreate lower plot for dsm
% dsm has to be 2 variables

problems  ={'smd1(1, 1, 1)'};
prob = eval(problems{1});
method = {'llmatchhyb', 'llmatchble', 'llmatcheim', 'llmatchswitch'};


k = 2: prob.n_lvar;
k = (k-1)/2;
xu = [0, k]; 
xu = [0,0];
create_primefigure(xu, prob)
seed  = 7 ;
num = prob.n_lvar;

%------hyb
savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_',method{1});
filename  = strcat(savepath, '\xl_', num2str(seed), '.csv');
xl = csvread(filename);
% xl = xl(1:end-1, :);

% n = size(xl, 1);
% xu = repmat(xu, n, 1);
% 
% [fl, fc] = prob.evaluate_l(xu, xl);
% create_algfigure(xl, fl, fc, prob, 'hyb');

%------ble
% savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_',method{2});
% filename  = strcat(savepath, '\xl_', num2str(seed), '.csv');
% xl = csvread(filename);
% xl = xl(1:end-1, :);


n = size(xl, 1);
figure(2);

for ii = 21:59
    ii
    xu_prog  = repmat( xu, ii, 1);
    xl_prog = xl(1:ii, :);
    xl_next = xl(ii+1, :);
    [fl, fc] = prob.evaluate_l(xu_prog, xl_prog);
    [fl_next, fc_next] = prob.evaluate_l(xu, xl_next);
    create_algfigure(xl_prog, fl, fc, prob, 'hyb', xl_next, fl_next);
end





% %------eim
% savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_',method{3})
% filename  = strcat(savepath, '\xl_', num2str(seed), '.csv');
% xl = csvread(filename);
% % xl = xl(1:end-1, :);
% 
% n = size(xl, 1);
% 
% [fl, fc] = prob.evaluate_l(xu, xl);
% create_algfigure(xl, fl, fc, prob, 'eim');



function create_algfigure(xu, fu, fc, prob, t, x_next, f_next)
% close;
clf;

[krg_obj, krg_con, ~] = update_surrogate(xu,fu, fc, str2func('normalization_z' ));

m = 101;
xu1 = linspace(prob.xl_bl(1), prob.xl_bu(1), m);
xu2 = linspace(prob.xl_bl(2), prob.xl_bu(2), m);
[xu1, xu2] = meshgrid(xu1, xu2);
f = zeros(m, m);
for i = 1:m
    for j = 1:m
        [kk,  mse]= predictor([xu1(i, j), xu2(i, j)], krg_obj{1});
        f(i, j) = denormzscore( fu,kk); 
    end
end
min(f(:));

% x_next = xu(1,:);
[f_predict, ~]= predictor(x_next, krg_obj{1});
f_predict =  denormzscore( fu, f_predict); 



surf(xu1, xu2, f);hold on;
scatter3(xu(:, 1), xu(:, 2), fu, 80, 'rx')
% scatter3(x_next(1), x_next(1), f_next, 80, 'bo', 'filled');
scatter3(x_next(1), x_next(2), f_predict, 80, 'mo', 'filled');
xlabel('xl1','FontSize', 16);
ylabel('xl2', 'FontSize', 16);
zlabel('fl', 'FontSize', 16);
colormap jet
shading interp
t = strcat(t, '   ', prob.name, 'k', num2str(prob.n_lvar));
title(t);
drawnow; 
pause(1);
end


 function f = denormzscore(trainy, fnorm)
[train_y_norm, y_mean, y_std] = zscore(trainy);
f = fnorm * y_std + y_mean;
 end

 
 function create_primefigure(xu, prob)
figure(1);
m = 101;
xu1 = linspace(prob.xl_bl(1), prob.xl_bu(1), m);
xu2 = linspace(prob.xl_bl(2), prob.xl_bu(2), m);
[xu1, xu2] = meshgrid(xu1, xu2);
f = zeros(m, m);
for i = 1:m
    for j = 1:m
        xl = [xu1(i, j), xu2(i,j)];
        f(i, j) = prob.evaluate_l(xu, xl);
    end
end
surf(xu1, xu2, f);hold on;
colormap jet
shading interp
title('real landscape');
xlabel('xl1','FontSize', 16);
ylabel('xl2', 'FontSize', 16);
 end