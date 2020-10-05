%% ulego_corerecreate

%% reproduce upper level plots
% read xu and calculate best matching xl
% recreate plot
% read xu and fu
% recreate plot
%
problems = {'smd5()'};
algs = {'hyb', 'ble', 'eim'};
seed = 1;

prob = eval(problems{1});

%---------------------hyb
x_savefile = strcat(pwd, '\result_folder\', prob.name,'_', num2str(prob.n_lvar), '_',  algs{1});
xu_savefile = strcat(x_savefile, '\xu_', num2str(seed), '.csv' );
xu = csvread(xu_savefile);

xl_savefile = strcat(x_savefile, '\xl_', num2str(seed), '.csv' );
xl = csvread(xl_savefile);
[fu, fc] = prob.evaluate_u(xu,xl);
create_algfigure(xu, fu, fc, prob, 'hyb recreate');

%------------------ble
x_savefile = strcat(pwd, '\result_folder\', prob.name,'_', num2str(prob.n_lvar), '_',  algs{2});
xu_savefile = strcat(x_savefile, '\xu_', num2str(seed), '.csv' );
xu = csvread(xu_savefile);

xl_savefile = strcat(x_savefile, '\fu_', num2str(seed), '.csv' );
fu = csvread(xl_savefile);
create_algfigure(xu, fu, fc, prob, 'ble recreate');

%-----------------------eim
x_savefile = strcat(pwd, '\result_folder\', prob.name,'_', num2str(prob.n_lvar), '_',  algs{3});
xu_savefile = strcat(x_savefile, '\xu_', num2str(seed), '.csv' );
xu = csvread(xu_savefile);

xl_savefile = strcat(x_savefile, '\xl_', num2str(seed), '.csv' );
xl = csvread(xl_savefile);
[fu, fc] = prob.evaluate_u(xu,xl);
create_algfigure(xu, fu, fc, prob, 'eim recreate');

%--------------------- prime
create_primefigure(xu, prob);

function create_algfigure(xu, fu, fc, prob, t)
figure
[krg_obj, krg_con, ~] = update_surrogate(xu,fu, fc, str2func('normalization_z' ));

xu1 = linspace(prob.xu_bl(1), prob.xu_bu(1), 100);
xu2 = linspace(prob.xu_bl(2), prob.xu_bu(2), 100);
[xu1, xu2] = meshgrid(xu1, xu2);
f = zeros(100, 100);
for i = 1:100
    for j = 1:100
        f (i, j)= predictor([xu1(i, j), xu2(i, j)], krg_obj{1});
        f(i, j) = denormzscore( fu, f(i, j)); 
    end
end
min(f(:))

surf(xu1, xu2, f);hold on;
scatter3(xu(:, 1), xu(:, 2), fu, 80, 'rx')
xlabel('xu1','FontSize', 16);
ylabel('xu2', 'FontSize', 16);
zlabel('fu', 'FontSize', 16);
colormap jet
shading interp
title(t);
end

function create_primefigure(xu, prob)
xl = [];
for ii = 1: size(xu, 1)
    xl_prime = matching_prime(xu(ii, :), prob);
    xl = [xl; xl_prime];
end

[fu, fc] = prob.evaluate_u(xu, xl);
[krg_obj, krg_con, ~] = update_surrogate(xu,fu, fc, str2func('normalization_z' ));

xu1 = linspace(prob.xu_bl(1), prob.xu_bu(1), 100);
xu2 = linspace(prob.xu_bl(2), prob.xu_bu(2), 100);
[xu1, xu2] = meshgrid(xu1, xu2);
f = zeros(100, 100);
for i = 1:100
    for j = 1:100
        f (i, j)= predictor([xu1(i, j), xu2(i, j)], krg_obj{1});
        f(i, j) = denormzscore( fu, f(i, j)); 
    end
end
min(f(:))
figure;
surf(xu1, xu2, f); hold on;
scatter3(xu(:, 1), xu(:, 2), fu, 80, 'yx')
xlabel('xu1','FontSize', 16);
ylabel('xu2', 'FontSize', 16);
zlabel('fu', 'FontSize', 16);
colormap jet
shading interp
title('prime');
end

function xl_prime =  matching_prime(xu, prob)
xl1 = ones(1, prob.q);
xu2 = xu(1, prob.p + 1:end);
xl2 = sqrt(abs(xu2));
xl_prime = [xl1, xl2];
end


 function f = denormzscore(trainy, fnorm)
[train_y_norm, y_mean, y_std] = zscore(trainy, 1, 1);
f = fnorm * y_std + y_mean;
end