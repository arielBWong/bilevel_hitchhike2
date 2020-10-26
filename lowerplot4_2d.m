function[] = lowerplot4_2d(prob, trainx,  trainy, trainc)
 
if size(trainy, 2) > 1 || size(trainx, 2) ~= 2
    error('not compatible for mo');
end


normhn = str2func('normalization_y');
[krg_obj, krg_con, info] = update_surrogate(trainx, trainy, trainc,normhn );

% figure(1);
% 
lb1 = prob.xl_bl(1);
ub1 = prob.xl_bu(1);

lb2 = prob.xl_bl(2);
ub2 = prob.xl_bu(2);

num_points = 101;

x1 = linspace(lb1, ub1, num_points);
x2 = linspace(lb2, ub2, num_points);
[x1, x2] = meshgrid(x1, x2);
f = zeros(num_points, num_points);







% 
% 
% for i = 1:num_points
%     for j = 1:num_points
%         [f(i, j), ~] = dace_predict([x1(i, j), x2(i, j)], krg_obj{1});
%         f(i, j) = denormzscore( trainy,  f(i, j));
%     end
% end
% min(f(:))
% surf(x1, x2, f); hold on;
% xlabel('x1', 'FontSize', 16);
% ylabel('x2', 'FontSize', 16);
% zlabel('f',  'FontSize', 16);
% colormap jet
% shading interp
% 

num_points = 101;
figure(2);
f = zeros(num_points, num_points);


for i = 1:num_points
    for j = 1:num_points
        f(i, j) = prob.evaluate_l([], [x1(i, j), x2(i, j)]);
    end
end
% 
% surfc(x1, x2, f); hold on;
% 
% xlabel('x1', 'FontSize', 16);
% ylabel('x2', 'FontSize', 16);
% zlabel('f',  'FontSize', 16);
% colormap jet
% shading interp
% title(prob.name,'FontSize', 18 );

contour(x1, x2, f); hold on;
colormap jet
title(prob.name,'FontSize', 18 );



end