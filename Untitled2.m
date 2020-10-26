
%%
prob = newBranin2();
train_xl = lhsdesign(20, 2,'criterion','maximin','iterations',1000);
train_xl = repmat(prob.xl_bl, 20, 1) ...
    + repmat((prob.xl_bu - prob.xl_bl),20, 1) .* train_xl;

a = 0;

train_fl = prob.evaluate_l([], train_xl);

 lowerplot4_2d(prob, train_xl, train_fl, []);

% 
% figure(1);
% 
% lb1 = prob.xl_bl(1);
% ub1 = prob.xl_bu(1);
% 
% lb2 = prob.xl_bl(2);
% ub2 = prob.xl_bu(2);
% 
% num_points = 101;
% 
% x1 = linspace(lb1, ub1, num_points);
% x2 = linspace(lb2, ub2, num_points);
% 
% [x1, x2] = meshgrid(x1, x2);
% f = zeros(num_points, num_points);
% 
% 
% for i = 1:num_points
%     for j = 1:num_points
%         f(i, j) = prob.evaluate_l([], [x1(i, j), x2(i, j)]);
%     end
% end
% 
% contour(x1, x2, f);
% colorbar;
% 
% a = 0