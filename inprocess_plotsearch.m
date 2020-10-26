function inprocess_plotsearch(fighn, prob, cons_hn, newx, trainx)
clf(fighn);


lb1 = prob.xl_bl(1);
ub1 = prob.xl_bu(1);

lb2 = prob.xl_bl(2);
ub2 = prob.xl_bu(2);

num_points = 101;

x1 = linspace(lb1, ub1, num_points);
x2 = linspace(lb2, ub2, num_points);
[x1, x2] = meshgrid(x1, x2);
f = zeros(num_points, num_points);


for i = 1:num_points
    for j = 1:num_points
        f(i, j) = prob.evaluate_l([], [x1(i, j), x2(i, j)]);
    end
end

contour(x1, x2, f); hold on;
colormap jet
fimplicit(cons_hn, [lb1, ub1, lb2, ub2], '--g','LineWidth',2 );
scatter(newx(1), newx(2), 'ro', 'filled');
scatter(trainx(:, 1), trainx(:, 2), 'bo');
scatter(prob.xprime(1), prob.xprime(2), 'ko', 'filled');
title(prob.name,'FontSize', 18 );
pause(1);



end