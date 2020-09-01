
%%
y   = linspace(-1, 2, 100);
f1 = (y-1).^2 +  y.^2;
f2 = (y-1).^2 * 2;
nd = [f1', f2'];
 nd_index = Paretoset(nd);
pf= nd(nd_index, :);
    

scatter(pf(:, 1), pf(:,2)); drawnow;
