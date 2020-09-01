
%%
y   = linspace(-1,2, 100);
f1 = y.^2;
f2 = (y-1).^2;

scatter(f1, f2); drawnow;
