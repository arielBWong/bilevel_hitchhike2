load('data')
x = data(:, 1);
y = data(:, 2);

normhn = str2func( 'normalization_z');
[krg_obj, krg_con, info] = update_surrogatedace(x, y, [],normhn);




xt = linspace(-5.12, 5.12, 100);
yt = dace_predict(xt', krg_obj{1});
yt = denormzscore(y, yt);

plot(xt', yt, 'r--'); hold on;
plot(x, y, 'bo');