%%
load('sig_gpr');
load('sig_dace');

min(sig_gpr)
min(sig_dace)

load('data')
x = data(:, 1);
y = data(:, 2);
normhn = str2func('normalization_z');
[krg_obj, krg_con, ~] = update_surrogate(x, y, [], normhn);


train_y_norm = normhn(y);
f_best = min(train_y_norm, [], 1);

prob = Shekel_curve();
rng(1, 'twister');

fighn = figure(1);



fithn = str2func('EIM_eval');
[new_xl, ~] = EIMnext_gpr(x, y, 5, -5, ...
         20, 20, [], fithn, normhn, fighn);