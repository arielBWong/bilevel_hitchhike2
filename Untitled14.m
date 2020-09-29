%% this script is to change the results of smd
prob=smd1();
tic;
 
 xu_start = [0.0007, 0.0046];
 [newxu, newxl, n_up, n_low] = blsovler(prob, xu_start, 20, 20, 20, 10, 200)
 toc;
