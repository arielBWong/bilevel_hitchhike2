%% problem test
prob = tp3();
xl =linspace(0, 1, 101);

f = prob.evaluate_l([], xl');
plot(xl, f); 

xu = [1,1 ];
% fm = prob.evaluate_l([], [-9]);
[match_xl, n_fev, flag] = llmatch(xu, prob, 20, 20, 'EIMnext', 100, 'EIM_eval', 1)
% [match_xl, n_fev, flag] =  llmatch_sao_archiveinsert(xu, prob,20, 100, 20, 1)

a = 0;