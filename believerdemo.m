%% 
 % prob = Forrestor();
prob = rastrigin(1,1);
rng(1, 'twister');
% [match_xl, n_fev, flag] = llmatch_sao_archiveinsert([0, 0], prob, 20, 5, 20, 20, 1, 'normalization_z');% 
[match_xl, n_fev, flag] = llmatch([0, 0], prob, 20, 20 , 'EIMnext', 5, 25, 'EIM_evaldace', 1,'normalization_z' );
% [match_xl, n_fev, flag] = llmatchgpr([0, 0], prob, 20, 20, 'EIMnext_gpr', 5, 25, 'EIM_eval', 1,'normalization_z' );

% llmatch_believeradapt([0, 0], prob, 20, 20, 'EIMnext', 5, 30,  'EIM_evaldace',  1,'normalization_z');

disp(match_xl);
f = prob.evaluate_l([], match_xl);
