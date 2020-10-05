function llmatch_globalmin_cmp(prob, match_method, seed)
% this function takes xu and problem as input
% go through different lower level match method
% return global match accuracy
% save to certain files
% input

prob = eval(prob);

% k = 2: prob.n_lvar;
% k = (k-1)/2;
% xu = [0, k];

% smd10
xu = [0, 0];

rng(seed, 'twister');

if ~contains(match_method, '_')  
    match_method =  str2func( match_method);
    match_method(xu, prob, 20, 20, 'EIMnext', 40, 'EIM_eval', seed);

elseif contains(match_method, 'hyb')
    match_method =  str2func( match_method);
    match_method(xu, prob, 20, 20, 'EIMnext', 20, 'EIM_eval', seed);
else    
    num_pop = 20;
    if contains(match_method, 'archive')
         num_gen =800;
         freq =20;
   
    else
        num_gen =40;
         freq =20;
    end
    
    match_method =  str2func( match_method);
    match_method (xu, prob,num_pop, num_gen, freq, seed);
end

end



