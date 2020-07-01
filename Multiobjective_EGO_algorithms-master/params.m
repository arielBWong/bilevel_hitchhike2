param.pop_size = 100;           %Population size
param.generations = 100;        %Number of generations

param.crossover_prob = 0.9;     %Crossover probability (close to 1 recommended) 
param.crossover_sbx_eta = 20;   %SBX crossover parameter (higher = steeper distribution)
param.mutation_prob = 0.1;      %Mutation probability (close to ~0-0.1 recommended) 
param.mutation_poly_eta = 30;   %Polynomial mutation parameter (higher = steeper distribution)

param.seed =2000;                %Run from same seed value will give same results
param.batch_mode = 0;           %Set 1 to suppress display during the run
param.max_fn_evals = 1e20;      %Max evaluations. Set 0 to disable. otherwise algo will run for min (popsize * generations, max_evals)
param.analysis_cache = 0;       %Implementation specific. No need to modify.

