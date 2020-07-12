function [pop] = reduce_pop(pop, popsize)
% Reduction is based on information of solutions that have been completely evaluated 


pop.X = pop.X(1: popsize, :);
pop.F = pop.F(1: popsize, :);
pop.C = pop.C(1: popsize, :);

return
