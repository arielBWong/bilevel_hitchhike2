%% plot median mse in every 10 new training data
%
% for a seed per problem
% (1) load xl, generate xu
% (2) iterate from  index 20 every 10 more until 60,
%       read xl, match xu, train dace, predict best solution
% (3) save mse to another matrix
% (4)

clearvars;
close all;

seedmax = 11;
          
 problems ={'smd10(1,1,1)'};     

                       
methods = {'llmatcheim',  'llmatcharchive',  'llmatchpop'};  % 'llmatchpop',
leg = {'EIM', 'BLE'};
np= length(problems);
nm = length(methods);


% lower level has fixed number of true evaluation before local search
% 60 (20, 30, 40, 50, 60)
% last one saved for local search final results
collectmatrix = zeros(, np*nm*2);

for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    nv = prob.n_lvar;
    
    xu = [1, 1];
    

end

 