% -----------------------------------------------------------------------------------------
% 1. The multiobjective EGO algorithm using EIM(expected improvement
%    matrix)-based criteria, which is significant cheaper-to-evaluate than the
%    state-of-the-art multiobjective EI criteria. For detailed description
%    about the EIM criteria, please refer to [1].
% 2. The dace toolbox [2] is used for building the Kriging models in the
%    implementations.
% 3. The non-dominated sorting method by Yi Cao [3] is used to identify the
%    non-dominated fronts from all the design points
% 4. The pareto fronts of the test problems are calculated using the code
%    in the PlatEMO toolbox[4].
% 5. The EIM criteria are maximized by DE [5] algorithm.
% -----------------------------------------------------------------------------------------
% [1]  D. Zhan, Y. Cheng, J. Liu, Expected Improvement Matrix-based Infill
%      Criteria for Expensive Multiobjective Optimization. IEEE Transactions
%      on Evolutionary Computation, 2017, 21 (6): 956-975.
% [2] Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging
%     Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%     Modelling, Technical University of Denmark, 2002. Available at:
%     http://www2.imm.dtu.dk/~hbn/dace/.
% [3] http://www.mathworks.com/matlabcentral/fileexchange/17251-
%      pareto-front.
% [4] N. Beume, C.M. Fonseca, M. Lopez-Ibanez, L. Paquete, J. Vahrenhold,
%     On the Complexity of Computing the Hypervolume Indicator, IEEE
%     Transactions on Evolutionary Computation 13(5) (2009) 1075-1082.
% [4] Tian, Y., R. Cheng, X. Zhang, and Y. Jin. PlatEMO: A MATLAB
%      Platform for Evolutionary Multi-Objective Optimization. IEEE
%      Computational Intelligence Magazine,2017, 12 (4):73-87.
% [5] K. Price, R. M. Storn, and J. A. Lampinen, Differential evolution:
%     a practical approach to global optimization: Springer Science & Business Media, 2006.
%     http://www.icsi.berkeley.edu/~storn/code.html
% -----------------------------------------------------------------------------------------
% zhandawei@swjtu{dot}edu{dot}cn
% 2017.05.03 initial creation
% 2018.03.19 update
% 2018.09.18 update
% 2018.11.28 update, use DE optimizer for finding EIM maximum
% 2019.04.03 update, calcuate the IGD value in each iteration
% -----------------------------------------------------------------------------------------
clearvars;close all;
addpath('test_problems','results')
seed = 1;
rng(seed, 'twister');
% settings of the problem
% ZDT1-4, ZDT6, DTLZ1-7
fun_name = 'ZDT3';
% number of objectives (should be 2 for ZDT problems)
num_obj = 2;
% number of design variables
num_vari = 6;
% infill criterion: 'EIM_Euclidean','EIM_Maximin','EIM_Hypervolume'
infill_name= 'EIM_Hypervolume';
%-------------------------------------------------------------------------
% number of initial design points
num_initial = 11*num_vari-1;
% the maximum allowed evaluations
max_evaluation = 100;
%-------------------------------------------------------------------------
% get the information about the problem
design_space = [zeros(1,num_vari);ones(1,num_vari)];
% the pareto fronts of the problem for calculating the IGD value
pareto_front = Calculate_Pareto_Front(fun_name, 10000, num_obj);
pareto_front = readtable('zdt3front.txt' );
pareto_front = pareto_front{:,:};
%-------------------------------------------------------------------------
% the intial design points, points sampled all at once
sample_x = repmat(design_space(1,:),num_initial,1) + repmat(design_space(2,:)-design_space(1,:),num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);

sample_y = feval(fun_name, sample_x, num_obj);
% scale the objectives to [0,1]
sample_y_scaled = (sample_y - repmat(min(sample_y),size(sample_y,1),1))./repmat(max(sample_y)-min(sample_y),size(sample_y,1),1);
%-------------------------------------------------------------------------
% initialize some parameters
evaluation = size(sample_x,1);
kriging_obj = cell(1,num_obj);
iteration = 0;
IGD = zeros(max_evaluation-num_initial+1,1);
ref_point = 1.1*ones(1,2);
hypervolume = zeros(max_evaluation-num_initial+1,1);
%-------------------------------------------------------------------------
% identify the non-dominated front
index = Paretoset(sample_y);
non_dominated_front = sample_y(index,:);
non_dominated_front_scaled = sample_y_scaled(index,:);

IGD(1) = mean(min(pdist2(pareto_front,non_dominated_front),[],2));
hypervolume(1) = Hypervolume(non_dominated_front,ref_point);

% print the IGD information
fprintf('----------------------------------------------------------------\n')
fprintf(' iteration: %d, evaluation: %d, , IGD: %0.4g \n', iteration, evaluation, IGD(1));
%-------------------------------------------------------------------------
% beginning of the iteration

eim_opt = [];
propose_x = [];
while evaluation <  max_evaluation
    % build the initial kriging model for each objective
    for ii=1:num_obj
        kriging_obj{ii} = dacefit(sample_x,sample_y_scaled(:,ii),'regpoly0','corrgauss',1*ones(1,num_vari), 0.001*ones(1,num_vari),1000*ones(1,num_vari));
    end
        
    % select updating points using the EIM criteria
    switch infill_name
        case 'EIM_Euclidean'
            infill_criterion = @(x)Infill_EIM_Euclidean(x, kriging_obj, non_dominated_front_scaled);
        case 'EIM_Maximin'
            infill_criterion = @(x)Infill_EIM_Maximin(x, kriging_obj, non_dominated_front_scaled);
        case 'EIM_Hypervolume'
            infill_criterion = @(x)Infill_EIM_Hypervolume(x, kriging_obj, non_dominated_front_scaled);
        otherwise
            error('you should select infill_name from EIM_Euclidean, EIM_Maximin, and EIM_Hypervolume');
    end
    
    
    
    [bestmem1, bestval1, nfeval1] = DE(infill_criterion, num_vari, design_space(1,:), design_space(2,:), 50, 50);
    [bestmem2, bestval2, nfeval2] = DE(infill_criterion, num_vari, design_space(1,:), design_space(2,:), 50, 50);
    [bestmem3, bestval3, nfeval3] = DE(infill_criterion, num_vari, design_space(1,:), design_space(2,:), 50, 50);
    [bestmem4, bestval4, nfeval4] = DE(infill_criterion, num_vari, design_space(1,:), design_space(2,:), 50, 50);
    
    bests = [bestmem1; bestmem2; bestmem3; bestmem4];
    bestval = [bestval1, bestval2, bestval3, bestval4];
    eim_opt = [eim_opt, min(bestval)]; 
    
    a = find(bestval == min(bestval));
    best_x = bests(a(1),:);
    if size(a)>1
        fprintf('multiple min f')
    end
    % csvwrite('best_x.csv', best_x)
    % y = Infill_EIM_Hypervolume(best_x, kriging_obj, non_dominated_front_scaled);
    propose_x = [propose_x; best_x];
    
    % add the new points to the design set
    if min(sum(sample_x - best_x, 2)) == 0
        same_x = sample_x;
        sample_y = sample_y;
        fprintf('got..')
    else
        sample_x = [sample_x;best_x];
        sample_y = [sample_y; feval(fun_name,best_x, num_obj)];        
    end
    
    
    
    sample_y_scaled = (sample_y - repmat(min(sample_y),size(sample_y,1),1))./repmat(max(sample_y)-min(sample_y),size(sample_y,1),1);
    evaluation = evaluation + size(best_x,1);
    iteration = iteration + 1;
    % update the non-dominated front
    index = Paretoset(sample_y);
    
    non_dominated_front = sample_y(index,:);
    non_dominated_front_scaled = sample_y_scaled(index,:);
    % csvwrite('non_dominated_front.csv', non_dominated_front)
    
    if num_obj == 2
        % scatter(sample_y(:,1), sample_y(:,2),'ro');
        clf('reset');
        scatter(pareto_front(:, 1), pareto_front(:, 2),'bo'); hold on;
        scatter(non_dominated_front(:,1), non_dominated_front(:,2),'ro', 'filled');        
        title(sprintf('iteration: %d, evaluations: %d',iteration,evaluation));drawnow;   
    end
    
   
    % following two lines only for test
    % pareto_front = readtable('zdt3front.txt' );
    % pareto_front = pareto_front{:,:};
    % normalized to 0-1 according to pf front
    pf_scaled = (pareto_front - repmat(min(pareto_front), size(pareto_front,1), 1))./ repmat(max(pareto_front)-min(pareto_front), size(pareto_front, 1), 1);
    nd_scaled = (non_dominated_front - repmat(min(pareto_front), size(non_dominated_front,1), 1) )./ repmat(max(pareto_front)-min(pareto_front), size(non_dominated_front, 1), 1);
    
    IGD(iteration+1) = mean(min(pdist2(pf_scaled,nd_scaled),[],2));
    hypervolume(iteration + 1) = Hypervolume(nd_scaled,ref_point);
    
    hv = hypervolume(iteration + 1);
    igd =  IGD(iteration+1);
    
    % IGD(iteration+1) = mean(min(pdist2(pareto_front,non_dominated_front),[],2));
    % hypervolume(iteration + 1) = Hypervolume(non_dominated_front,ref_point);
    % print the IGD information
    fprintf(' iteration: %d, evaluation: %d, hv: %0.4g, IGD: %0.4g \n', iteration, evaluation, hypervolume(iteration + 1) , IGD(iteration+1));
end


m = mean(eim_opt);
s = std(eim_opt);
fprintf('eim optimization, mean %.4g, std %.4g first %.4g \n', m, s, eim_opt(1));
csvwrite('propose_x.csv', propose_x);

out = [hv, igd];
fname = string(seed) + 'out.csv';
csvwrite(fname, out);

% save the data of the run
%save('../results/data.mat','sample_x','sample_y','IGD');
% plot the iteration history with respect to the IGD value
%figure;
%plot((0:max_evaluation - num_initial)',IGD,'ro-');
%xlabel('iteration');ylabel('IGD');
%title(sprintf('EIM-%s-EGO on %d-d %d-objective %s problem',infill_name(5:end),num_vari,num_obj,fun_name));
%saveas(gcf, '../results/IGD_interation_history.png');
%saveas(gcf, '../results/IGD_interation_history.fig');
% plot the final non-dominated front points
pareto_front = readtable('zdt3front.txt' );
pareto_front = pareto_front{:,:};
if num_obj == 2
    figure;
    
    scatter(pareto_front(:, 1), pareto_front(:, 2),'d')
    hold on
    scatter(non_dominated_front(:,1), non_dominated_front(:,2),'ro');
    xlabel('f1');ylabel('f2');
    title(sprintf('EIM-%s-EGO on %d-d %d-objective %s problem',infill_name(5:end),num_vari,num_obj,fun_name));
    s = 'non_dominated_solutions'+string(seed) + '.png';
    saveas(gcf, s);
    %    saveas(gcf, '../results/non_dominated_solutions.fig');
elseif num_obj == 3
    figure;
    scatter3(non_dominated_front(:,1), non_dominated_front(:,2),non_dominated_front(:,3),'ro');
    xlabel('f1');ylabel('f2');zlabel('f3');
    title(sprintf('EIM-%s-EGO on %d-d %d-objective %s problem',infill_name(5:end),num_vari,num_obj,fun_name))
    saveas(gcf, '../results/non_dominated_solutions.png');
    saveas(gcf, '../results/non_dominated_solutions.fig');
end
