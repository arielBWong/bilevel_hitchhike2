function empirical_pfc(problem, seed)


global xu_g
global xl_g
global ll_ln
xu_g =[];
xl_g =[];
ll_ln = 0;


tic;
fprintf(problem);

prob = eval(problem);
lnumgen =50;
lnumpop = 50;


funh_obj = @(x)obj(x, prob, lnumpop,lnumgen);
funh_con = @(x)con(x, prob,  lnumpop, lnumgen);

num_xvar = size(prob.xu_bl, 2);
initmatrix = [];
param = struct();
param.gen= 100;
param.popsize = 100;
gsolver(funh_obj, num_xvar, prob.xu_bl, prob.xu_bu, initmatrix, funh_con, param,prob,xu_g, xl_g);

% xu _g and xl_g stores all the solutions in the whole process
% therefore, final results can use xu_g and xl_g to generate
finalresults_process(xu_g, xl_g, prob, seed);
clear global
toc;
end



function  fu = obj(xu, prob, num_pop, num_gen)

global xu_g
global xl_g
global ll_ln

fu = [];
n = size(xu, 1);
for i = 1:n
    xu_single = xu(i, :);
    xl_single = check_exist(xu_single);
    if isempty(xl_single)
        [xl_single, fev] = llmatch_ea(xu_single, prob,num_pop, num_gen);
        
        xu_g = [xu_g; xu_single];
        xl_g = [xl_g; xl_single];
        ll_ln = ll_ln + fev;
    end
    [fu_single, ~] = prob.evaluate_u(xu_single, xl_single);
    fu = [fu; fu_single];
end
end


function [c] = con(xu, prob,  num_pop, num_gen)

global xu_g
global xl_g
global ll_ln

n = size(xu, 1);
c = [];
for i = 1:n
    xu_single = xu(i, :);
    xl_single = check_exist(xu_single);
    if isempty(xl_single)
        [xl_single, fev] = llmatch_ea(xu_single, prob,num_pop, num_gen);
        % fprintf('con llmatch is called %d\n', n);
        xu_g = [xu_g; xu_single];
        xl_g = [xl_g; xl_single];
        ll_ln = ll_ln + fev;
    end
    [~, c_single] = prob.evaluate_u(xu_single, xl_single);
    c = [c; c_single];
end
end



function xl = check_exist(xu)
%
xu = round(xu, 10);
global xu_g
global xl_g
if isempty(xu_g)
    xl = [];
else
    xu_g = round(xu_g, 10);
    diff = xu_g - xu;
    ind = sum(diff, 2) == 0;
    
    if sum(ind)>0           % should be only 1, if exists
        % fprintf('found xu in global save %d\n', sum(ind));
        
        if sum(ind)>1
            xl = xl_g(ind(1), :);
            fprintf('there is repeated x more than once');
        else
            xl = xl_g(ind, :);
        end
    else
        xl = [];
    end
end
end



function finalresults_process(xu_g, xl_g, prob, seed)
% this function extract the nd front for upper mo or best solution for
% upper so
% all assume lower problem is so
% usage
%   xu_g; global variable in main process for saving xu found
%   xl_g: global variable in main process for saving xl found
%  prob: problem instance
% output
%   best solutions saved in result folder
global ll_ln
[fu, cu] = prob.evaluate_u(xu_g, xl_g);
[~, cl] = prob.evaluate_l(xu_g, xl_g);

% if lowe level is constraint problem
% modify  fu according to lower feasibility
if ~isempty(cl)
    numins = size(cl, 1);
    cl(cl<=0) = 0;
    feasi_flag = sum(cl, 2) == 0;
    for ii = 1:numins
        fu = llfeasi_modify(fu, feasi_flag, ii);
    end
end
nfev_u = size(xu_g, 1);
nfev_l = ll_ln +nfev_u * 20 * 10;
perfrecord_umoc(xu_g, fu, cu, prob, seed, 'emp_pf', nfev_u, nfev_l);
end

