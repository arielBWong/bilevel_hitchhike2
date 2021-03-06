function ul_llea(prob_str, seed)
% this is a comparison algorithm for ulego
% it uses true evaluations on both levels to search for optimum solution
% usage
%       prob, seed
% output
%       statistical results saved in csv file
%--------------------------------------------

%  algorithm steps
% (1) evolution on upper level, but fitness needs to be an evolution on
% lower level, external archive is needed to handle lower level matching
% search.
% (2) fitness
%--- this function should take a population of xu, and return its objective
% value
% ---this function should also return a constraint value
% ---the handling of lower level search needs
% (3) llmatch_ea
%-- llmatch _ea takes in a xu and returns a xl
%--
%
num_pop = 20;
num_gen = 20;

prob = eval(prob_str);

global xu_g
global xl_g
global ll_ln
xu_g =[];
xl_g =[];
ll_ln = 0;

funh_obj = @(x)obj(x, prob, 20,4);
funh_con = @(x)con(x, prob,  20, 4);

num_xvar = size(prob.xu_bl, 2);
initmatrix = [];
param = struct();
param.gen= 6;
param.popsize = 20;


gsolver(funh_obj, num_xvar, prob.xu_bl, prob.xu_bu, initmatrix, funh_con, param,prob,xu_g, xl_g);

% xu _g and xl_g stores all the solutions in the whole process
% therefore, final results can use xu_g and xl_g to generate
finalresults_process(xu_g, xl_g, prob, seed)
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
        xl = xl_g(ind, :);
        if sum(ind)>1
            error('there cannot be repeat x more than once')
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
nfev_l = ll_ln +nfev_u * 20 * 4;
perfrecord_umoc(xu_g, fu, cu, prob, seed, 'ea_ea', nfev_u, nfev_l);
end