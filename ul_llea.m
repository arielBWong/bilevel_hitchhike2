function ul_llea(prob, seed)
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

global xu_g
global xl_g
global ll_ln
xu_g =[];
xl_g =[];
ll_ln = 0;

funh_obj = @(x)obj(x, prob, num_pop, num_gen);
funh_con = @(x)con(xu, prob,  num_pop, num_gen);

gsolver(funh_obj, num_xvar, lb, ub, initmatrix, funh_con, param);

% xu _g and xl_g stores all the solutions in the whole process
% therefore, final results can use xu_g and xl_g to generate

finalresults_process(xu_g, xl_g, prob, seed)
end

function  fu = obj(xu, prob, num_pop, num_gen)
xl = check_exist(xu);
if isempty(xl)
     [xl, fev] = llmatch_ea(xu, prob,num_pop, num_gen);
     global xu_g
     global xl_g
     xu_g = [xu_g; xu];
     xl_g = [xl_g; xl];
      ll_ln = ll_ln + fev;
end
[fu, ~] = prob.evaluate_u(xu, xl);
end


function [c] = con(xu, prob,  num_pop, num_gen)
xl = check_exist(xu);

if isempty(xl)
    [xl, fev] = llmatch_ea(xu, prob,num_pop, num_gen);
    % fprintf('con llmatch is called %d\n', n);
    global xu_g
    global xl_g
    xu_g = [xu_g; xu];
    xl_g = [xl_g; xl];
    ll_ln = ll_ln + fev;
end
[~, c] = prob.evaluate_u(xu, xl);

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

[fu, cu] = prob.evaluate_u(xu_g, xl_g);
[~, cl] = prob.evaluate_l(xu_g, xl_g);

% if lowe level is constraint problem
% modify  fu according to lower feasibility
if ~isempty(cl) 
    numcon = size(cl, 2);
    numins = size(cl, 1);
    cl(cl<=0) = 0;
    feasi_flag = sum(cl, 2) == 0;
    for ii = 1:numins
            fu = llfeasi_modify(fu, feasi_flag, ii);
    end
end
perfrecord_umoc(xu_g, fu, cu, prob, seed, 'ea_ea');
end