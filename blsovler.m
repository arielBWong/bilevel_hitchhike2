function [xu_end, xl_end, n_up, n_low] = blsovler(prob, xu_start, num_pop, ...
    num_gen, inisize_l, numiter_l, penalityf)
% this function performs local search on upper level problem
% a wrapper around lower ego is provided to local solver
% usage
% input     
%           prob                       :  bilevel problem
%           xu_start                  :  starting point
%           num_pop              :  lower ego parameter
%           num_gen               :  lower ego parameter
%           inisize_l                  :  lower ego parameter
%           numiter_l               :  lower ego parameter
%           penalty_f                : upper level f value for infeasible xl
%                                              match 
%
% output    
%           xu_end                  :optimization result on xu
%           xl_end                    :optimiation result on matching xl
%           n_up                       :count on number function evaluation on ul
%           n_low                     :total count of number function evaluation on ll
%--------------------------------------------------------------------------
global xu_g
global xl_g
global ll_n
xu_g =[];
xl_g =[];
ll_n = 0;

%local search with sqp
fmin_obj = @(x)blobj(x, prob, num_pop, num_gen, inisize_l, numiter_l, penalityf);
fmin_con = @(x)blcon(x, prob, num_pop, num_gen, inisize_l, numiter_l);
opts = optimset('fmincon');
opts.Algorithm = 'sqp';
opts.Display = 'off';
opts.MaxFunctionEvaluations = 20;
[xu_end, ~, ~, output] = fmincon(fmin_obj, xu_start, [], [],[], [],  ...
    prob.xu_bl, prob.xu_bu, fmin_con,opts);

%return xu with corresponding xl
xl_end = check_exist(xu_end);
n_up = output.funcCount;
n_low = ll_n;

clear global
end



function fu =  blobj(xu, prob, num_pop, num_gen, inisize_l, numiter_l, penaltyf)
% to improve efficiency check existing match
xl = check_exist(xu, prob);

if isempty(xl)
    [xl, n, flag] =  llmatch(xu, prob, num_pop, num_gen,'EIMnext',...
    numiter_l, 'EIM_eval', 1);

    global xu_g
    global xl_g
    global ll_n
    ll_n = ll_n + n;
    xu_g = [xu_g; xu];
    xl_g = [xl_g; xl];
    if ~flag  % infeasible 
        fu = penaltyf;
        return 
    end
end

[fu, ~] = prob.evaluate_u(xu, xl);

% fu should be modified according to feasibility of lower constraints

% following efficiency can be improved
[~, lc] = prob.evaluate_l(xu, xl);
global ll_n
ll_n = ll_n + 1;
lc(lc<1e-6) = 0;    % constraint tolerance
if sum(lc)>0         % cope with infeasible lower  
    fu = penaltyf;
end


end

function [c, ceq] = blcon(xu, prob,  num_pop, num_gen, inisize_l, numiter_l)
xl = check_exist(xu, prob);

if isempty(xl)
    [xl, n, ~] =  llmatch(xu, prob, num_pop, num_gen,'EIMnext',...
    numiter_l, 'EIM_eval', 1);
    % fprintf('con llmatch is called %d\n', n);
    global xu_g
    global xl_g
    global ll_n
    ll_n = ll_n + n;
    xu_g = [xu_g; xu];
    xl_g = [xl_g; xl];
end
[~, c] = prob.evaluate_u(xu, xl);
ceq = [];
end

function xl = check_exist(xu, prob)

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
            name = prob.name;
            fprintf('there cannot be repeat x more than once, prob, %s \n', name);
            nn = 1:size(xu_g, 1);
            archive_xug = nn(ind);
            ind = archive_xug(1);
        end
        xl = xl_g(ind, :);
    else
        xl = [];
    end
end
end