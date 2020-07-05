function [xu_end, xl_end, n_up, n_low] = blsovler(prob, xu_start, num_pop, ...
    num_gen, inisize_l, numiter_l)
% this function performs local search on upper level problem
% a wrapper around lower ego is provided to local solver
% usage
% input     
%           prob                       :bilevel problem
%           xu_start                  :starting point
%           num_pop              :lower ego parameter
%           num_gen               :lower ego parameter
%           inisize_l                  :lower ego parameter
%           numiter_l              :lower ego parameter
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
fmin_obj = @(x)blobj(x, prob, num_pop, num_gen, inisize_l, numiter_l);
fmin_con = @(x)blcon(x, prob, num_pop, num_gen, inisize_l, numiter_l);
opts = optimset('fmincon');
opts.Algorithm = 'sqp';
opts.Display = 'off';
opts.MaxFunctionEvaluations = 20;
[xu_end, endfu, flag, output] = fmincon(fmin_obj, xu_start, [], [],[], [],  ...
    prob.xu_bl, prob.xu_bu, fmin_con,opts);

%return xu with corresponding xl
xl_end = check_exist(xu_end);
n_up = output.funcCount;
n_low = ll_n;

clear global
end

function fu =  blobj(xu, prob, num_pop, num_gen, inisize_l, numiter_l)
% to improve efficiency check existing match
xl = check_exist(xu);

if isempty(xl)
    [xl, n, ~] = llmatch(xu, prob,num_pop, num_gen,inisize_l, numiter_l);

    global xu_g
    global xl_g
    global ll_n
    ll_n = ll_n + n;
    xu_g = [xu_g; xu];
    xl_g = [xl_g; xl];
end

[fu, ~] = prob.evaluate_u(xu, xl);
end

function [c, ceq] = blcon(xu, prob,  num_pop, num_gen, inisize_l, numiter_l)
xl = check_exist(xu);

if isempty(xl)
    [xl, n, ~] = llmatch(xu, prob,num_pop, num_gen,inisize_l, numiter_l);
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

function xl = check_exist(xu)

xu = round(xu, 10);

global xu_g
global xl_g
if isempty(xu_g)
    xl = [];
else
    xu_g = round(xu_g, 10);
    diff = xu_g - xu;
    ind = sum(diff, 2) == 0;
    
    if sum(ind)>0 % should be only 1, if exists
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