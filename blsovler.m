function [xu_end, xl_end, n_up, n_low] = blsovler(prob, xu_start, num_pop, num_gen, inisize_l, numiter_l)
% this function runs local search on upper level problem
% usage
% input     prob        :bilevel problem
%           xu_start    :starting point
%
% output    xu_end      :optimization result
%--------------------------------------------------------------------------
global xu_g
global xl_g
global ll_n
ll_n = 0;

fmin_obj = @(x)blobj(x, prob, num_pop, num_gen, inisize_l, numiter_l);
fmin_con = @(x)blcon(x, prob, num_pop, num_gen, inisize_l, numiter_l);
opts = optimset('fmincon');
opts.Algorithm = 'sqp';
opts.MaxFunctionEvaluations = 20;
[xu_end, endfu, flag, output] = fmincon(fmin_obj, xu_start, [], [],[], [],  ...
    prob.xu_bl, prob.xu_bu, fmin_con,opts);

xl_end = check_exist(xu);



clear global
end

function fu =  blobj(xu, prob, num_pop, num_gen, inisize_l, numiter_l)
[xl, n, ~] = llmatch(xu, prob,num_pop, num_gen,inisize_l, numiter_l);
fprintf('obj is called %d', n);
[fu, ~] = prob.evaluate_u(xu, xl);
end

function [c, ceq] = blcon(xu, prob,  num_pop, num_gen, inisize_l, numiter_l)
[xl, n, ~] = llmatch(xu, prob,num_pop, num_gen,inisize_l, numiter_l);
fprintf('con is called %d', n);
[~, c] = prob.evaluate_u(xu, xl);
ceq = [];
end

function xl = check_exist(xu)
nvar = size(xu, 2);
xu = round(xu, 10);
xu_g = round(xu_g, 10);
diff = xu_g - xu;
ind = sum(diff, 2) == 0;

if sum(ind)>0 % should be only 1, if exists
    fprintf('found xu in global save %d', sum(ind));
    xl = xl_g(ind);
    if sum(ind)>1
        error('there cannot be repeat x more than once')
    end
else
    xl = [];
end
end