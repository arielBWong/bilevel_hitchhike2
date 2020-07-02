function [best_x,best_f, s] = localsolver_startselection(x, fc, ff)
% this function select a starting point for local search
% need documentation and unit test
% need to check compatibility with non constraint probs
if size(ff, 2)>0
    error('this function does not compatible with mo problem')
end
if isempty(fc) % non constraint problem
    [best_f, i] = min(ff);
    best_x = x(i, :);
    s = true;
else
    num_con = size(fc, 2);
    index_c = sum(fc <= 0, 2) == num_con;
    if sum(index_c) == 0 % no feasible
        sum_c = sum(fc, 2);
        [~, i] = min(sum_c);
        best_x = x(i, :);
        best_f = ff(i, :);
        s = false;
    else % has feasible
        feasi_ff = ff(index_c, :);
        feasi_x = x(index_c, :);
        [~, i] = min(feasi_ff);
        best_x = feasi_x(i, :);
        best_f = feasi_ff(i, :);
        s = true;
    end
end
end