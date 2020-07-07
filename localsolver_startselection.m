function [best_x,best_f, best_c, s] = localsolver_startselection(x,  ff, fc)
% this function select a starting point for a local search algorithm
% usage
% input:
%           x                       : 2d matrix of design variables,
%                                              (original range)
%           ff                      : objective function values of x
%                                              (original range)
%           fc                      : constraint values of x
%                                              (original range)
% output:
%           best_x              : one x instance with smallest objective value
%                                               or with smallest feasible objective value
%                                               or with smallest const value if no
%                                                               feasible solution exists
%           best_f              : corresponding obj value of best_x
%           best_c              : corresponding con value of best_x
%           s                       : flag indicate whether there is feasible solution
%--------------------------------------------------------------------------
if size(ff, 2)> 1
    error('this function does not compatible with mo problem')
end
if size(ff, 2) == 1 % single objective
    if isempty(fc) % non constraint problem
        [best_f, i] = min(ff);
        best_x = x(i, :);
        s = true;
        best_c = [];
    else
        num_con = size(fc, 2);
        index_c = sum(fc <= 0, 2) == num_con;
        if sum(index_c) == 0 % no feasible
            sum_c = sum(fc, 2);
            [~, i] = min(sum_c);
            best_x = x(i, :);
            best_f = ff(i, :);
            best_c = fc(i, :);
            s = false;
        else % has feasible
            feasi_ff = ff(index_c, :);
            feasi_x = x(index_c, :);
            feasi_fc = fc(index_c, :);
            [~, i] = min(feasi_ff);
            best_x = feasi_x(i, :);
            best_f = feasi_ff(i, :);
            best_c = feasi_fc(i, :);
            s = true;
        end
    end
else % multiple objective
    % for multiple objective 
    % starting point should be ff that has best hv contribution
end
end