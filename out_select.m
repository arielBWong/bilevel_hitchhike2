function [best_x,best_f, best_c, s, index] = out_select(xu,  xl, prob)
% this function select a starting point for a local search algorithm
% usage
% input:
%           x                       : 2d matrix of design variables,
%                                              (original range)
%           
% output:
%           best_x                  : one x instance with smallest objective value
%                                               or with smallest feasible objective value
%                                               or with smallest const value if no
%                                                               feasible solution exists
%           best_f                  : corresponding obj value of best_x
%           best_c                  : corresponding con value of best_x
%           s                       : flag indicate whether there is feasible solution
%           index                   : index of output x in archive
%--------------------------------------------------------------------------
[uf, uc] = prob.evaluate_u(xu, xl);
[lf, lc]= prob.evaluate_l(xu, xl); 
c = [uc, lc];

if ~isempty(c)              % constraint problem
     num_con = size(c, 2);
        index_c = sum(c <= 0, 2) == num_con;
        if sum(index_c) == 0 % no feasible, return f with smallest constraints
            sum_c = sum(c, 2);
            [~, i] = min(sum_c);
            best_x = xu(i, :);
            best_f = uf(i, :);
            best_c = uc(i, :);
            s = false;
            index = i;
        else % has feasible, return feasible smallest f
            feasi_ff = uf(index_c, :);
            feasi_x = xu(index_c, :);
            feasi_fc = uc(index_c, :);
            [~, i] = min(feasi_ff);
            best_x = feasi_x(i, :);
            best_f = feasi_ff(i, :);
            best_c = feasi_fc(i, :);
            s = true;
            [~, index] = ismember(best_x, xu);
        end
else                        % unconstraint problem
       [best_f, i] = min(uf);
        best_x = xu(i, :);
        s = true;
        best_c = [];
        index = i;
    
end

end