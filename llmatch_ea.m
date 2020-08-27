function[match_xl, n_fev] = llmatch_ea(xu, prob, num_pop, num_gen)
% this function use gsolver to find a matching xl for given xu
%   usage
%       input
%           prob:  bilevel problem
%           num_pop: ea parameter population size
%           num_gen: ea parameter generation size
%           xu: upper x to be matched
%       ouput
%           match_xl: upper level x for xu
%           n_fev: number of function evaluation in local search
%----------
funh_obj = @(x)llobj(x, xu, prob);
funh_con = @(x)llcon(x, xu,  prob);

param = struct;
param.gen = num_gen;
param.popsize = num_pop;
num_xvar = length(prob.xl_bl);
initmatrix =[];

[bestx, bestf, bestc, archive] = gsolver(funh_obj, num_xvar,  prob.xl_bl, prob.xl_bu, initmatrix, funh_con, param);
% check bestc is feasible or not
% 
% if size(bestf, 1) == 1
%     [bestx, bestf, bestc] = multimodal_handle(archive.pop_last.X, archive.pop_last.F, archive.pop_last.C);
%     if size(bestx,1)>1
%         fprintf('there is multiple minimum exists with same f but different x');
%     end
% end


num_con = size(bestc, 2);
bestc(bestc <= 0) = 0;
s = sum(bestc, 2) == 0;

% what if there is multiple local minimum ?


% follow with local solver
[match_xl, n_fev] = losolver(bestx, s, bestf, bestc, prob,  xu);

end

%objective wrapper for both local and global search
function f = llobj(xl, xu, prob)
[f, ~] = prob.evaluate_l(xu, xl);
end

%constraint wrapper for global search
function [c]  = llcon(xl, xu, prob)
[~, c] = prob.evaluate_l(xu, xl);
end

% exception handle: rear situation where there is multiple minimum
% and they are the same
% ** only applicable on single objective lower level problem
function [bestxs, bestfs, bestcs] = multimodal_handle(popx, popf, popc)
ndig = 10;
bestxs = popx(1, :);  bestfs = popf(1, :);
if ~isempty(popc)
    bestcs = popc(1, :);
else
    bestcs = [];
end
n= size(popx, 1);
% (1) loop through population except the first individual
bestf = round(popf(1,:), ndig);
for i = 2:n
    % (2) check whether current individual f is the same as previous one
    if round(popf(i,:),10) == bestf
        % (3) if same f , check whether same x,
        % add current x to bestxs, after round to 10 digs
        % if they are all different then, this x is a new minimum
        % steps of the following line: stack -round-unique-size
        stack = round([bestxs; popx(i,:)], ndig);
        if size(unique(stack, 'row'), 1)== size(bestxs, 1) + 1
            bestxs= [bestxs; popx(i, :)];
            bestfs = [bestfs; popf(i, :)];
            if ~isempty(popc)
                bestcs =[bestcs; popc(i, :)];
            end
        end
    else
        %(4) if not same f,  means same f check is done break and return
        % with current return list
        return;       
    end
end
end

