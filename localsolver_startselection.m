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
%           best_x                  : one x instance with smallest objective value
%                                               or with smallest feasible objective value
%                                               or with smallest const value if no
%                                                               feasible solution exists
%           best_f                  : corresponding obj value of best_x
%           best_c                  : corresponding con value of best_x
%           s                       : flag indicate whether there is feasible solution
%--------------------------------------------------------------------------
if size(ff, 2) == 1 % single objective
    if isempty(fc) % non constraint problem
        [best_f, i] = min(ff);
        best_x = x(i, :);
        s = true;
        best_c = [];
    else
        num_con = size(fc, 2);
        index_c = sum(fc <= 0, 2) == num_con;
        if sum(index_c) == 0 % no feasible, return f with smallest constraints
            sum_c = sum(fc, 2);
            [~, i] = min(sum_c);
            best_x = x(i, :);
            best_f = ff(i, :);
            best_c = fc(i, :);
            s = false;
        else % has feasible, return feasible smallest f
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
    % input ff is all ffs ulego have got, new starting point should be
    % the solution that contriubte most to hv compared to the rest of the ffs
    % it means this ff locates on a nd area where there are not many points
    % reference point should be just 1.1 in normalized space
    [i, s] = hvcontribution_selectf(ff, fc);
    best_x = x(i, :);
    best_f = ff(i, :);
    best_c = fc(i, :);
end
end


% auxiliary function on selecting f with most contribution on hv
function [bestf_id, s] = hvcontribution_selectf(ff, fc)
% this function selects one feasible f that contribute most in hv
% has compatibility with constrain and non constraint problems
% usage
%   input
%   output
%---------------------------------------------------------------------------
num_con = size(fc, 2);
num_obj = size(ff, 2);
num_ff = size(ff, 1);
ref_point = ones(1, num_obj) * 1.1;

if num_con == 0
    ff_norm = (ff -repmat(min(ff), num_ff, 1))./repmat(max(ff)-min(ff),num_ff,1);
    nd_index = Paretoset(ff_norm);
    nd_front = ff_norm(nd_index, :);
    ff_max =  maxcontribution_select(nd_front, ref_point);
    % max_id's original ff location
    ff_norm = sum(ff_norm - ff_max, 2);
    bestf_id = find(ff_norm < 1e-10);
    s = true;
    return
end
%- consider constrainted problem
index_c = sum(fc <= 0, 2) == num_con;
if sum(index_c) ~=0   % feasible solution exists 
    % normalize with all ff, reference point is set to 1.1
    ff_norm = (ff -repmat(min(ff), num_ff, 1))./repmat(max(ff)-min(ff),num_ff,1);
    
    feasi_fnorm = ff_norm(index_c, :);
    nd_index = Paretoset(feasi_fnorm);
    nd_front = feasi_fnorm(nd_index, :);
    ff_max =  maxcontribution_select(nd_front, ref_point);
    
    % max_id's original ff location
    ff_norm = sum(ff_norm - ff_max, 2);
    bestf_id = find(ff_norm < 1e-10);
    s = true;
else
    % feasible solution does not exist
    % choose smallest feasibility
    sum_c = sum(fc, 2);
    [~, bestf_id] = min(sum_c);   
    s = false;
end

end


function ff_max = maxcontribution_select(nd_front, ref_point)
% this function takes as input a list of nd_front, and a reference point
% gives back a ff instance, which singlely contributes most to hv values
% usage
%   input
%   ouput
%--------------------------------------------------------------------------
num_nd = size(nd_front, 1);
hv_leaveoneout = zeros(num_nd, 1);

hv_fullhouse = Hypervolume(nd_front,ref_point);

for i = 1:num_nd
    leave1out_nd = nd_front;
    leave1out_nd(i, :) = [];
    hv_leaveoneout(i) = Hypervolume(leave1out_nd,ref_point);
end

hv_contri = hv_fullhouse - hv_leaveoneout;

[~, max_id] = max(hv_contri);
ff_max = nd_front(max_id, :);
end