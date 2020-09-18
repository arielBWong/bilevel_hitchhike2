function[best_x, info] =Believer_nextupper(train_x, train_y, xu_bound, xl_bound, ...
    num_pop, num_gen, train_c, fitnesshn, normhn)
% method of using one by one infill to generate next point
% normalization on objectives
% compatible with multiple objective
% usage:
%
% input:
%        train_x                - design variables
%                                           1/2d array: (num_samples, num_varibles)
%        train_y                - objective values
%                                           1/2d array: (num_samples, num_objectives)
%        xu_bound          - upper bound of train_x
%                                           1d array
%        xl_bound           - lower bound of train_x
%                                           1d array
%        num_pop          - EIM optimization parameter
%        num_gen          - EIM optimization parameter
%        train_c                - constraints values
%                                           1/2d array: (num_samples, num_constraints)
%        fitnesshn           -fitness valuation handle (callback) for ea process of
%                                           proposing next poing
%         normhn             -normalization handle on y
% output:
%       best_x                 - proposed next x to be evaluated by EIM
%       info                     - returned information for functor caller to recreate
%                                   - or check information
%                                   - info.krg
%                                   - info.krgc
%                                   - info.train_xmean
%                                   - info.train_ymean
%                                   - info.train_xstd
%                                   - info.train_ystd
%                                   - info.info.eim_normf
%--------------------------------------------------------------------------
