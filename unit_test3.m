%%
% This unittest is for 1 objectives 2 variables and 2 constraints
clearvars;
close all;
% (1) intialize training data
% (2) finding xl for initialization data
%     (2-1) for each xu find xl. input: xu, problem, outupt: match_xl
% (3) use xu and fu build surrogate 
% (4) eim propose next newxu
%      (4-1) match xl for xu
% (5) conduct bilevel local search
% (6) conduct final hybrid ll search
