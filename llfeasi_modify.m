function fu = llfeasi_modify(fu, feasi_list, ind)  % unit test needed
% this function is to varify fu wrt the feasibility of xl
% if xl is infeasible, fu is modified to a higher value
% so that this point is not preferred in later search
% consider range(1:ind) to deal with both one instance and
% a list of instances
% input 
%   fu                                   : list of original fu values
%                                                   1/2d column form
%   feasi_list                          : list of flags refering whether
%                                                   matching xl is feasible (true/false)
%   ind                                : indicator(int) which row/instance we are  dealing with
% ouput
%   fu                                  : possibly modified fu
%-----------------------------------------------------------------------------------
% use those feasible fu before current match
% to determine modified fu value
feasi_xl = feasi_list(1:ind-1) == true; 
if sum(feasi_xl) >0                                         % if there eixts feasible solution in previous list to set modified fu                                                                      
    if ~feasi_list(ind)                                         % current flag is false
        f = max(fu(feasi_xl, :),[], 1);                         % multiple objective compatible
        f =  f + 1;
        % all previous infeasible fu needs update
        % modify_ind = feasi_list(1:ind) == false;
        % fu(modify_ind, :) = f;
        %---compatible with mo assignment
        for k = 1:ind
            if feasi_list(k) == false
                fu(k, :) = f;
            end
        end
        %------------------------------------
    end
end
end