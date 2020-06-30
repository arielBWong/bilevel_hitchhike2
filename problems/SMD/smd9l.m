
function [functionValue, equalityConstrVals, inequalityConstrVals] = smd9l(xu,xl)

    r = floor(length(xu)/2);
    p = length(xu) - r;
    q = length(xl) - r;
    
    xu1 = xu(1:p);
    xu2 = xu(p+1:p+r);

    xl1 = xl(1:q);
    xl2 = xl(q+1:q+r);

    functionValue = sum((xu1).^2) ...
                    + sum((xl1).^2) ...
                    + sum((xu2 - log(1+xl2)).^2);

    functionValue = -functionValue;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the constraints here
    inequalityConstrVals(1) = sum(xl1.^2)+sum(xl2.^2) - floor(sum(xl1.^2)+sum(xl2.^2)+0.5);
    inequalityConstrVals = - inequalityConstrVals;
    equalityConstrVals = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%