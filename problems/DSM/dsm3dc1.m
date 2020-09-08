classdef dsm3dc1
    properties
        p ;
        q;
        n_lvar;
        n_uvar;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        name = 'dsm3dc1';
        uopt = NaN;
        lopt = NaN; % double check needed
    end
    methods
        function obj = dsm3dc1(k)
            obj.p = k;
            obj.q = k;
            
            % level variables
            obj.n_lvar = obj.q;
            obj.n_uvar = obj.p;
            
            % bounds
            % init bound upper level
            obj.xu_bl = [0, ones(1, k-1) * (-k)];
            obj.xu_bu = [0.5, ones(1, k-1) * k];
            
            
            % init bound lower level
            obj.xl_bl = ones(1, obj.q) * (-k);
            obj.xl_bu = ones(1, obj.q) * k;
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            %-obj
            r = 0.1;
            tao = -1;
            
            p3 = tao* sum((xl(:, 2:obj.n_lvar) - xu(:, 2:obj.n_uvar)) .^ 2, 2);
            
            p2 = 2: obj.n_lvar;
            p2 =( p2 - 1) /2;
            p2 =  sum((xu(:, 2:obj.n_lvar) - p2) .^2 , 2);
            
            p1 = pfshape_line(xu, r);
             
            f(:, 1) = p1(:, 1) + p2 + p3 ;
            f(:, 2) = p1(:, 1) + p2 + p3;
            
            
            
            %-cie
            c = constraint1_u(xu, xl);
            
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            
            p2 = sum(( xl - xu) .^2, 2);
            %-obj
            p3 = 10 * abs(sin(pi/obj.n_lvar .* (xl(:, 2:obj.n_lvar) - xu(:, 2:obj.n_uvar))));
            f(:, 1) = p2 + sum(p3, 2);
            
            %-cie
            c = constraint1_l(xu, xl);
            
        end
        
        function pf = upper_pf(obj, num_point)
           
            sep = pi/(2 *(num_point-1));
            pf = [1, 0];
            
            deg = 0;
            for i = 1:num_point-1
                one = [cos(deg + i * sep), sin(deg + i * sep)];
                pf = [pf; one];
            end
            
            
        end
    end
end
