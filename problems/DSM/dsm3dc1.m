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
        r;
    end
    methods
        function obj = dsm3dc1(k1, k2)
            obj.r = 0.5;
            obj.p = k1;
            obj.q = k2;
            
            % level variables
            obj.n_lvar = obj.q;
            obj.n_uvar = obj.p;
            
            % bounds
            % init bound upper level
            obj.xu_bl = [0, ones(1, obj.p-1) * (-obj.p)];
            obj.xu_bu = [0.5, ones(1, obj.p-1) * obj.p];
            
            
            % init bound lower level
            obj.xl_bl = ones(1, obj.q) * (-obj.q);
            obj.xl_bu = ones(1, obj.q) * obj.p;
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            %-obj
            obj.r = 0.1;
            tao = 1;
            
            p3 = tao* sum((xl(:, 2:obj.n_lvar) - xu(:, 2:obj.n_uvar)) .^ 2, 2);
            
            p2 = 2: obj.n_lvar;
            p2 =( p2 - 1) /2;
            p2 =  sum((xu(:, 2:obj.n_lvar) - p2) .^2 , 2);
            
            p1 = pfshape_line(xu, obj.r);
             
            f(:, 1) = p1(:, 1) + p2 + p3 ;
            f(:, 2) = p1(:, 2) + p2 + p3;
            
            
            %-cie
            c = constraint1_u(xu);
            
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            
            p2 = sum(( xl - xu) .^2, 2);
            %-obj
            p3 = 1 * abs(sin(pi/obj.n_lvar .* (xl(:, 2:obj.n_lvar) - xu(:, 2:obj.n_uvar))));
            f(:, 1) = p2 + sum(p3, 2);
            
            %-cie
            c = constraint1_l(xu, xl);
            
        end
        
        function pf = upper_pf(obj, num_point)
           
            [pf, ~] = UniformPoint(num_point,2);
            pf = pf * 0.5   - [0, 0.4];
            
            
        end
    end
end
