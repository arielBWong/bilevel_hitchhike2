classdef Gpc
    properties
        n_var;
        xl;
        xu;
        n_con;
        n_obj;
        ref;
        name;
    end
    methods
        function obj = Gpc()
            obj.n_var = 2;
            obj.xl = [-2, -2];
            obj.xu = [2, 2];
            obj.n_con = 2;
            obj.n_obj = 1;
            obj.ref = [];
            obj.name = 'Gpc';
        end
        function [f, con] = evaluate(obj, x)
            x1 = x(:,1); x2 = x(:, 2);
            
            A = 19 - 14 * x1 + 3 * x1.^2 - 14 * x2 + 6 * x1 .* x2 + 3 * x2  .^ 2;
            B = 18 - 32 * x1 + 12 * x1.^2 + 48 * x2 - 36 * x1 .* x2 + 27 * x2 .^ 2;
            
            f = (1 + A .* (x1 + x2 + 1).^2) * (30 + B .* (2 * x1 - 3 * x2).^2);
            
            g1 = -3 * x1 + (-3 * x2).^3;
            g2 = x1 - x2 - 1;
            
            con = [g1, g2];
            
        end
        
        function [f, con] = evaluate_l(obj, xu, x)
            x1 = x(:,1); x2 = x(:, 2);
            
            A = 19 - 14 * x1 + 3 * x1.^2 - 14 * x2 + 6 * x1 .* x2 + 3 * x2  .^ 2;
            B = 18 - 32 * x1 + 12 * x1.^2 + 48 * x2 - 36 * x1 .* x2 + 27 * x2 .^ 2;
            
            f = (1 + A .* (x1 + x2 + 1).^2) * (30 + B .* (2 * x1 - 3 * x2).^2);
            
            g1 = -3 * x1 + (-3 * x2).^3;
            g2 = x1 - x2 - 1;
            
            con = [g1, g2];
            
        end
    end
end