classdef Mystery
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
        function obj = Mystery()
            obj.n_var = 2;
            obj.xl = [0, 0];
            obj.xu = [5, 5];
            obj.n_con = 1;
            obj.n_obj = 1;
            obj.ref = [];
            obj.name = 'Mystery';
        end
        function [f, con] = evaluate(obj, x)
            x1 = x(:,1); x2 = x(:, 2);
            
            part1 = 0.01 * (x2 - x1.^2).^2;
            part2 = (1 - x1).^2;
            part3 = 2 * (2 - x2).^2;
            part4 = 7 * sin(0.5 * x1) * sin(0.7 .* x1 .* x2);

            f = 2 + part1 + part2 + part3 + part4;

            con = -sin(x1 - x2 - pi/8);
     
        end
    end
end