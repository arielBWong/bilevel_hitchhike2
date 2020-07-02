classdef smd9
    properties
        p = 1;
        q = 2;
        r = 1;
        n_lvar;
        n_uvar;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        name;
    end
    methods
        function obj = smd9(p, q, r)
            if nargin == 3
                obj.p = p;
                obj.q = q;
                obj.r = r;
            else
                obj.p = 1;
                obj.q = 2;
                obj.r = 1;
            end
            obj.name = 'SMD9';

            % level variables
            obj.n_lvar = obj.q + obj.r;
            obj.n_uvar = obj.p + obj.r;
            
            %init bound upper level
            xu_bl_1 = ones(1, obj.p)* (-5.0);
            xu_bu_1 = ones(1, obj.p) * 10.0;
            xu_bl_2 = ones(1, obj.r)* (-5.0);
            xu_bu_2 = ones(1, obj.r) * 1.0;
            obj.xu_bl = [xu_bl_1, xu_bl_2];
            obj.xu_bu = [xu_bu_1, xu_bu_2];
            
            % init bound lower level
            xl_bl_1 = ones(1, obj.q) * (-5.0);
            xl_bu_1 = ones(1, obj.q) * 10.0;
            xl_bl_2 = ones(1, obj.r) * (-1.0 + 1e-10);
            xl_bu_2 = ones(1, obj.r) * (-1.0 + exp(1));
            obj.xl_bl = [xl_bl_1, xl_bl_2];
            obj.xl_bu = [xl_bu_1, xl_bu_2];
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            xu1 = xu(1:obj.p);
            xu2 = xu(obj.p+1:obj.p+obj.r);
            
            xl1 = xl(1:obj.q);
            xl2 = xl(obj.q+1:obj.q+obj.r);
            
            functionValue = sum((xu1).^2) ...
                - sum((xl1).^2) ...
                + sum((xu2).^2) - sum((xu2 - log(1+xl2)).^2);
            
            f = functionValue;
            
            inequalityConstrVals(1) = sum(xu1.^2)+sum(xu2.^2) - floor(sum(xu1.^2)+sum(xu2.^2)+0.5);
            c = - inequalityConstrVals;
            
        end
        
        function[f, c] = evaluate_l(obj, xu, xl)
            xu1 = xu(:, 1:obj.p);
            xu2 = xu(:, obj.p+1:obj.p+obj.r);
         
            xl1 = xl(:, 1:obj.q);
            xl2 = xl(:, obj.q+1:obj.q+obj.r);
            
            functionValue = sum((xu1).^2, 2) ...
                + sum((xl1).^2, 2) ...
                + sum((xu2 - log(1+xl2)).^2, 2);
            
            f = functionValue;
            
            
            %Write the constraints here
            inequalityConstrVals = sum(xl1.^2, 2)+sum(xl2.^2, 2) - floor(sum(xl1.^2, 2)+sum(xl2.^2, 2)+0.5);
            c = - inequalityConstrVals;
            
        end
    end
end