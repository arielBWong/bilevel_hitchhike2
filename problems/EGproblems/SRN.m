classdef SRN
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
        function obj = SRN()
            obj.n_var = 2;
            obj.xl = [-20,-20];
            obj.xu = [20, 20];
            obj.n_con = 2;
            obj.n_obj = 2;
            obj.ref = [250, 50];
            obj.name = 'SNR';
        end
        function [f, con] = evaluate(obj, x)
            % ----------------------------------------------------------------------------
            %Deb, K.,Pratap, A.,Agarwal, S., et al. A fast and elitist multiobjective
            %genetic algorithm: NSGA-II. IEEE Transactions on Evolutionary Computation
            % 2002, 6 (2): 182-197.
            % x1, x2 = [-20,20].
            % ----------------------------------------------------------------------------
            x1 = x(:,1); x2 = x(:, 2);
            f(:, 1) = 2+ (x1 - 2).^2 + (x2 - 1).^2;
            f(:, 2) =9*x1 - (x2 - 1).^2;
            con(:, 1) = x1.^2 + x2.^2 -225;
            con(:, 2) = x1 - 3*x2 +10;
            
            
        end
    end
end
