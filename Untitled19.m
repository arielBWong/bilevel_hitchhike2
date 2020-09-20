k =3;
prob = dsm1dc2(k, k);
% xl21 = linspace(-k, -0.5, 50);
% xl22 = linspace(1.5, k, 50);
% xl2 = [xl21, xl22];
% xl2 = linspace(-k, k, 100);
% xu = [0, 0.5];
% xl1 = linspace(-k, k, 100);

xl2 = linspace(-k, k, 100);
xl3 = linspace(-k, k, 100); 
 [xl2, xl3] = meshgrid(xl2, xl3);
 f = zeros(100, 100);
 for i = 1:100
     for j = 1:100
         f2 = (xl2(i, j) + 3).^2 + 10 * abs(sin(pi/3*(xl2(i, j) + 3)));
         f3 =  (xl3(i, j) + 3 ).^2 + 10 * abs(sin(pi/3*(xl3(i, j) +3)));
         f (i, j)= f2 + f3;
     end
 end
min(f(:))
surf(xl2, xl3, f);
xlabel('xl2','FontSize', 16);
ylabel('xl3', 'FontSize', 16);
zlabel('fl', 'FontSize', 16);
colormap jet
shading interp

% 
% k=2;
% xl1 = linspace(-k, k, 100);
% 
% xl2 = linspace(-k, k, 100);
%  [xl1, xl2] = meshgrid(xl1, xl2);
%  
%   f = zeros(100, 100);
%    for i = 1:100
%      for j = 1:100
%         
%          f (i, j)=(xl1(i, j)-0)^2 +  (xl2(i, j) - 0.5).^2 + 10 * abs(sin(pi/2*(xl2(i, j)- 0.5)));
%       end
%    end
% surf(xl1, xl2, f);
% xlabel('xl1','FontSize', 16);
% ylabel('xl2', 'FontSize', 16);
% zlabel('fl', 'FontSize', 16);
% colormap jet
% shading interp

% 
% f =  (xl2 - 0.5).^2 + 10 * abs(sin(pi/3*(xl2- 0.5)));
% plot(xl2, f);


% [mg1, mg2] = meshgrid(xl1, xl2);
% f = zeros(100, 100);
% for i =1:100
%     for j = 1:100
%         f(i, j) = xl2 - 
%         
%     end
% end

% 
% surf(mg1, mg2, f);
% xlabel('xl1');
% ylabel('xl2');
% 

function f = ll_eval(xu, xl)

p2 = sum(( xl(:, 1) - xu(:, 1)) .^2, 2);
%-obj
% p3 = 1 * abs(sin(pi/obj.n_lvar .* (xl(:, 2:obj.n_lvar) - xu(:, 2:obj.n_uvar))));
p3  = ll_p3(xu, xl);
f(:, 1) = p2 + sum(p3, 2);
end

function p3 = ll_p3(xu, xl)
n_lvar = size(xl, 2);
n_uvar = size(xu, 2);
p3 = 10 * abs(sin(pi/n_lvar.* (xl(:, 2:n_lvar) - xu(:, 2:n_uvar))));
end
            