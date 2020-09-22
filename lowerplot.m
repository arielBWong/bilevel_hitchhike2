prob = smd9(1, 1, 1);
xu = [0, 0];
xl1 = linspace(-5, 10, 100);
xl2 = linspace(-1.0 + 1e-10, -1.0 + exp(1), 100);

 [xl1, xl2] = meshgrid(xl1, xl2);
 
  f = zeros(100, 100);
 for i = 1:100
     for j = 1:100
%   f = sum((xu1).^2, 2) ...
%                 + sum((xl1).^2, 2) ...
%                 + sum((xu2 - log(1+xl2)).^2, 2);
%             
             f (i, j) = xu(1)^2 +  xl1(i, j).^2  +  (xu(2) - log(1+ xl2(i, j))).^2;
       
     end
 end
 
 min(f(:))
surf(xl1, xl2, f);
xlabel('xl1','FontSize', 16);
ylabel('xl2', 'FontSize', 16);
zlabel('fl', 'FontSize', 16);
colormap jet
shading interp
