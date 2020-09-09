% k = 3;
% xu2 = linspace(-k, k, 100);
% xu1 = 0.5 - xu2.^2;
% xu1(xu1 < 0) = 0;
% plot(xu2, xu1);




% xu1 = linspace(0, 0.5, 100);
% f1 = 1.1 *(1-cos(pi * xu1));
% f2 = 1.1* (1-sin(pi * xu1));
%
% plot(f1, f2);

r = 0.1;
num_point = 100;
sep = pi/(2 *(num_point-1));
pf = [0,  (1+r)];

deg = 0;
for i = 1:num_point-1
    one = [(1+r) *sin(deg + i * sep), (1+r) * cos(deg + i * sep)];
    pf = [pf; one];
end

scatter(pf(:, 1), pf(:, 2));