k = 3;
xu2 = linspace(-k, k, 100);
xu1 = 0.5 - xu2.^2;
xu1(xu1 < 0) = 0;
plot(xu2, xu1);




% xu1 = linspace(0, 0.5, 100);
% f1 = 1.1 *(1-cos(pi * xu1));
% f2 = 1.1* (1-sin(pi * xu1));
% 
% plot(f1, f2);