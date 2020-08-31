r = 0.1;
num_point = 100;
sep = pi/(2 *(num_point-1));
pf = [0,  (1+r)];

deg = 0;
for i = 1:num_point-1
    one = [(1+r) * sin(deg + i * sep), (1+r) * cos(deg + i * sep)];
    pf = [pf; one];
end
