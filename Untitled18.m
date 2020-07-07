%% number calcuation to think about whether nomalization on nd front is monotonic

ref_point = [1.1, 1.1];
f = [1, 0];
f = [f; 1/2.0, 1.];
% f = [f; 3/4.0, 1/2.0];
n = size(f, 1);
f_norm = (f - (min(f)))./(max(f) - min(f));
f_norm =(f - repmat(min(f), n, 1))./(repmat(max(f), n, 1) - repmat(min(f), n, 1));
hv1 = Hypervolume(f_norm,ref_point); 
disp(hv1);

f = [f; 0, 10];
n = size(f, 1);
f_norm =(f - repmat(min(f), n, 1))./(repmat(max(f), n, 1) - repmat(min(f), n, 1));
hv2 = Hypervolume(f_norm,ref_point); 
disp(hv2);
