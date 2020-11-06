function test2
clear all;close all;clc;
x_trg= [0, 0.5, 1]'; no_trials=10; no_samples=20;
f_trgo=forrester(x_trg);
f_trg = f_trgo;
[f_trg, m, s] = zscore(f_trgo);
% mdl = fitrgp(x_trg,f_trg,'Standardize',1,'BasisFunction','none','KernelFunction','ardsquaredexponential','FitMethod','exact','PredictMethod','exact');
% mdl_kernel_parameters=mdl.KernelInformation.KernelParameters;
% mdl = fitrgp(x_trg,f_trg,'Standardize',1,'BasisFunction','none','KernelFunction','ardsquaredexponential','FitMethod','none','Sigma',1e-30,'PredictMethod','exact','Kernelparameters',mdl_kernel_parameters');
mdl = fitrgp(x_trg, f_trg, 'Standardize',1,'BasisFunction','none','KernelFunction','ardsquaredexponential','FitMethod','exact','Sigma',1e-30,'PredictMethod','exact');


x_tst = linspace(0, 1, 100);
f_tst = forrester(x_tst');

f_pre = predict(mdl, x_tst');
f_pre = denormzscore(f_trgo, f_pre);

plot(x_tst, f_tst); hold on;
plot(x_tst, f_pre); hold on;
plot(x_trg, f_trgo, 'ro'); hold on;

%plot(x_trg, f_trg,'ro');
end

function f=forrester(x_trg)
f= (6 .* x_trg - 2).^2 .* sin(12 .* x_trg -4);
return
end