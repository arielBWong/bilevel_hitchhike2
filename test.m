function test
clear all;close all;clc;
x_trg= [0, 0.5, 1]';no_trials=10;no_samples=20;
f_trg=forrester(x_trg);
for i=1:no_trials
    mdl = fitrgp(x_trg,f_trg,'Standardize',1,'BasisFunction','none','KernelFunction','ardsquaredexponential','FitMethod','exact','PredictMethod','exact');
    mdl_kernel_parameters=mdl.KernelInformation.KernelParameters;
    mdl = fitrgp(x_trg,f_trg,'Standardize',1,'BasisFunction','none','KernelFunction','ardsquaredexponential','FitMethod','none','Sigma',1e-30,'PredictMethod','exact','Kernelparameters',mdl_kernel_parameters');
    
    % Do EI maximization using local search multistart(no_trials)
    
    x_trg=[x_trg;x_new];
    f_new=forrester(x_new);
    f_trg=[f_trg;f_new];
end
mdl = fitrgp(x_trg,f_trg,'Standardize',1,'BasisFunction','none','KernelFunction','ardsquaredexponential','FitMethod','exact','PredictMethod','exact');
mdl_kernel_parameters=mdl.KernelInformation.KernelParameters;
mdl = fitrgp(x_trg,f_trg,'Standardize',1,'BasisFunction','none','KernelFunction','ardsquaredexponential','FitMethod','none','Sigma',1e-30,'PredictMethod','exact','Kernelparameters',mdl_kernel_parameters');


xx=linspace(0,1,100);
ff=forrester(xx);
plot(xx,ff,'b-');hold on;
ffp=predict(mdl,xx');
plot(xx,ffp,'r--');
plot(x_trg,f_trg,'bo');
fp_trg=predict(mdl,x_trg);
disp(fp_trg-f_trg);
return

function f=forrester(x_trg)
f= (6 .* x_trg - 2).^2 .* sin(12 .* x_trg -4);
return


