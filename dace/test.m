clear all; close all;clc;
addpath(genpath(strcat(pwd,strcat(filesep,'dace'))));
% xtrg=rand(100,4);
% f1trg=sum(xtrg.^2,2);
% f2trg=sum(xtrg,2);
% model=dace_train(xtrg,[f1trg f2trg]);
% xtst=rand(200,4);
% f1tst=sum(xtst.^2,2);
% f2tst=sum(xtst,2);
% [fprd,sigmaprd]=dace_predict(xtst,model);
% figure(1);plot(fprd(:,1),f1tst,'ro');
% figure(2);plot(fprd(:,2),f2tst,'ro');

xtrg=linspace(0,1,3);  
ytrg=((6*xtrg-2).^2).*sin(12*xtrg-4);
model=dace_train(xtrg',ytrg');
%model.
xtst=linspace(0,1,100);
ytst=((6*xtst-2).^2).*sin(12*xtst-4);
[yprd,sigmaprd]=dace_predict(xtst',model);
figure(1);
plot(xtst,yprd,'r-');hold on;
plot(xtst,ytst,'b-');
plot(xtrg,ytrg,'bo');
plot(xtst,yprd+sigmaprd,'k-');
plot(xtst,yprd-sigmaprd,'k-');

%% This is perfectly fine for unconstrained SO problems
fmin=min(ytrg);
for i=1:length(xtst)
    [mu,sigma]=dace_predict(xtst(i),model);
    muF=mu(:,1:size(fmin,1));sigmaF=sigma(:,1:size(fmin,1));
    u = (fmin - muF);
    z = u./sigmaF;
    f1a = u.*normcdf(z, 0, 1);
    f1a(sigmaF==0,1) = 0;
    f1b = sigmaF.*normpdf(z, 0, 1);
    f = -(f1a + f1b);
    EI(i)=f;
end
figure(2);plot(xtst,EI,'b-');
% Next sampling point should be xtst(b), where [~,b]=min(EI);

% xtrg,ytrg, Build model
% Solve an optimization problem which is find the location x such that EI
% is miximization
% Initialize a populaltion of solutions, use model to calculate the EI
% correspoinding to all solutions in the initilzed population, use EA to
% maximize EI
% For that solution, please do a real evaluation to get the real f value.
% Retrain the model and continue the process

%% Constrained problem
% Is there any feasible solution so far. xtrg and ftrg, g1trg,g2trg etc etc

% Yes there is atleast a single feasible solution, then it is simple: You need fmin as the obj of the feasible solution
% For every x, you can also compute a term called probability of
% feasibility.
% Please maximize EI*PF

% If no feasible solution so far: Maximize PF only





