function Create_plots(prob,pop,archive,evals,param,track)
if(param.strategy==2)
    str='NSGA-II';
else
    str='MF';
end
F=archive.sols(:,prob.nx+2:prob.nx+prob.nf+1);
id=find(sum(isnan(F),2)==0);
fnd=F(id,:);
[fronts,~,~] = nd_sort(fnd,[1:size(fnd,1)]');
ids=fronts(1).f;
figure(1);
if(prob.nf==2)
    plot(fnd(ids,1),fnd(ids,2),'ro');
    set(gca,'fontsize',16);box on; grid on;
    xlabel('Objective 1'); ylabel('Objective 2');
    saveas(gcf,strcat(pwd,filesep,'Results',filesep,num2str(param.prob_name),filesep,str,filesep,'Fig_1_',num2str(param.prob_name),'.fig'));
    saveas(gcf,strcat(pwd,filesep,'Results',filesep,num2str(param.prob_name),filesep,str,filesep,'Fig_1_',num2str(param.prob_name),'.eps'),'epsc');
end

if(prob.nf==3)
    plot3(fnd(ids,1),fnd(ids,2),fnd(ids,3),'ro');
    set(gca,'fontsize',16);box on; grid on;
    xlabel('Objective 1'); ylabel('Objective 2');zlabel('Objective 3');
    saveas(gcf,strcat(pwd,filesep,'Results',filesep,num2str(param.prob_name),filesep,str,filesep,'Fig_1_',num2str(param.prob_name),'.fig'));
    saveas(gcf,strcat(pwd,filesep,'Results',filesep,num2str(param.prob_name),filesep,str,filesep,'Fig_1_',num2str(param.prob_name),'.eps'),'epsc');
end

%% Figure 2 Plotting (Stacked bar plot of number of objective evals in each gen)
figure(2);
bar(evals,'stacked');
xlabel('Generation Number'); ylabel('Objective Evaluations');
set(gca,'fontsize',16);box on; grid on;
h= legend('Objective 1','Objective 2','Location','Best');
set(h,'fontsize',10);
saveas(gcf,strcat(pwd,filesep,'Results',filesep,num2str(param.prob_name),filesep,str,filesep,'Fig_2_',num2str(param.prob_name),'.fig'));
saveas(gcf,strcat(pwd,filesep,'Results',filesep,num2str(param.prob_name),filesep,str,filesep,'Fig_2_',num2str(param.prob_name),'.eps'),'epsc');

return