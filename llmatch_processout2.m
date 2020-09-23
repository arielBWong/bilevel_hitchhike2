%% evaluate 100 median prediction and real function 
% conduct statistic analysis
% corr(,,'type', 'kendall')
%
% for a seed per problem
% (1) load xl, generate xu
% (2) iterate from  index 20 every 10 more until 60,
%       read xl, match xu, train dace, predict best solution
% (3) save mse to another matrix
% (4)

clearvars;
close all;



seedmax = 1;
% problems ={'dsm1(3, 3)', 'dsm1d(3, 3)','dsm1dc1(3, 3)','dsm1dc2(3, 3)', ...
%                          'dsm2(3, 3)', 'dsm2d(3, 3)','dsm2dc1(3, 3)','dsm2dc2(3, 3)',...
%                            'dsm3(3, 3)', 'dsm3d(3, 3)','dsm3dc1(3, 3)','dsm3dc2(3, 3)'};


%           problems = { 'dsm1(2,2)', 'dsm1d(2,2)','dsm1dc1(2,2)','dsm1dc2(2,2)',...
%          'dsm2(2,2)', 'dsm2d(2,2)','dsm2dc1(2,2)','dsm2dc2(2,2)',...
%               'dsm3(2,2)', 'dsm3d(2,2)','dsm3dc1(2,2)','dsm3dc2(2,2)' }; % change p3 term back to scale 10
          
%            problems = { 'dsm1(4,4)', 'dsm1d(4,4)','dsm1dc1(4,4)','dsm1dc2(4,4)',...
%          'dsm2(4,4)', 'dsm2d(4,4)','dsm2dc1(4,4)','dsm2dc2(4,4)',...
%               'dsm3(4,4)', 'dsm3d(4,4)','dsm3dc1(4,4)','dsm3dc2(4,4)' }; % change p3 term back to scale 10
          
% problems = { 'dsm1(5,5)', 'dsm1d(5,5)','dsm1dc1(5,5)' };        % ,'dsm1dc2(5,5)'
problems = { 'smd10(1,1,1)' };       
 %  problems = { 'dsm1(2, 2)', 'dsm1d(2, 2)','dsm1dc1(2, 2)' };  
 % problems = { 'dsm1(4, 4)', 'dsm1d(4, 4)','dsm1dc1(4, 4)' };  
                       
methods = {'llmatcheim',  'llmatchble'};  % 'llmatchpop',
leg = {'EIM', 'BEL'};
np= length(problems);
nm = length(methods);


% lower level has fixed number of true evaluation before local search
% 60 (20, 30, 40, 50, 60)
% last one saved for local search final results
collectmatrix = zeros(6, np*nm*2);

for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    
   upper_bound = prob.xl_bu;
   lower_bound = prob.xl_bl;
   init_size = 100;
   sample_xl = lhsdesign(init_size,prob.n_lvar,'criterion','maximin','iterations',1000);
   sample_xl = repmat(lower_bound, init_size, 1) ...
                + repmat((upper_bound - lower_bound), init_size, 1) .* sample_xl;
    
    nv = prob.n_lvar;
    
%     k = 2: prob.n_lvar;
%     k = (k-1)/2;
%     xu = [0, k];
xu=[1,1];
    
    fig1 = gcf;


    x = 1:6;
    barplot = [];
    barplotdev = [];
    for jj = 1:nm
        method = methods{jj};
        num = length(prob.xl_bl);
        savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_',method);
        collectionpermethod= zeros( 6, seedmax);
        for kk = 1:seedmax
            savename = strcat(savepath, '\xl_', num2str(kk), '.csv' )
            xl = csvread(savename);
            
            % every 10 training data make a prediction
            for mm = 1:5
                xl_sep = xl(1: (mm-1) * 10 + 20, :);
                mse = prediction(xu, xl_sep, prob);
                % mse = collection_prediction(sample_xl, xu, xl_sep, prob);
                
               %  kend =  kendall(sample_xl, xu, xl_sep, prob);
                collectionpermethod(mm, kk) = mse;
            end
            
            % put the final result of search to the end
            savename = strcat(savepath, '\fl_', num2str(kk), '.csv' );
            fl = csvread(savename);
            collectionpermethod(end, kk) = fl(end, 1);
        end
        
        % finish one method, fill to the biggest matrix
        % use mean and std first
        mean_permethod  = mean(collectionpermethod, 2);
        std_permethod = std(collectionpermethod, [], 2);
        barplot = [barplot, mean_permethod];
        barplotdev = [barplotdev, std_permethod];
        
        collectmatrix(:, (ii-1)*nm*2 + (jj-1) * 2 + 1) = mean_permethod;
        collectmatrix(:, (ii-1)*nm*2+ (jj-1) *2 + 2)  = std_permethod;
    end
    
  %  ---------plots--
 % for each problem plot three methods with error bar
    bar(x, barplot, 'FaceColor','flat');  hold on;
   
    hBar = bar(barplot, 0.8);
    colors = [[0.8500 0.3250 0.0980]; [0 0.4470 0.7410];  [0.9290 0.6940 0.1250] ];
    
    pause(2); % allow offset calculated
    for k1 = 1:size(barplot,2)
        ctr(k1,:) = hBar(k1).XData + hBar(k1).XOffset;
        ydt(k1,:) = hBar(k1).YData;
        hBar(k1).FaceColor  =colors(k1, :);
    end
    errorbar(ctr, ydt, barplotdev', '.k');

    xlabel('training data size', 'FontSize', 16);
    ylabel('mean mse on gobal optimal solution', 'FontSize', 16);
    t = strcat(prob.name, ' k ', num2str(nv)); 
    title(t,  'FontSize', 16);
    set(gca, 'XTickLabel',{'20','30','40','50','60','local search'}, 'FontSize', 12);
    legend(hBar, leg{1}, leg{2}, 'FontSize', 14); %methods{3}, leg{3}, 
    
    savename = strcat(pwd, '\result_folder\plots3\', prob.name, '_llmatchperformance.fig');
    savefig(savename);
    savename = strcat(pwd, '\result_folder\plots3\', prob.name, '_llmatchperformance.png');
    saveas(fig1, savename);
    close(fig1);
    %---------plots--
  
    
end

 % collection 
   % save into csv
savepath = strcat(pwd, '\result_folder\plots3\dsm_llres3.csv');
fp=fopen(savepath,'w');

fprintf(fp, 'problem_method, ');
fprintf(fp, '20, ');
fprintf(fp, '30, ');
fprintf(fp, '40, ');
fprintf(fp, '50, ');
fprintf(fp, '60, ');
fprintf(fp, '61\n ');

% construct index 
index = cell(1,  np * nm * 2);
nn = 1;
for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    for jj = 1:nm
        index{nn} = strcat(prob.name, '_', methods{jj}, '_mean');
        nn = nn + 1;
        index{nn} = strcat(prob.name, '_', methods{jj}, '_std');
        nn = nn + 1;
    end
end
% format header
collectmatrix = collectmatrix';
for  ii =1:  np * nm * 2    
    fprintf(fp, '%s, ', index{ii} );
    for jj = 1: 6
        fprintf(fp, '%f, ', collectmatrix(ii, jj));
    end
    fprintf(fp, '\n');
end
fclose(fp);

function mse = prediction(xu, xl, prob)
% re-train model to get prediction on global optimum

num_xl = size(xl, 1);
xu = repmat(xu, num_xl, 1);
[fl, fc] = prob.evaluate_l(xu, xl);
[krg_obj, krg_con, info] = update_surrogate(xl, fl, fc);

% first three problem lower value stays the same
if ~contains(prob.name, 'c2' )
    xl_prime = xu(1,:);
    xl_prime = [2, pi/4];
    [fl_prime, ~] = prob.evaluate_l(xu(1, :), xl_prime);
    
%     if abs(fl_prime) >0.001
%         error('lower f prime value is wrong');
%     end
%     
    % check prediction
    fpred = dace_predict(xl_prime, krg_obj{1});
    fpred = denormzscore( fl, fpred);
    
    mse = (fpred - fl_prime)^2;
else
    xl_prime=xu(1, :);
    % xl_prime(2) = -0.5;     % if the lower level equation term p3 is  multiplied by 1, then there are two lower level optimum exists
     %  xl_prime(2) = -1.5;     % if the lower level equation term p3 is  multiplied by 10, then there are two lower level optimum exists
    xl_prime(2) = -0.5;     % if the lower level equation term p3 is  multiplied by 10, then there are two lower level optimum exists
    [fl_prime, ~] = prob.evaluate_l(xu(1, :), xl_prime);
    
    
    % check prediction
    fpred = dace_predict(xl_prime, krg_obj{1});
    fpred = denormzscore( fl, fpred);
    
    mse = (fpred - fl_prime)^2;
    
end




end

function kendallstat = kendall(xlsample, xu, xl, prob)
% there are 100 points in xlsample
% these xlsamples are the same every time this method called
num_xl = size(xl, 1);
xu = repmat(xu, num_xl, 1);
[fl, fc] = prob.evaluate_l(xu, xl);
[krg_obj, krg_con, info] = update_surrogate(xl, fl, fc);

xu = xu(1, :);
num_sample = size(xlsample, 1);
xu =  repmat(xu, num_sample, 1);

% sample real value
[samplefl, ~] = prob.evaluate_l(xu, xlsample);

% predict the sampled data
samplepred = dace_predict(xlsample, krg_obj{1});
samplepred = denormzscore( fl, samplepred);

% get kendall stat
% samplefl = sort(samplefl);
% samplepred = sort(samplepred);

kendallstat = corr(samplefl, samplepred, 'type', 'Kendall');
end

function err = collection_prediction(xlsample, xu, xl, prob)
% there are 100 points in xlsample
% these xlsamples are the same every time this method called
num_xl = size(xl, 1);
xu = repmat(xu, num_xl, 1);
[fl, fc] = prob.evaluate_l(xu, xl);
[krg_obj, krg_con, info] = update_surrogate(xl, fl, fc);

xu = xu(1, :);
num_sample = size(xlsample, 1);
xu =  repmat(xu, num_sample, 1);

% sample real value
[samplefl, ~] = prob.evaluate_l(xu, xlsample);

% predict the sampled data
samplepred = dace_predict(xlsample, krg_obj{1});
samplepred = denormzscore( fl, samplepred);
err = immse(samplepred, samplefl);
end


function f = denormzscore(trainy, fnorm)
[train_y_norm, y_mean, y_std] = zscore(trainy, 1, 1);
f = fnorm * y_std + y_mean;
end

%-----auxiliary function ---
function [krg_obj, krg_con, info] = update_surrogate(trainx, trainy, trainc)
% this function updates train kriging model from train x and train y
%
%
if size(trainy,2)>1
    error(' following zscore norm only applies to single objective problems');
end
[train_y_norm, y_mean, y_std] = zscore(trainy, 1, 1);
num_obj = size(trainy, 2);
krg_obj = cell(1, num_obj);
num_vari = size(trainx, 2);
for ii = 1:num_obj
    % kriging_obj{ii} = dace_train(train_x_norm,train_y_norm(:,ii));
    krg_obj{ii} = dacefit(trainx,train_y_norm(:,ii),...
        'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  % for test
end

info = struct();
info.ymean = y_mean;
info.ystd = y_std;

% deal with constraints
if ~isempty(trainc)
    num_con = size(trainc, 2);
    krg_con = cell(1, num_con);
    [train_c_norm, c_mean, c_std] = zscore(trainc, 1, 1);
    
    for ii = 1:num_con
        krg_con{ii} = dacefit(trainx, train_c_norm(:,ii),...
            'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));  % for test
    end
    info.cmean = c_mean;
    info.cstd = c_std;
else
    krg_con = [];
    info.cmean = [];
    info.cstd = [];
end

end