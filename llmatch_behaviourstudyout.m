%% lower level investigation
clearvars;
close all;

seedmax = 29;
problems = {'smd1()', 'smd2()','smd3()', 'smd4()',  'smd5()',   'smd6()', 'smd7()', 'smd8()',  'smd9()',   'smd10()', 'smd11()', 'smd12()',...
    'dsm1(2,2)','dsm1(3,3)', 'dsm1(4,4)','dsm1(5,5)','dsm1dc1(2,2)','dsm1dc1(3,3)', 'dsm1dc1(4,4)',  'dsm1dc1(5,5)'};
problems ={
     'tp3(5,5)'...
    };
% 'tp3(5,5)','tp5(5,5)','tp6(5,5)','tp7(5,5)','tp8(5,5)','tp9(5,5)'
% 'tp3(4,4)','tp5(4,4)','tp6(4,4)','tp7(4,4)','tp8(4,4)','tp9(4,4)'
% 'tp3(3,3)','tp5(3,3)','tp6(3,3)','tp7(3,3)','tp8(3,3)','tp9(3,3)'
% 'tp3(2,2)','tp5(2,2)','tp6(2,2)','tp7(2,2)','tp8(2,2)','tp9(2,2)'

k = 5;
init_size = 11 * k - 1;
methods = {'llmatchapt', 'llmatcheim',  'llmatchble'};  % 'llmatchpop', 'llmatcheim',  'llmatchble',  'llmatchapt'
% leg = {'EIM', 'BEL', 'GEN'};

% methods = {'llmatch', 'llmatch_sao_archiveinsert'};  % 'llmatchpop','llmatcheim', 'llmatcheimfix',,  'llmatch_hyb'
leg = { 'HYB','EIM', 'BLE'}; % , 'HYB'
np  = length(problems);
nm  = length(methods);
infill_size = 300;

num_experiments = 100; % size xu in its file

% ------read in file
out = cell(1, np);



for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    out{ii} = cell(1, nm);
    
    for jj = 1:nm
        method = methods{jj};
        num = length(prob.xl_bl);
        %savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_llstudy_',method);
        savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_',method);
        savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_',method, '_init_', num2str(init_size));
        
        for kk = 1:seedmax
            % read algorithm fl
            savename = strcat(savepath, '\fl_', num2str(kk), '.csv' )
            fl = csvread(savename);
            
            % create prime fl
            fl_prime = [];
            % savename = strcat(savepath, '\xu_', num2str(kk), '.csv');
            % xu = csvread(savename);
            
            % xu1 = xu(1,:);
            % [match_xl, n_fev, flag] = match_method (xu1, prob,20, num_gen, freq, seed);
            
            savename = strcat(savepath, '\xl_', num2str(kk), '.csv');
            xl_search = csvread(savename);
            xl_prime = [];
            %         for mm = 1:size(xu, 1)
            %             if contains(prob.name, 'SMD')
            %                 xl = prob.get_xlprime(xu(mm, :));
            %             else
            %                 % xl = xu(mm, :);
            %                 xl = prob.xl_prime;
            %             end
            %             fl_p = prob.evaluate_l(xu(mm, :), xl);
            %             fl_prime = [fl_prime; fl_p];
            %         end
            % deltaF = abs(fl - fl_prime);
            temp =  fl(init_size:end-1)';
            for mm = 1:size(temp, 2)
                out{ii}{jj}(kk, mm) = min(temp(1:mm));
            end
        end
        % out{ii}(jj, :) = deltaF;
    end
    
    
end
xbase = init_size : init_size + infill_size;


% medianmatrix = zeros(np, nm);
% stdmatrix = zeros(np, nm);
% meanmatrix = zeros(np, nm);
% for ii = 1:np
%     for jj = 1:nm
%         medianmatrix(ii, jj) = median(out{ii}(jj, :));
%         meanmatrix(ii, jj) = mean(out{ii}(jj, :));
%         stdmatrix(ii, jj) = std(out{ii}(jj, :));
%     end
% end

meanfl = zeros(nm, infill_size + 1);
stdfl = zeros(nm, infill_size + 1);
for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    for jj = 1:nm
        meanfl(jj, :) = mean(out{ii}{jj}, 1);
        stdfl(jj, :) = std(out{ii}{jj}, 1);
    end
    
    figure(ii);
    xbase = init_size : init_size + infill_size;
    x = [xbase, fliplr(xbase)];
    plot(xbase, meanfl(1, :), 'r', 'LineWidth',2); hold on;
    plot(xbase, meanfl(2, :), 'k', 'LineWidth',2); hold on;
    plot(xbase, meanfl(3, :), 'b', 'LineWidth',2); hold on;
    
    y1 = meanfl(1, :) + stdfl(1,:);
    y2 = meanfl(1, :) - stdfl(1,:);
    y = [y1, fliplr(y2)];
    fill(x, y, 'r', 'FaceAlpha', 0.1, 'EdgeColor','none');
    
    y1 = meanfl(2, :) + stdfl(2,:);
    y2 = meanfl(2, :) - stdfl(2,:);
    y = [y1, fliplr(y2)];
    fill(x, y, 'k', 'FaceAlpha', 0.4, 'EdgeColor','none');
    
    y1 = meanfl(3, :) + stdfl(3,:);
    y2 = meanfl(3, :) - stdfl(3,:);
    y = [y1, fliplr(y2)];
    fill(x, y, 'b', 'FaceAlpha', 0.3, 'EdgeColor','none');
    
    numk = prob.n_lvar;
    legend( leg{1}, leg{2},leg{3}, 'FontSize', 14);
    t = strcat( prob.name , '  converge k', num2str(numk), ' init size ', num2str(init_size));
    title(t,  'FontSize', 16);
    a = 0;
    

    savename = strcat(pwd, '\result_folder\', prob.name,'_', num2str(numk), '_converge.fig');
    savefig(savename);
    savename = strcat(pwd, '\result_folder\', prob.name,'_', num2str(numk), '_converge.png');
    saveas(figure(ii), savename);
    close(figure(ii));
end


%
% % to csv
% %-----mean to csv
% foldername = strcat(pwd, '\result_folder\matchxl_invest');
% n = exist(foldername);
% if n ~= 7
%     mkdir(foldername)
% end
%
%
%
% savename = strcat(foldername,'\matchxl2_mean.csv');
% fp=fopen(savename,'w');
% fprintf(fp, 'problem_method, ');
% for jj = 1:nm
%     fprintf(fp,'%s, ', leg{jj} );
% end
% fprintf(fp,'\n');
%
% for ii = 1:np
%     prob = problems{ii};
%     prob = eval(prob);
%     name = strcat( prob.name, '_', num2str(prob.n_lvar));
%     fprintf(fp, '%s,', name);
%     for jj = 1:nm
%
%         st = strcat(num2str(meanmatrix(ii, jj)), char(177) , num2str(stdmatrix(ii, jj)));
%         fprintf(fp, '%s,', st);
%     end
%     fprintf(fp, '\n');
% end
% fclose(fp);
%
%
%
% savename = strcat(foldername,'\matchxl2_std.csv');
% fp=fopen(savename,'w');
% fprintf(fp, 'problem_method, ');
% for jj = 1:nm
%     fprintf(fp,'%s, ', leg{jj} );
% end
% fprintf(fp,'\n');
%
% for ii = 1:np
%     prob = problems{ii};
%     prob = eval(prob);
%     name = strcat( prob.name, '_', num2str(prob.n_lvar));
%     fprintf(fp, '%s,', name);
%     for jj = 1:nm
%
%         st = num2str(stdmatrix(ii, jj));
%         fprintf(fp, '%s,', st);
%     end
%     fprintf(fp, '\n');
% end
% fclose(fp);
%
%
% % to csv
% %-----median to csv
% foldername = strcat(pwd, '\result_folder\matchxl_invest');
% n = exist(foldername);
% if n ~= 7
%     mkdir(foldername)
% end
%
% savename = strcat(foldername,'\matchxl2_median.csv');
% fp=fopen(savename,'w');
% fprintf(fp, 'problem_method, ');
% for jj = 1:nm
%     fprintf(fp,'%s, ', leg{jj} );
% end
% fprintf(fp,'\n');
%
% for ii = 1:np
%     prob = problems{ii};
%     prob = eval(prob);
%     name = strcat( prob.name, '_', num2str(prob.n_lvar));
%     fprintf(fp, '%s,', name);
%     for jj = 1:nm
%         fprintf(fp, '%f,', medianmatrix(ii, jj));
%     end
%     fprintf(fp, '\n');
% end
% fclose(fp);
%

















function mse = prediction(xu, xl, prob)
% re-train model to get prediction on global optimum
% xl is for training
% xl_prime is the real global optimal
num_xl = size(xl, 1);
xu = repmat(xu, num_xl, 1);
[fl, fc] = prob.evaluate_l(xu, xl);
[krg_obj, krg_con, info] = update_surrogate(xl, fl, fc);
% first three problem lower value stays the same
if ~contains(prob.name, 'c2' )
    if contains(prob.name, 'SMD')
        xl_prime = prob.xl_prime;
    else
        xl_prime = xu(1,:);
        
    end
    [fl_prime, ~] = prob.evaluate_l(xu(1, :), xl_prime);
    
    if abs(fl_prime) >0.001
        error('lower f prime value is wrong');
    end
    
    % check prediction
    fpred = dace_predict(xl_prime, krg_obj{1});
    fpred = denormzscore( fl, fpred);
    
    mse = abs(fpred - fl_prime);
else
    xl_prime=xu(1, :);
    % xl_prime(2) = -0.5;     % if the lower level equation term p3 is  multiplied by 1, then there are two lower level optimum exists
    %  xl_prime(2) = -1.5;     % if the lower level equation term p3 is  multiplied by 10, then there are two lower level optimum exists
    xl_prime(2) = -0.5;     % if the lower level equation term p3 is  multiplied by 10, then there are two lower level optimum exists
    [fl_prime, ~] = prob.evaluate_l(xu(1, :), xl_prime);
    
    
    % check prediction
    [fpred, ~] = dace_predict(xl_prime, krg_obj{1});
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


function dsm_rebuild(xu, xl_all, prob)
% use post process to rebuild plot

if prob.n_lvar > 2
    error('process not compatible with more than 2 variables');
end
xl = xl_all(1:end-1, :); % last one is returned matching xl

num_xl = size(xl, 1);
xu = repmat(xu, num_xl, 1);
[fl, fc] = prob.evaluate_l(xu, xl);
[krg_obj, krg_con, info] = update_surrogate(xl, fl, fc);


funh_obj = @(x)llobj(x, krg_obj);
funh_con = @(x)llcon(x, krg_con);

param.gen       = 100;
param.popsize   = 100;
[~,~,~, archive] = gsolver(funh_obj, prob.n_lvar,  prob.xl_bl, prob.xl_bu, [], funh_con, param);
bestx = archive.pop_last.X(1, :);




fig1 = gcf;

k =2;
xl1 = linspace(-k, k, 100);
xl2 = linspace(-k, k, 100);
[xl1, xl2] = meshgrid(xl1, xl2);
f = zeros(100, 100);
%  for i = 1:100
%      for j = 1:100
%          f (i, j) = (xl1(i, j) - 0 ).^2 + (xl2(i, j) - 0.5).^2 + 10 * abs(sin(pi/2*(xl2(i, j) - 0.5)));
%      end
%  end

for i = 1:100
    for j = 1:100
        [f(i, j), ~] = dace_predict([xl1(i, j), xl2(i, j)], krg_obj{1});
        f(i, j) = denormzscore( fl, f(i, j));
    end
end
min(f(:))
surf(xl1, xl2, f); hold on;
xlabel('xl1', 'FontSize', 16);
ylabel('xl2', 'FontSize', 16);
zlabel('fl',  'FontSize', 16);
colormap jet
shading interp

scatter3(xl(:, 1), xl(:, 2), fl(:, 1), 80,  'yx')
[ret, ~] = prob.evaluate_l(xu(1, :), xl_all(end, :));
scatter3(xl_all(end,1), xl_all(end, 2),ret, 80, 'ro', 'filled' );

[surr_search, ~] = prob.evaluate_l(xu(1, :), bestx);
scatter3(bestx(1), bestx(2),surr_search, 80, 'yo', 'filled' );
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


% believer objectives
function  f = llobj(x, kriging_obj)
num_obj = length(kriging_obj);   % krg cell array?
num_x = size(x, 1);
f = zeros(num_x, num_obj);
for ii =1:num_obj
    [f(:, ii), ~] = dace_predict(x, kriging_obj{ii});
end
end


% believer constraints
function c = llcon(x, krging_con)
num_con = length(krging_con);
num_x = size(x, 1);
c = zeros(num_x, num_con);
for ii =1:num_con
    [c(:, ii), ~] = dace_predict(x, krging_con{ii});
end
end




