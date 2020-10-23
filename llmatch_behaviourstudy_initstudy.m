%%
clearvars;
close all;
seedmax = 29;

k = 3;
init_size = 5* k + 1;

problems ={
    'smd1()'...
    };
% 'tp3(5,5)','tp5(5,5)','tp6(5,5)','tp7(5,5)','tp8(5,5)','tp9(5,5)'
% 'tp3(4,4)','tp5(4,4)','tp6(4,4)','tp7(4,4)','tp8(4,4)','tp9(4,4)'
% 'tp3(3,3)','tp5(3,3)','tp6(3,3)','tp7(3,3)','tp8(3,3)','tp9(3,3)'
% 'tp3(2,2)','tp5(2,2)','tp6(2,2)','tp7(2,2)','tp8(2,2)','tp9(2,2)'




methods = {'llmatcheim',  'llmatchble', 'llmatchapt'};  % 'llmatchpop',
leg = { 'EIM','BLE', 'HYB'};
np  = length(problems);
nm  = length(methods);

max_iter = 300;


plotmatrix_2eim = zeros(np, max_iter  + 1); % 10 beats eim,  0 no significant
plotmatrix_2ble = zeros(np, max_iter  + 1); % -10 beats ble, 0 no significant

for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    num = prob.n_lvar;
    name  = prob.name;
    
    algs_results = cell(1, nm);
    
    for jj = 1:nm
        algs_results{jj} = zeros(seedmax, max_iter+ 1);
        savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num), '_', methods{jj},  '_init_', num2str(init_size));
        for kk = 1: seedmax
            filename  = strcat(savepath, '\fl_', num2str(kk), '.csv')
            fl = csvread(filename);
            tmp = fl(init_size: init_size + max_iter)';
            for nn = 1 : max_iter+1
                algs_results{jj}(kk, nn) =  min(tmp(1:nn));
            end
            
        end
    end
    
   
    
    % ---------
    
%     eim_median = zeros(1, max_iter + 1);
%     ble_median = zeros(1, max_iter + 1);
%     hyb_median = zeros(1, max_iter + 1);
%     for nn = 1:max_iter+1
%         eim = algs_results{1}(:, nn);
%         ble = algs_results{2}(:, nn);
%         hyb = algs_results{3}(:, nn);
%         
%         eim_median(nn) = median(eim);
%         ble_median(nn) = median(ble);
%         hyb_median(nn) = median(hyb);
%     end
    
%     leg = { 'EIM','BLE', 'HYB'}; % , 'HYB'
% 
%     plot(eim_median(1:100), 'r', 'LineWidth', 2); hold on;
%     plot(ble_median(1:100), 'b', 'LineWidth', 2); hold on;
%     plot(hyb_median(1:100), 'k', 'LineWidth', 2); hold on;
%     legend( leg{1}, leg{2},leg{3}, 'FontSize', 14);
    
    
    % run significant test on each iteration
    pp = 1;
    for nn = 1:max_iter+1
        eim = algs_results{1}(:, nn);
        ble = algs_results{2}(:, nn);
        hyb = algs_results{3}(:, nn);
        
        %         % test whether equal median in median ble test
        %         [p,h] = ranksum(eim, ble);
        %         if h == 0
        %             plotmatrix(ii, pp) = 0;
        %         else
        %             % test increase in eim test
        %             [p,h] = ranksum(ble, eim,'tail','left');
        %             if h == 0
        %                 plotmatrix(ii, pp) = -10;
        %             else
        %                 % test increase in ble
        %                 [p,h] = ranksum(eim, ble,'tail','left');
        %                 if h == 0
        %                     plotmatrix(ii, pp) = 10;
        %                 else
        %                     plotmatrix(ii, pp) = 0;
        %                 end
        %             end
        %         end
        
        % test increase in eim test
        [p,h] = ranksum(ble, eim);
        if h == 0
            plotmatrix_2eim(ii, pp) = 1;
        else
            % test increase in eim test
            [p,h] = ranksum(ble, eim, 'tail','left');
            if h == 0
                plotmatrix_2eim(ii, pp) = 10;
            else
                plotmatrix_2eim(ii, pp) = 5;
            end
        end
        % test increase in ble test
        [p,h] = ranksum(eim, ble);
        if h == 0
            plotmatrix_2eim(ii, pp) = -1;
        else
            % test increase in eim test
            [p,h] = ranksum(eim, ble,'tail','left');
            if h == 0
                plotmatrix_2ble(ii, pp) = -10;
            else
                plotmatrix_2ble(ii, pp) = -5;
            end
        end
        
        pp = pp + 1;
    end
end
plotmatrix = [plotmatrix_2eim; plotmatrix_2ble];
% plot
figure(1);
x = init_size: max_iter;

clims = [-10 10];
imagesc(plotmatrix(:, 1:pp-1), clims);
% Initialize a color map array of 256 colors.
colorMap = jet();
% Apply the colormap and show the colorbar
colormap(colorMap);

grid on;
colorbar;
yticks = 1:2;
yticklabels = {'compare2eim', 'compare2hyb'};
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
t = strcat(name, ' k=',num2str(k), ' init ', num2str(init_size), ' 300 iteration', ...
    ' red: beat eim, blue: beat ble');
title(t);
xlabel('iteration');

