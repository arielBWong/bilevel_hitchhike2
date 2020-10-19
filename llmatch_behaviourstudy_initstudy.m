%%
clearvars;
close all;
seedmax = 29;
problems ={
        'tp3(10, 10)','tp5(10, 10)','tp6(10, 10)','tp7(10, 10)','tp8(10, 10)','tp9(10, 10)'...
    };
% 'tp3(5,5)','tp5(5,5)','tp6(5,5)','tp7(5,5)','tp8(5,5)','tp9(5,5)'
% 'tp3(4,4)','tp5(4,4)','tp6(4,4)','tp7(4,4)','tp8(4,4)','tp9(4,4)'
% 'tp3(3,3)','tp5(3,3)','tp6(3,3)','tp7(3,3)','tp8(3,3)','tp9(3,3)'
% 'tp3(2,2)','tp5(2,2)','tp6(2,2)','tp7(2,2)','tp8(2,2)','tp9(2,2)'




methods = {'llmatcheim',  'llmatchble'};  % 'llmatchpop',
leg = { 'EIM','BEL'};
np  = length(problems);
nm  = length(methods);

max_iter = 300;
init_size = 11;

plotmatrix = zeros(np, max_iter  + 1); % 10 ble wins, -10 eim wins, 0 no significant

for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    num = prob.n_lvar;
    
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
    
    % run significant test on each iteration
    pp = 1;
    for nn = 1:max_iter+1
        eim = algs_results{1}(:, nn);
        ble = algs_results{2}(:, nn);
        
        % test whether equal median in median ble test
        [p,h] = ranksum(eim, ble);
        if h == 0
            plotmatrix(ii, pp) = 0;
        else
            % test increase in eim test
            [p,h] = ranksum(ble, eim,'tail','left');
            if h == 0
                plotmatrix(ii, pp) = -10;
            else
                % test increase in ble
                [p,h] = ranksum(eim, ble,'tail','left');
                if h == 0
                    plotmatrix(ii, pp) = 10;
                else
                    plotmatrix(ii, pp) = 0;
                end
            end
        end
        pp = pp + 1;
    end
end

% plot
figure(1);
x = init_size: max_iter;

clims = [-10 10];
imagesc(plotmatrix(:, 1:pp-1));
% Initialize a color map array of 256 colors.
colorMap = jet();
% Apply the colormap and show the colorbar
colormap(colorMap);
grid on;
colorbar;
title('k=5 problem significant test, 11 init, 300 iterations');

