function perf_record(prob, fu, cu, fl, cl, n_up, n_low, seed, method)
%
%
%-create save folder
savepath = strcat(pwd, '\result_folder\', prob.name, '_', method);
n = exist(savepath);
if n ~= 7
    mkdir(savepath)
end
%-create save matrix
savematrix = zeros(3,3);
savematrix(1, 1) = abs(fu - prob.uopt);
savematrix(1, 2) = abs(fl - prob.lopt);
savematrix(2, 1) = n_up;
savematrix(2, 2) = n_low;

if ~isempty(cu)
    savematrix(3, 1) = cu;
else
    savematrix(3, 1) = -1;
    
end
if ~isempty(cl)
    savematrix(3, 2) = cl;
else
    savematrix(3, 2) = -1;
end
%-create save file name
savename = strcat(savepath, '\out_', num2str(seed),'.csv');
%-save and done
csvwrite(savename, savematrix);
end

