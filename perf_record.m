function perf_record(prob, fu, cu, fl, cl, n_up, n_low, seed)
%                                   
%
%-create save folder
savepath = strcat(pwd, '\result_folder\', prob.name );
n = exist(savepath);
if n ~= 7
mkdir(savepath)
end
%-create save matrix
savematrix = zeros(3,3);
savematrix(1,1) = abs(fu - prob.uopt);
savematrix(1, 2) = abs(fl - prob.lopt);
savematrix(2, 1) = n_up;
savematrix(2, 2) = n_low;
savematrix(3, 1) = cu;
savematrix(3, 2) = cl;
%-create save file name
savename = strcat(savepath, '\acc_con_fea_', num2str(seed),'.csv');
%-save and done
csvwrite(savename, savematrix);
end

