function  lowersave(x, f, c, prob, seed, method,  varargin)

num = length(prob.xl_bl);
savepath = strcat(pwd, '\result_folder\', prob.name, '_', num2str(num) ,'_',method);
% savepath = strcat(pwd, '\result_folder\', prob.name, '_',method);
if ~isempty( varargin)
    init_size =  varargin{1};
    savepath = strcat(savepath, '_init_', num2str(init_size));
end


n = exist(savepath);
if n ~= 7
    mkdir(savepath)
end

savename_xu = strcat(savepath, '\xl_', num2str(seed),'.csv');
savename_fu = strcat(savepath, '\fl_', num2str(seed),'.csv');
savename_fc = strcat(savepath, '\cl_', num2str(seed),'.csv');

csvwrite(savename_xu, x);
csvwrite(savename_fu, f);
csvwrite(savename_fc, c);

end