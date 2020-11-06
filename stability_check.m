
function flag = stability_check(trainx, trainy, trainc, krg_obj, krg_con)
flag = true;
n = length(krg_obj);
for i = 1:n
    [y, sig] = dace_predict(trainx, krg_obj{i});
     y =  denormzscore(trainy, y);
     deltaf = abs(y - trainy)
     
     if max(deltaf) >1
         flag = false;
         return
     end    
end

if ~isempty(krg_con)
    nc = length(krg_con);
    for i = 1:nc
        [c, sig] = dace_predict(trainx, krg_con{i});       
        deltac = abs(c - trainc);
        
        if max(deltac) >1
            flag = false;
            return
        end
    end
end

end
