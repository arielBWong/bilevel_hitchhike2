function updated_order=Sparse_selection(id_fronts,f_fronts,guarantee)
all_ids=cell2mat(id_fronts);
id_frontone=all_ids(1:length(cell2mat(id_fronts(1))));
ndset=cell2mat(f_fronts(1));
if(size(ndset,1)>guarantee)
    [~,order] = Choose_grid_set(ndset,guarantee);
    remaining=setdiff([1:length(id_frontone)],order);
    if ~isempty(remaining)
        updated_order=[all_ids(order);all_ids(remaining) ;all_ids(guarantee+1:end)];
    else
        updated_order=[all_ids(order);all_ids(numreq+1:end)];
    end
else
    updated_order=all_ids;
end
return
