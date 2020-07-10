function [who,which_obj]=Select(mu,sigma,prob)
f=mu(1:size(mu,1)/2,:);
[fronts,~,~] = nd_sort(f,[1:size(f,1)]');
nd_ids=fronts(1).f;
for i=1:size(mu,1)/2
    for j=1:length(nd_ids)
        p_obj(i,j,1:size(mu,2))= Compute_prob_dom_obj(mu(size(mu,1)/2+i,:),mu(nd_ids(j),:),sigma(size(mu,1)/2+i,:),sigma(nd_ids(j),:),prob);
    end
end

who=[];which_obj=[];
for i=1:length(nd_ids)
    for k=1:size(mu,2)
        t=p_obj(:,i,k);
        l=reshape(t,size(mu,1)/2,1);
        % Among these can only select the best that has so far not been
        % evaluated
        tmp=sigma(size(mu,1)/2+1:size(mu,1),k)~=0;
        if(sum(tmp~=0))
            [~,b]=max(l.*tmp);
            who=[who b];
            which_obj=[which_obj k];
        end
    end
end
who=who+size(mu,1)/2;
return



