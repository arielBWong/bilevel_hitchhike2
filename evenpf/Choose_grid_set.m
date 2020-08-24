function [finalsetnd,finalset] = choose_grid_set(ndset,numreq)
% Ensuring the first selected point is the extremity
[~,b]=min(ndset(:,1));
idall = [b setdiff(randperm(size(ndset,1)),b)];

if size(ndset,1) >= numreq
    finalset = idall(1);
    leftset = setdiff((1:size(ndset,1)),finalset);
    while numel(finalset) < numreq
        mindist = zeros(numel(leftset),2);
        for i = 1:numel(leftset)
            dist = [];
            dist = (sqrt(sum((repmat(ndset(leftset(i),:),numel(finalset),1) - ndset(finalset,:)).^2,2)))';
            mindist(i,:) = [leftset(i) nanmin(dist)];
        end
        [~,idmax] = max(mindist(:,2));
        finalset = [finalset mindist(idmax,1)];
        leftset = setdiff((1:size(ndset,1)),finalset);
    end
    finalsetnd = unique(ndset(finalset,:),'rows','stable');
else
    finalsetnd = unique(ndset,'rows','stable');
end
end