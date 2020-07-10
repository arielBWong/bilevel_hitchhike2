function [mu,sigma]=Compute_mu_sigma(pop,prob,param,archive)
mu=zeros(size(pop.X,1),prob.nf+prob.ng);
sigma=zeros(size(pop.X,1),prob.nf+prob.ng);

for i=1:size(pop.X,1)
    tmp=pop.F(i,:);
    
    % These f's need to be predicted
    to_fill=find(isnan(tmp));
    already_filled=setdiff(1:prob.nf,to_fill);
    muf(i,already_filled)=tmp(already_filled);
    sigmaf(i,already_filled)=0;aug_id=[];
    
    if ~isempty(to_fill)
        % These are f's to be augmented with
        aug_id=setdiff(1:prob.nf,to_fill);
        
        % Filling the f's for this solution
        for k=1:length(to_fill)
            req_ids=[to_fill(k) aug_id];
            % Extracting all solutions which have the above complete information
            req_cols=repmat(1+prob.nx,1,length(req_ids))+req_ids;
            aug_cols=repmat(1+prob.nx,1,length(aug_id))+aug_id;
            ids=find(sum(~isnan(archive.sols(:,[req_cols aug_cols])),2)==length(req_cols)+length(aug_cols));
            
            %% Nearest Neighbour Search by normalized dimensions
            x_arc= unique(archive.sols(ids,2:prob.nx+1),'rows','stable');
            x_arc_n = Norm_dim(x_arc,prob.range);
            X_n = Norm_dim(pop.X(i,:),prob.range);
            arc_n = Norm_dim(archive.sols(:,2:prob.nx+1),prob.range);
            
            NNids = knnsearch(x_arc_n,X_n,'K',prob.nx*param.neighbour);
            
            x_NN_n= x_arc_n(NNids,:);
            [~,arc_ids] = ismember(x_NN_n,arc_n,'rows');
            
            X=[archive.sols(arc_ids,2:prob.nx+1) archive.sols(arc_ids,aug_cols)];
            F=archive.sols(arc_ids,to_fill(k)+prob.nx+1);
            
            %% Fitting the model for f
            modelname = 'krig';
            model = feval([modelname,'_train'],X,F);
            [muf(i,to_fill(k)),sigmaf(i,to_fill(k))] = feval([modelname,'_predict'],[pop.X(i,:) pop.F(i,aug_id)],model);
        end
    end
end

mu =muf; sigma = sigmaf;
return
