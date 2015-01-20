function [weights, AESS, sumSqrWeights,ttime,ESS] = predia_weight_matrix(ctrl, prior_data,obs_data, obs_err_std)

% version 1 / Jan 15 / AGeiges WNowak

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttime(1)= cputime;

if ~strct_flag_check(ctrl,'error_marginalizatons')
    marg_factor = 1;
else
    marg_factor = 2;
end

[n_dim_data, n_mc  ] = size(prior_data);
[n_dim_obs  ,n_meas] = size(obs_data);

if n_dim_data ~= n_dim_obs
    error('Dimension of data disagree')
else
    n_dim = n_dim_data;
    clear n_dim_data n_dim_obs;
end


%% Normalization of data
for i = 1:n_dim
    prior_data(i,:) = prior_data(i,:)  ./ obs_err_std(i) ./ sqrt(2*marg_factor);
    obs_data (i,:)  = obs_data(i,:)    ./ obs_err_std(i) ./ sqrt(2*marg_factor);
end

weights     = zeros(n_meas,n_mc);

ttime(1) = cputime -ttime(1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttime(2)= cputime;

%% ALTERNATIVE CALCULATION 1 (normally fasters)
for j=1:n_mc
    for i=1:n_dim
        weights(:,j) = weights(:,j) + ((prior_data(i,j)-obs_data(i,:)).^2)';
    end
end

% %% ALTERNATIVE CALCULATION 2
% for i=1:n_dim
%     weights = weights + ((repmat(prior_data(i,1:n_mc),n_meas,1) - repmat(obs_data(i,:)',1,n_mc)).^2);
% end
% 
% %% ALTERNATIVE CALCULATION 3
% for j=1:n_dim
%     [XMKX1 XMKX2]   = ndgrid(obs_data(j,:),prior_data(j,:));
%     weights         = weights + reshape((XMKX1(:)-XMKX2(:)).^2,n_meas,[]);
% end

%% Singel call of the exponent
weights   = exp(-weights);


ttime(2) = cputime - ttime(2);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% POSTPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttime(3)= cputime;
%%
if isfield(ctrl,'bayes_evidence')
    
    bayes_evidence = weights;
    
end

if isfield(ctrl,'diag_zero') && ctrl.diag_zero
    
    weights(sub2ind(size(weights),1:size(obs_data,2),1:n_meas)) = 0;
    
end

%% Normalizing of weights

sumWeights = sum(weights,2);
sumWeights(sumWeights< 2 * eps) = 1;
weights        = spdiags(1./sumWeights,0,n_meas,n_meas)*weights;

%% avoiding -inf errror when one weight is equal to one
sumSqrWeights = sum(weights.^2,2);
del_idx = sumSqrWeights> 0.999;
n_del = sum(del_idx);
if n_del > 0
    sumSqrWeights(del_idx) = NaN;
    weights(del_idx,:) = [];
    disp(['### Warning: ' num2str(n_del) ' Realizations deleted ###'])
end
idx_ = true(size(sumSqrWeights));
idx_(del_idx) = 0;
AESS = mean(1./sumSqrWeights(idx_));
ESS =1./sumSqrWeights;

ttime(3) = cputime - ttime(3);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
