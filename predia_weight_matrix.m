function [weights, ESS, sumSqrWeights,ttime] = predia_weight_matrix(ctrl, prior_data,obs_data, obs_err_std)

% CORE EVALUATON OF THE WEIGHTING MATRIX

% INPUT:            NAME                                        DIMENSION
% ===================================================================================
% ctrl              Control structure containing various info   STUCTURE
%     .no_err_marg  flag if marginalizing over obs errro        1
%                   (see Leube et al. 2012, WRR)   
% prior_data        data sample                                 DIM:N_MC
% obs_data          sample of observations                      DIM:N_MEAS
% obs_err_std       standart deviation of measurement error     DIM:1
%
% OUTPUT:           NAME                                        DIMENSION
% ===================================================================================
% weights           weighting matrix that contain the           N_MEAS:N_MC
%                   weights proportinal to the likelyhood of 
%                   each sample to represent the given obs.
%                   given the observation data sample
% ESS               Effective sample size for each condition    1:N_MEAS
% sumSqrWeights     Sum of squared weights for postprocessing   N_MEAS:1
% ttime             Time vector or init, main calculatio
%                   and post processing                         1:3
%
%                   sample (given data)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%
ttime(1)= cputime;

if strct_flag_check(ctrl,'no_err_marg')
    marg_factor = 1; % no additional smoothing of the likelood function
else
    marg_factor = 2; % additional smoothing of the likelood function by factor of 2
end

% determination of dimension sizes
[n_dim_data, n_mc  ] = size(prior_data);
[n_dim_obs  ,n_meas] = size(obs_data);

if strct_flag_check(ctrl,'n_mc')
   n_mc = ctrl.n_mc;
end

if n_dim_data ~= n_dim_obs
    error('Dimension of data disagree')
else
    n_dim = n_dim_data;
    clear n_dim_data n_dim_obs;
end

weights     = zeros(n_meas,n_mc);

ttime(1) = cputime -ttime(1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttime(2)= cputime;

% Normalization of data with observation errors
for i = 1:n_dim
    prior_data(i,:) = prior_data(i,:)  ./ obs_err_std(i) ./ sqrt(2*marg_factor);
    obs_data (i,:)  = obs_data(i,:)    ./ obs_err_std(i) ./ sqrt(2*marg_factor);
end

if strct_flag_check(ctrl,'no_err_marg')
    % adding measurement error in the simulated values in case that not
    % marinalized analytically
    prior_data = prior_data + randn(size(prior_data));
    
end

for j=1:n_mc
    for i=1:n_dim
        weights(:,j) = weights(:,j) + ((prior_data(i,j)-obs_data(i,:)).^2)';
    end
end

% Singel call of the exponent
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
    weights(del_idx,:) = NaN;
    if ~strct_flag_check(ctrl,'no_warning')
        warning(['### Warning: ' num2str(n_del) ' measurement realizations deleted ###'])
    end
end
idx_ = true(size(sumSqrWeights));
idx_(del_idx) = 0;
AESS = mean(1./sumSqrWeights(idx_));
ESS  = 1./sumSqrWeights';

ttime(3) = cputime - ttime(3);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
