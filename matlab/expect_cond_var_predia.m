function [E_cond_var,ESS, cond_var] = expect_cond_var_predia(ctrl, prior_data,obs_data, obs_err_std,pred_data,prior_weight)
% Author: Andreas Geiges
% E-Mail: andreas.geiges@uni-stuttgart.de
% Date:  05/2015

%% WRAPPER:
% Evalulation of the expected conditional variance of the prediction data
% for prior data given the observation data

% INPUT:            NAME                                        DIMENSION
% ===================================================================================
% ctrl              Control structure containing various info   STUCTURE
%     .no_err_marg  flag if marginalizing over obs errro        1
%                   (see Leube et al. 2012, WRR)
%     .n_para       number of parallel computations for         1
%                   memory management
%     .sys.memory   memory of the system in in bytes            1
%     .warn_ESS     minimal required ESS                        1
%     .no_warning   do not display convergence warnings         1
%
% prior_data        data sample                                 DIM:N_MC
% obs_data          sample of observations                      DIM:N_MEAS
% obs_err_std       standart deviation of measurement error     DIM:1
%
% pred_data         prediction data                             DIM:N_MC
% pred_err_std      standart deviation of prediction error      DIM:1
%
% OUTPUT:           NAME                                        DIMENSION
% ===================================================================================
% E_cond_var        expected conditional prediction variance    1:N_MEAS
%                   given the observation data sample
% ESS               Effective sample size for each condition    1:N_MEAS
%                   sample (given data)

%% When using this code please cite the original paper:
% BIBTeX format:
%@article{Leube2012,
%	author = "P. Leube and A. Geiges and W. Nowak",
%	doi = "10.1029/2010WR010137",
%	journal = "Water Resources Research",
%	note = "{W02501}",
%	number = "2",
%	title = "{Bayesian assessment of the expected data impact on prediction confidence in optimal sampling design}",
%	volume = "48",
%	year = 2012
%}


%% INIT
[n_dim_data, n_mc  ] = size(prior_data);
[n_dim_obs  ,n_meas] = size(obs_data);

if ~isfield(ctrl,'warn_ESS')
    ctrl.warn_ESS = 50;
end

if n_dim_data ~= n_dim_obs
    error('Dimension of data disagree')
end

if ~strct_bool_check(ctrl,'n_para')
    ctrl.n_para = 1; % no paralell computing
end

% calculation in how many parts the calculation needs to be split in order
% to match the computer memory
[n_split, part_start, part_end, n_part] = get_n_splits(ctrl, n_mc,n_meas,ctrl.n_para);

if n_split > 1 && ctrl.n_para == 1
    disp(['splitting in ' num2str(n_split) ' pieces'])
end

%% Calculation loop
if n_split > 1 && ctrl.n_para == 1
    fprintf(['Eval split: 0/' num2str(n_split)])
end
for t = 1:n_split
    if n_split > 1 && ctrl.n_para == 1
        if num2str(n_split) > 9
            fprintf(['\b\b\b\b' num2str(t) '/' num2str(n_split)])
        else
            fprintf(['\b\b\b' num2str(t) '/' num2str(n_split)])
        end
    end
    obs_data_part = obs_data(:,part_start(t):part_end(t));
    
    % calculation of weighting matrix
    [weights, ESS, sumSqrWeights,ttime] = predia_weight_matrix(ctrl, prior_data,obs_data_part, obs_err_std,prior_weight);
    
    % calculation of weighted variance
    cond_var(part_start(t):part_end(t))      = weighted_cond_var(ctrl, weights,sumSqrWeights,pred_data,'weights');
    
end
if n_split > 1 && ctrl.n_para == 1
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
end

cond_var = cond_var(~isnan(cond_var));
E_cond_var = mean(cond_var);

if min(ESS) < ctrl.warn_ESS
    n_crit = sum(ESS < ctrl.warn_ESS);
    if ~strct_bool_check(ctrl,'no_warning')
        disp(['ESS lower than ' num2str(ctrl.warn_ESS) ' in ' num2str(n_crit/n_meas*100) '% of obs. realizations'])
    end
end
