function [cond_var,ESS] = expect_cond_ent_predia(ctrl, prior_data,obs_data, obs_err_std,pred_data,pred_err_std)

% version 1 / Jan 15 / AGeiges WNowak

%% INIT
[n_dim_data, n_mc  ] = size(prior_data);
[n_dim_obs  ,n_meas] = size(obs_data);
n_eval_pts = 1000;
ctrl.n_mc = n_mc;

if n_dim_data ~= n_dim_obs
    error('Dimension of data disagree')
end

if ~strct_flag_check(ctrl,'n_para')
    ctrl.n_para = 1; % no paralell computing
end

% calculation in how many parts the calculation needs to be split in order
% to match the computer memory
[n_split, part_start, part_end, n_part] = get_n_splits(ctrl, n_mc,n_meas,ctrl.n_para);

% numerical calculation of prediction pde
pred_data_cell{1} = pred_data; % source pts
pred_data_cell{2} = pred_err_std; % kernel width
pred_data_cell{3} = linspace(min(pred_data-5.*pred_err_std), max(pred_data+ 5.*pred_err_std), n_eval_pts); % density evaluation pts

%% Calculation loop
for t = 1:n_split
    
    obs_data_part = obs_data(:,part_start(t):part_end(t));
    
    % calculation of weighting matrix
    [weights, AESS, sumSqrWeights,ttime,ESS] = predia_weight_matrix(ctrl, prior_data,obs_data_part, obs_err_std);
    
    % calculation of weighted variance
    cond_var(part_start(t):part_end(t))      = weighted_cond_ent(ctrl, weights,sumSqrWeights,pred_data_cell,'weights');
    
end

if min(ESS) < 100
    n_crit = sum(ESS < 100);
    warning(['Effective sample size is ' num2str(n_crit) 'X close to critical value for proper computation of a variance measure'])
end
