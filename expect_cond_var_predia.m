function [cond_var,ESS] = expect_cond_var_predia(ctrl, prior_data,obs_data, obs_err_std,pred_data,pred_err_std)

% version 1 / Jan 15 / AGeiges WNowak

%% INIT
[n_dim_data, n_mc  ] = size(prior_data);
[n_dim_obs  ,n_meas] = size(obs_data);

if n_dim_data ~= n_dim_obs
    error('Dimension of data disagree')
end

if ~strct_flag_check(ctrl,'n_para')
    ctrl.n_para = 1; % no paralell computing
end

% calculation in how many parts the calculation needs to be split in order
% to match the computer memory
n_split = get_n_split(ctrl, n_mc,n_meas,ctrl.n_para);
if n_split>1
    split = fix(n_meas/n_split);
    part_start(1) = 1;
    part_end(1)   = split;
    n_part(1)     = split;
    for  t = 2:n_split
        part_start(t) = part_end(t-1)+1;
        if t == n_split ,
            part_end(t)   = n_meas;
        else
            part_end(t) = part_end(t-1)+split;
        end;
        n_part(t) =  part_end(t) - part_start(t) +1;
    end
    
else
    part_start = 1;
    part_end   = n_meas;
end


for t = 1:n_split
    
    %         prior_data_part = prior_data(:,part_start(t):part_end(t));
    
    %     if ctrl.rand_meas
    %     end
    obs_data_part = obs_data(:,part_start(t):part_end(t));
    
    [weights, AESS, sumSqrWeights,ttime,ESS] = predia_weight_matrix(ctrl, prior_data,obs_data_part, obs_err_std);
    
    cond_var(part_start(t):part_end(t)) = weighted_cond_var(ctrl, weights,sumSqrWeights,pred_data,'weights');
    
end

