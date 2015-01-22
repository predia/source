%% CONDITIONAL VARINANCE
%
% version 1 / Jan 15 / AGeiges WNowak
%
% calculates the conditional variance form the given PreDIAMatrix

function varout = weighted_cond_var(ctrl, weight_matrix,sumSqrWeights,targ_data,delete_name)

%% describption of phi_arguments
% phi_arguments{1} = data_meas;

%%
if exist('delete_name','var')
    assignin('caller', delete_name, []);
end

% nan_idx = find(isnan(sumSqrWeights));
% weight_matrix(nan_idx,:) = [];
% sumSqrWeights(nan_idx)   = [];
% sumSqrWeights(isnan(sumSqrWeights)) = [];

if isfield(ctrl,'transposed') && ctrl.transposed
    varout = (1./(1-sumSqrWeights)) .* ( targ_data.^2 * weight_matrix - ( targ_data* weight_matrix).^2) ;
else
    varout = repmat((1./(1-sumSqrWeights)),1,size(targ_data,1)) .* ( weight_matrix *  targ_data.^2' - (weight_matrix *  targ_data').^2) ;
end


% end