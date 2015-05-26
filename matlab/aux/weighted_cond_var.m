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


if strct_bool_check(ctrl,'no_sample_variance') 
    varout =  weight_matrix *  targ_data.^2' - (weight_matrix *  targ_data').^2;
else
    varout = repmat((1./(1-sumSqrWeights)),1,size(targ_data,1)) .* ( weight_matrix *  targ_data.^2' - (weight_matrix *  targ_data').^2) ;
end


% end