%% CONDITIONAL ENTROPY
%
% version 1 / Jan 15 / AGeiges WNowak

function [h_zlx] = weighted_cond_ent(ctrl, weight_matrix,~,targ_data,delete_name)

%% describption of targ_data
% targ_data{1} = data;
% targ_data{2} = data_error_var;
% targ_data{3} = data_meas;

if exist('delete_name','var')
    assignin('caller', delete_name, []);
end
n_dim = size(targ_data,1);
for pred = 1:n_dim
%     if strcmp(func2str(ctrl.kde_estimator),'kde_gpuml')
%         h_zlx = zeros(size(weight_matrix,1),1);
%         for i = 1:size(weight_matrix,1)
%             [~, h_zlx(i)] = kde_gpuml(targ_data{1}(pred,1:ctrl.n_mc),targ_data{3}, ctrl.kde_eps,weight_matrix(i,:),(targ_data{2}).^(.5),[],[]);
%         end
%     else
        if n_dim == 1
            % delete matrix in this workspace 
            [~, h_zlx(:,pred)] = ctrl.kde_estimator(targ_data{1}(pred,1:ctrl.n_mc),targ_data{3}(pred,:), ctrl.kde_eps,weight_matrix,targ_data{2},[],[],'PreDIAMatrix');
        else
            [~, h_zlx(:,pred)] = ctrl.kde_estimator(targ_data{1}(pred,1:ctrl.n_mc),targ_data{3}(pred,:), ctrl.kde_eps,weight_matrix,targ_data{2},[],[]);
        end
    end
end


% [pdf, h_zlx(:,pred)] = kde_figtree(targ_data{1}(pred,1:ctrl.n_mc),targ_data{3}(pred,:), ctrl.kde_eps,PreDIAMatrix,(targ_data{2}).^(.5),[],[],'PreDIAMatrix');
% H_zlx = mean(H_zlx);
% hold on
%%
% scatter(targ_data{3}(pred,:),pdf(:,245))
