%% TEST FILE FOR ILLUSTRATION

% This file show how the internal predia files are used in a test example.

clear all
n_mc = 20000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMPLE NON LINEAR MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% with relevant parameters 1 and 7
% other parameter are only pseudo parameter
% parameter 10 is independent to the output
input = randn(1,n_mc);
%simluation different less dependent input parameter
for i= 2:9
    input(i,:) = input(i-1,:) + randn(1,n_mc).*0.5;
end
input(10,:) =randn(1,n_mc);

% inputs are observable quanties, for which an measurement error needs to
% be defined for the latter data worth analysis
meas_err_std(1:4)  = 0.8;
meas_err_std(5:10) = 0.3;

% model prediction
output = (input(1,:)+1).^2 - (.5 .* input(7,:) + randn(1,n_mc).* .25).^3 + randn(1,n_mc).* 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    clf
    hold on
    n_plot = 1000;
for i = 1:10

    scatter(input(11-i,1:n_plot),output(1:n_plot))
end

[X1,X2] = meshgrid(-5:.1:5,-5:.1:5);

fid_out = fopen('input.data','w');
fwrite(fid_out,input(:),'double');
fclose(fid_out);

fid_out = fopen('output.data','w');
fwrite(fid_out,output,'double');
fclose(fid_out);

%%
figure(4)
clf
surf(X1,X2,(X1+1).^2 - (.5 .* X2).^3)
[X,Y,Z] = sphere(10);
X = X./10;
Y = Y./10;

hold on
for i = 1:1000
    surf(input(1,i)+X,input(5,i)+Y,output(i)+Z)
end
shading flat
light
view(-190,37)
% scatter3(input(1,:),input(5,:),output,50,output,'filled')
%% single measurement
ctrl = [];
% ctrl.no_err_marg = 0;


% Evaluation of the prior variance of the model output
prior_var = var(output);

% Evaluation of the posterior variance of the model output conditional on the input data 1 to 10
for i= 1:10
    tic
    [cond_var(i),ESS] = expect_cond_var_predia(ctrl, input(i,:),input(i,1:2000), meas_err_std(i),output);
    toc
end

% Plotting
figure(2)
hold on
plot(prior_var- cond_var)
ylim([0 5])
xlabel('given input')
ylabel('variance reduction')
legend('one obs.')
grid on



%% two measurements


for i= 1:10
    tic
    comb_input = input([1 i],:);  % fixed first observation in combination of any other observations
    [cond_var_comb2(i),ESS] = expect_cond_var_predia(ctrl, comb_input ,comb_input(:,1:2000), meas_err_std([1,i]) ,output);
    toc
end

figure(2)
hold on

plot(1:10,prior_var- cond_var_comb2(1:10),'r')
xlabel('given input')
ylabel('variance reduction')

legend({'one obs.','two obs.'})

return
%%

[weights, ESS, sumSqrWeights,ttime] = predia_weight_matrix(ctrl,comb_input(:,1:10),comb_input(:,1:2), meas_err_std([1,i]))
weights
tmp_cond_var = weighted_cond_var(ctrl, weights,sumSqrWeights,output(1:10))


[weights, ESS, sumSqrWeights,ttime] = predia_weight_matrix(ctrl,comb_input(:,1:20000),comb_input(:,1:2000), meas_err_std([1,i]));
mean(weighted_cond_var(ctrl, weights,sumSqrWeights,output(1:20000)))