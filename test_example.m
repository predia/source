%% TEST FILE FOR ILLUSTRATION

% This file show how the internal predia files are used in a test example.

clear all
n_mc = 30000;


%%% SIMPLE NON LINEAR MODEL %%%%%%%%%%%%
% with relevant parameters 1 and 7
% other parameter are only pseudo parameter
% parameter 10 is independent to the output


input = randn(1,n_mc);
%simluation different less dependent input parameter
for i= 2:9
    input(i,:) = input(i-1,:) + randn(1,n_mc).* 0.5;
end
input(10,:) =randn(1,n_mc);

output = (input(1,:)+1).^2 - (.5 .* input(7,:) + randn(1,n_mc).* .5).^3 + randn(1,n_mc).* 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    clf
    hold on
for i = 1:10

    scatter(input(11-i,:),output)
end

[X1,X2] = meshgrid(-5:.1:5,-5:.1:5)
%%
figure(4)
clf
surf(X1,X2,(X1+1).^2 - (.5 .* X2).^3)
hold on
scatter3(input(1,:),input(5,:),output,50,output,'filled')
%% single measurement
ctrl = [];
ctrl.err_marg = 0;
for i= 1:10
    [cond_var(i),ESS] = expect_cond_var_predia(ctrl, input(i,:),input(i,1:2000), 0.5,output,0.2);
end
prior_var = var(output)

figure(2)
clf
plot(prior_var- cond_var)
ylim([0 5])
legend('one obs.')

%% two measurements

for i= 1:10
    comb_input = [input(1,:) ; input(i,:)]; % fixed first observation
    [cond_var_comb2(i),ESS] = expect_cond_var_predia(ctrl, comb_input ,comb_input(:,1:2000), [0.5; 0.5] ,output,0.2);
end

%%
figure(2)
hold on

plot(1:10,prior_var- cond_var_comb2(1:10),'r')
xlabel('given input')
ylabel('variance reduction')

legend({'one obs.','two obs.'})