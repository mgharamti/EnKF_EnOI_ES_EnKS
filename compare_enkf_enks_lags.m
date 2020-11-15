clear 
close all
clc

rng(1)

%% Test a large range of lags and see how the EnKS performs relative to the EnKF. 
tag      = 0; 
eqn      = 'advection';
ens_size = [20, 100];
LAG      = 1:10;
N        = length(ens_size);
L        = length(LAG);


fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          1. EnKF_DA !!!           \n')
fprintf(     '========= ************** ========= \n\n')

for e = 1:N
    [~, ~, RMSA_EnKF(:, e)]= EnKF(tag, eqn, ens_size(e)); %#ok
end


fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          2. EnKS_DA !!!           \n')
fprintf(     '========= ************** ========= \n\n')

cc = zeros(L, 1);
for e = 1:N
    for k = 1:L
        tic
        [~, ~, RMSA_EnKS(:, k, e)]= EnKS(tag, eqn, LAG(k), ens_size(e)); %#ok
        cc(k, e) = toc;
    end
end


fprintf( '\n\n========= ********************** ========= \n')
fprintf(     '          :OUTPUT VISUALIZATION:           \n')
fprintf(     '========= ********************** ========= \n\n')

enkf = mean(RMSA_EnKF, 1);
enks = squeeze(mean(RMSA_EnKS, 1));

%% 
figure( 'uni','pi','pos',[300, 700, 650, 500] )

plot(LAG, ones(k, 1) * mean(enkf(1)), '-k', 'LineWidth', 2); hold on
plot(LAG, enks(:, 1),'-b', 'LineWidth', 2); grid on

plot(LAG, ones(k, 1) * mean(enkf(2)), '--k', 'LineWidth', 2); hold on
plot(LAG, enks(:, 2),'--b', 'LineWidth', 2); grid on

set(gca, 'FontSize', 16, 'XLim', [LAG(1), LAG(L)])
xlabel('Smoother Lag', 'FontSize', 18)
ylabel('Posterior RMSE', 'FontSize', 18)

title('Performance Comparison for different Ensemble Size', 'FontSize', 20)

Lg = legend([ 'EnKF, N: ' num2str(ens_size(1)) ], ...
       [ 'EnKS, N: ' num2str(ens_size(1)) ], ...
       [ 'EnKF, N: ' num2str(ens_size(2)) ], ...
       [ 'EnKS, N: ' num2str(ens_size(2)) ], 'Location', 'NorthWest', 'FontSize', 16);
title(Lg, eqn)