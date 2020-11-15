%
% SCRIPT_NAME: 
%   DA_EnKF_EnOI_ES_EnKS: [1] Ensemble Kalman Filter - [2] Ensemble Optimal 
%   Interpolation - [3] Ensemble Smoother - [4] Ensemble Kalman Smoother  
%
% Syntax:  
%   EnKF_EnOI_ES_EnKS.m 
%
% Inputs:
%          (1) eqn [1/2: diffusion/advection]
%          (2) tag [0/1: Don't/Do show statistics for each Scheme] 
%          (3) lag-0 for EnKS
%          (4) N: ensemble size
%
% Outputs: 
%          (1) var02: Root-Mean-Squared-Errors for Free Run
%          (2) var03: Root-Mean-Squared-Errors for Forecast
%          (3) var04: Root-Mean-Squared-Errors for Analysis
%          (4) var05: Average-Ensemble-Spread
%          (5) var06: Diagnostics for Var #1
%          (6) var07: Diagnostics for Var #2
%
% Example: 
%   An example is readily available 
%
% Author:       Mohamad El Gharamti
% Work address: National Center for Atmospheric Research
% Email:        gharamti@ucar.edu 
% Website:      http://dart.ucar.edu/
% Modified:     Nov. 14, 2020


Input0= 'Choose model form: [diffusion=1, advection=2]: ';
Input1= 'Display Output for each scheme? [yes=1/no=0]: ';
Input2= 'Length of the lagged smoothing window? [1/2/3/...]: ';
Input3= 'Ensemble size? [10/20/30/...]: ';

eqn= input(Input0); if isempty(eqn), eqn= 1;  end
tag= input(Input1); if isempty(tag), tag= 0;  end
lag= input(Input2); if isempty(lag), lag= 3;  end
N  = input(Input3); if isempty(N  ), N  = 100; end

switch eqn 
    case 1
        eqn = 'diffusion';
        
    case 2
        eqn = 'advection';
end          


fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          1. EnKF_DA !!!           \n')
fprintf(     '========= ************** ========= \n\n')

[~, RMSF_EnKF, RMSA_EnKF, AESP_EnKF, EnV1_EnKF, EnV2_EnKF] = EnKF(tag, eqn, N);


fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          2. EnOI_DA !!!           \n')
fprintf(     '========= ************** ========= \n\n')

[~, RMSF_EnOI, RMSA_EnOI, EnV1_EnOI, EnV2_EnOI] = EnOI(tag, eqn);


fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          3. ES_DA   !!!           \n')
fprintf(     '========= ************** ========= \n\n')

[RMSE, RMSF_ES, RMSA_ES, AESP_ES, EnV1_ES, EnV2_ES] = ES(tag, eqn, N);


fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          4. EnKS_DA !!!           \n')
fprintf(     '========= ************** ========= \n\n')

[~, RMSF_EnKS, RMSA_EnKS, AESP_EnKS, EnV1_EnKS, EnV2_EnKS] = EnKS(tag, eqn, lag, N);


fprintf( '\n\n========= ********************** ========= \n')
fprintf(     '          :OUTPUT VISUALIZATION:           \n')
fprintf(     '========= ********************** ========= \n\n')

%% 
figure('pos', [300, 700, 950, 540])

plot(RMSE,      '-','Color', [.7, .7, .7], 'LineWidth', 4) ; hold on
plot(RMSA_EnOI, '-xk', 'LineWidth', 2, 'MarkerSize', 8) ; hold on
plot(RMSA_ES  , '-xb', 'LineWidth', 2, 'MarkerSize', 8) ; grid on
plot(RMSA_EnKF, '-xr', 'LineWidth', 2, 'MarkerSize', 8) ;
plot(RMSA_EnKS, '-xg', 'LineWidth', 2, 'MarkerSize', 8) ;
set( gca,'FontSize',18 ); 

xlabel('DA Cycles', 'FontSize', 18)
ylabel('Posterior RMSE', 'FontSize', 18)

title( [ 'Averages: Free = '        sprintf( '%.3f', mean( RMSE ) ), ...
         ', EnOI = '      sprintf( '%.3f', mean( RMSA_EnOI ) ), ...
         ', ES = '   sprintf( '%.3f', mean( RMSA_ES ) ) ...
         ', EnKF = ' sprintf( '%.3f', mean( RMSA_EnKF ) ), ...
         ', EnKS = ' sprintf( '%.3f', mean( RMSA_EnKS ) ) ], ...
         'FontSize', 20)
     
lgd = legend('Free Run', 'EnOI', 'ES', 'EnKF','EnKS', 'Location', 'Best'); 

title(lgd, {eqn, [ 'Ensemble Size: ' num2str(N) ], ['lag: ' num2str(lag)]})