%
% SCRIPT_NAME: 
%   DA_EnKF_EnOI_ES_EnKS: [1] Ensemble Kalman Filter - [2] Ensemble Optimal 
%   Interpolation - [3] Ensemble Smoother - [4] Ensemble Kalman Smoother  
%
% Syntax:  
%   EnKF_EnOI_ES_EnKS.m 
%
% Inputs:
%   model: 1D-heat model (simple)
%   input: (1) tag [0/1: Don't/Do show statistics for each Scheme] 
%          (2) lag-0 for EnKS
% Outputs: (1) var01: Reference Trajectory for state
%          (2) var02: Root-Mean-Squared-Errors for Free Run
%          (3) var03: Root-Mean-Squared-Errors for Forecast
%          (4) var04: Root-Mean-Squared-Errors for Analysis
%          (5) var05: Average-Ensemble-Spread
%          (6) var06: Diagnostics for Var #1
%          (7) var07: Diagnostics for Var #2
%          (8) var08: Diagnostics for Var #3
%
% Example: 
%   An example is readily available 
%
% Author:       Mohamad El Gharamti
% Work address: Nansen Environmental and Remote Sensing Center
% Email:        mohamad.gharamti@nersc.no 
% Website:      http://www.nersc.no/
% Modified:     20-Apr-2016



Input1= 'Display Output for each scheme? [yes=1/no=0]: ';
Input2= 'Length of the lagged smoothing window? [1/2/3/...]: ';
Input3= 'Ensemble size? [10/20/30/...]: ';

tag= input(Input1); if isempty(tag), tag= 0;  end
lag= input(Input2); if isempty(lag), lag= 5;  end
N  = input(Input3); if isempty(N  ), N  = 20; end



fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          1. EnKF_DA !!!           \n')
fprintf(     '========= ************** ========= \n\n')

[~, ~, RMSF_EnKF, RMSA_EnKF, AESP_EnKF, EnV1_EnKF, EnV2_EnKF, EnV3_EnKF] = EnKF(tag, N);


fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          2. EnOI_DA !!!           \n')
fprintf(     '========= ************** ========= \n\n')

[~, ~, RMSF_EnOI, RMSA_EnOI, EnV1_EnOI, EnV2_EnOI, EnV3_EnOI] = EnOI(tag, N);


fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          3. ES_DA   !!!           \n')
fprintf(     '========= ************** ========= \n\n')

[Ur, RMSE, RMSF_ES, RMSA_ES, AESP_ES, EnV1_ES, EnV2_ES, EnV3_ES] = ES(tag, N);


fprintf( '\n\n========= ************** ========= \n')
fprintf(     '          4. EnKS_DA !!!           \n')
fprintf(     '========= ************** ========= \n\n')

[~, ~, RMSF_EnKS, RMSA_EnKS, AESP_EnKS, EnV1_EnKS, EnV2_EnKS, EnV3_EnKS] = EnKS(tag,lag, N);


fprintf( '\n\n========= ********************** ========= \n')
fprintf(     '          :OUTPUT VISUALIZATION:           \n')
fprintf(     '========= ********************** ========= \n\n')

%% 
figure('pos', [300, 700, 950, 540])

plot(RMSE,      '-','Color', [.7, .7, .7], 'LineWidth', 2) ; hold on
plot(RMSA_EnOI, '-xk', 'LineWidth', 2, 'MarkerSize', 8) ; hold on
plot(RMSA_ES  , '-xb', 'LineWidth', 2, 'MarkerSize', 8) ; grid on
plot(RMSA_EnKF, '-xr', 'LineWidth', 2, 'MarkerSize', 8) ;
plot(RMSA_EnKS, '-xg', 'LineWidth', 2, 'MarkerSize', 8) ;
set( gca,'FontSize',18 ); 

xlabel('DA Cycles', 'FontSize', 18)
ylabel('Analysis RMSE', 'FontSize', 18)

title( [ 'Averages: Free = '        sprintf( '%.3f', mean( RMSE ) ), ...
         ', EnOI = '      sprintf( '%.3f', mean( RMSA_EnOI ) ), ...
         ', ES = '   sprintf( '%.3f', mean( RMSA_ES ) ) ...
         ', EnKF = ' sprintf( '%.3f', mean( RMSA_EnKF ) ), ...
         ', EnKS = ' sprintf( '%.3f', mean( RMSA_EnKS ) ) ], ...
         'FontSize', 20)
     
lgd = legend('Free Run', 'EnOI', 'ES', 'EnKF','EnKS', 'Location', 'Best'); 

title(lgd, [ 'Ensemble Size: ' num2str(N) ])