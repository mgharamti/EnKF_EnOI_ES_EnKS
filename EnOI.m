function [RMSE, RMSF, RMSA, EnV1, EnV2] = EnOI(tag, eqn)

if nargin < 1
    tag = 1;
    eqn = 'diffusion';
end

rng( 'default' )

% System Config
model = sys_conf(eqn);

% Truth & obs
[model, XR, Y, F] = pmo(model);

% DA configure 
[model, Y, RMSE, F, Xa, XP] = da_init(model, XR, F, Y, 1);
 

% Static "historical" ensemble
Ns = 1000;  % no. of static states

Xo = XP(:, randperm(model.Cy, Ns));
Ao = (Xo - repmat(mean(Xo, 2), 1, Ns)) / sqrt(Ns - 1);
Co = (Y.frwd_operator * Ao) * (Y.frwd_operator * Ao)' + Y.obs_err_var;

[U1, S1, V1] = svd(Co);

Fi = find(cumsum(diag(S1) / sum(S1(:))) > 0.99, 1);
Ci = V1(:, 1:Fi) * diag(1 ./ diag(S1(1:Fi, 1:Fi))) * U1(:, 1:Fi)';
 
Go = Ao * ((Y.frwd_operator * Ao)' * Ci); 

% Diagnostic Variables:
RMSF = zeros(Y.da_Cy, 1); 
RMSA = zeros(Y.da_Cy, 1); 
EnV1 = zeros(Y.da_Cy, 1);
EnV2 = zeros(Y.da_Cy, 1);


% DA loop
t = 1;

while t <= Y.da_Cy
    
    % Propagation
    Xf = Xa; 
    for i = 1:Y.freq_obs_time
        f  = F(:, Y.freq_obs_time * (t-1) + i);
        Xf = forecast_model(Xf, model, f) + Y.mod_err * randn(model.Nx, 1);
    end
    D = Y.data(:, t) - Y.frwd_operator * Xf;
    
    % Update
    disp( [ 'Update at DA step: ' num2str(t) ] );
    Xa = Xf + Go * D; 
    
    % Calculations
    RMSF(t) = 1 / sqrt(model.Nx) * norm(Xf - XR(:, t * Y.freq_obs_time)); 
    RMSA(t) = 1 / sqrt(model.Nx) * norm(Xa - XR(:, t * Y.freq_obs_time)); 
    EnV1(t) = Xf(Y.vars_to_diag(1));
    EnV2(t) = Xf(Y.vars_to_diag(2));
    
    t = t + 1 ;
end


%% Plotting:
if tag

    figure( 'uni','pi','pos',[300, 700, 1000, 400] )

    h1 = plot(1:Y.da_Cy, RMSE, 'k', 'LineWidth', 2) ; hold on
    h2 = plot(1:Y.da_Cy, RMSF, 'b', 'LineWidth', 2) ; grid on
    h3 = plot(1:Y.da_Cy, RMSA, 'r', 'LineWidth', 2) ; 
    
    xlabel('DA steps', 'FontSize', 18)
    ylabel('Statistics', 'FontSize', 18)
    
    set(gca, 'FontSize', 16); 
    
    Lg = legend( [h1, h2, h3]       , ...
            'RMS: Free-Run'         , ...
            'Prior RMSE'            , ...
            'Posterior RMSE'        , 'Location', 'Best'); 
    
    title( [ 'Averages: RMSE (free-run): ' sprintf('%.2f', mean(RMSE)) ...
             ', Prior RMSE: ' sprintf('%.2f', mean(RMSF)), ...
             ', Post RMSE: '  sprintf('%.2f', mean(RMSA)) ], 'FontSize', 20)


    figure('uni','pi','pos',[300, 700, 1300, 450])

    subplot(121)
    hM = plot(1:Y.da_Cy, EnV1                             , '-.r', 'LineWidth', 2); hold on
    hR = plot(1:Y.da_Cy, XR(Y.vars_to_diag(1), Y.da_steps), '-b' , 'LineWidth', 2); grid on
    
    xlabel('DA steps', 'FontSize', 18)
    title(['Observed Variable: #' num2str(Y.vars_to_diag(1)) ...
        ', Average Abs. Bias: ' sprintf('%.3f', mean(abs(EnV1'  - XR(Y.vars_to_diag(1), Y.da_steps)))) ], 'FontSize', 20)
    
    
    set(gca, 'FontSize', 16)
    
    legend([hR, hM], 'Reference (truth)', 'Ensemble Mean', 'Location', 'Best')

    subplot(122)
    hM = plot(1:Y.da_Cy, EnV2                             , '-.r', 'LineWidth', 2); hold on
    hR = plot(1:Y.da_Cy, XR(Y.vars_to_diag(2), Y.da_steps), '-b' , 'LineWidth', 2); grid on
    
    xlabel('DA steps', 'FontSize', 18)
    title(['Unobserved Variable: #' num2str(Y.vars_to_diag(2)) ...
        ', Average Abs. Bias: ' sprintf('%.3f', mean(abs(EnV2'  - XR(Y.vars_to_diag(2), Y.da_steps)))) ], 'FontSize', 20)
    
    set(gca, 'FontSize', 16)
    
    legend([hR, hM], 'Reference (truth)', 'Ensemble Mean', 'Location', 'Best')

    
end
