function [RMSE, RMSF, RMSA, AESP, EnV1, EnV2] = ES(tag, eqn, Ne)

if nargin < 1
    tag = 1;
    eqn = 'diffusion';
    Ne  = 80;
end

rng( 'default' )

% System Config
model = sys_conf(eqn);

% Truth & obs
[model, XR, Y, F] = pmo(model);

% DA configure 
[model, Y, RMSE, F, Xf] = da_init(model, XR, F, Y, Ne);

X_tel = zeros(model.Nx * Y.da_Cy, Ne);
S_tel = zeros(Y.num_obs_space * Y.da_Cy, Ne);
R_tel = zeros(Y.num_obs_space * Y.da_Cy, Y.num_obs_space * Y.da_Cy);
D_tel = zeros(Y.num_obs_space * Y.da_Cy, Ne);

% Diagnostic Variables:
RMSF = zeros(Y.da_Cy, 1); 
RMSA = zeros(Y.da_Cy, 1); 
AESP = zeros(Y.da_Cy, 1);
EnV1 = zeros(Y.da_Cy, Ne);
EnV2 = zeros(Y.da_Cy, Ne);

% DA loop
t  = 1;

while t <= Y.da_Cy
    
    % Propagation
    for i = 1:Y.freq_obs_time
        f = F(:, Y.freq_obs_time * (t-1) + i);
        
        for e = 1:Ne
            Xf(:, e) = forecast_model(Xf(:, e), model, f) + Y.mod_err .* randn(model.Nx, 1);
        end
    end
    
    S = Y.frwd_operator * (Xf - repmat(mean(Xf, 2), 1, Ne)); 
    D = Y.data(:, t) * ones(1, Ne) + ( (Y.frwd_operator * Y.obs_err_sd) * ones(1, Ne)) .* randn(Y.num_obs_space, Ne) - Y.frwd_operator * Xf;
    
    X_tel((t-1) * model.Nx + 1 : t * model.Nx, :) = Xf; 
    S_tel((t-1) * Y.num_obs_space + 1 : t * Y.num_obs_space, :) = S;
    R_tel((t-1) * Y.num_obs_space + 1 : t * Y.num_obs_space, (t - 1) * Y.num_obs_space + 1 : t * Y.num_obs_space) = Y.obs_err_var;
    D_tel((t-1) * Y.num_obs_space + 1 : t * Y.num_obs_space, :) = D;
    
    t = t + 1 ;
end

% Update
disp( [ 'Update at DA step: ' num2str(t) ] );

C_tel = S_tel * S_tel' + (Ne - 1) * R_tel;

[U1, S1, V1] = svd(C_tel);

Fi     = find( cumsum(diag(S1) / sum(S1(:))) > 0.99, 1);
Ci_tel = V1(:, 1:Fi) * diag(1 ./ diag(S1(1:Fi, 1:Fi))) * U1(:, 1:Fi)';

X5 = eye(Ne) + S_tel' * Ci_tel * D_tel; %Evensen-Style
X  = X_tel * X5;
    
    
% Calculations
for t = 1:Y.da_Cy
    RMSF(t)    = 1/sqrt(model.Nx) * norm(mean(X_tel((t - 1) * model.Nx + 1 : t * model.Nx,: ), 2) - XR(:, t * Y.freq_obs_time)); 
    RMSA(t)    = 1/sqrt(model.Nx) * norm(mean(X((t - 1) * model.Nx + 1 : t * model.Nx,: ), 2) - XR(:, t * Y.freq_obs_time)); 
    AESP(t)    = sqrt(1 / model.Nx * sum(var(X_tel((t - 1) * model.Nx + 1 : t * model.Nx,: ), 0, 2))); 
    EnV1(t, :) = X_tel((t - 1) * model.Nx + Y.vars_to_diag(1), :);
    EnV2(t, :) = X_tel((t - 1) * model.Nx + Y.vars_to_diag(2), :);
end


%% Plotting:
if tag

    figure( 'uni','pi','pos',[300, 700, 1000, 400] )

    h1 = plot(1:Y.da_Cy, RMSE, 'k', 'LineWidth', 2) ; hold on
    h2 = plot(1:Y.da_Cy, RMSF, 'b', 'LineWidth', 2) ; grid on
    h3 = plot(1:Y.da_Cy, RMSA, 'r', 'LineWidth', 2) ; 
    h4 = plot(1:Y.da_Cy, AESP, 'g', 'LineWidth', 2) ; 
    
    xlabel('DA steps', 'FontSize', 18)
    ylabel('Statistics', 'FontSize', 18)
    
    set(gca, 'FontSize', 16); 
    
    Lg = legend( [h1, h2, h3, h4]   , ...
            'RMS: Free-Run'         , ...
            'Prior RMSE'            , ...
            'Posterior RMSE'        , ...
            'Prior Spread'          , 'Location', 'Best'); 
    title(Lg, [ 'Ensemble size: ' num2str(Ne) ])
    
    title( [ 'Averages: RMSE (free-run): ' sprintf('%.2f', mean(RMSE)) ...
             ', Prior RMSE: ' sprintf('%.2f', mean(RMSF)), ...
             ', Post RMSE: '  sprintf('%.2f', mean(RMSA)) ], 'FontSize', 20)


    figure('uni','pi','pos',[300, 700, 1300, 450])

    subplot(121)
    hE = plot(1:Y.da_Cy, EnV1                             , '--' , 'Color', [0.75, 0.75, 0.75]); hold on 
    hM = plot(1:Y.da_Cy, mean(EnV1, 2)                    , '-.r', 'LineWidth', 2); grid on
    hR = plot(1:Y.da_Cy, XR(Y.vars_to_diag(1), Y.da_steps), '-b' , 'LineWidth', 2);
    
    xlabel('DA steps', 'FontSize', 18)
    title(['Observed Variable: #' num2str(Y.vars_to_diag(1)) ...
        ', Average Abs. Bias: ' sprintf('%.3f', mean(abs(mean(EnV1, 2)'  - XR(Y.vars_to_diag(1), Y.da_steps)))) ], 'FontSize', 20)
    
    
    set(gca, 'FontSize', 16)
    
    legend([hR, hE(1), hM], 'Reference (truth)', 'Ensemble Members', 'Ensemble Mean', 'Location', 'Best')

    subplot(122)
    hE = plot(1:Y.da_Cy, EnV2                             , '--' , 'Color', [0.75, 0.75, 0.75]); hold on 
    hM = plot(1:Y.da_Cy, mean(EnV2, 2)                    , '-.r', 'LineWidth', 2); grid on
    hR = plot(1:Y.da_Cy, XR(Y.vars_to_diag(2), Y.da_steps), '-b' , 'LineWidth', 2);
    
    xlabel('DA steps', 'FontSize', 18)
    title(['Unobserved Variable: #' num2str(Y.vars_to_diag(2)) ...
        ', Average Abs. Bias: ' sprintf('%.3f', mean(abs(mean(EnV2, 2)'  - XR(Y.vars_to_diag(2), Y.da_steps)))) ], 'FontSize', 20)
    
    set(gca, 'FontSize', 16)
    
    legend([hR, hE(1), hM], 'Reference (truth)', 'Ensemble Members', 'Ensemble Mean', 'Location', 'Best')

    
end
