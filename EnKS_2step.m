function [RMSE, RMSF, RMSA, AESP, EnV1, EnV2]= EnKS_2step(tag,FL, Ne)

if nargin < 1
    tag = 1;
    FL  = 1;
    eqn = 'advection';
    Ne  = 100;
end

rng( 'default' )

% System Config
model = sys_conf(eqn);

% Truth & obs
[model, XR, Y, F] = pmo(model);

% DA configure 
[model, Y, RMSE, F, Xa] = da_init(model, XR, F, Y, Ne);

% Diagnostic Variables:
RMSF = zeros(Y.da_Cy, 1); 
RMSA = zeros(Y.da_Cy, 1); 
AESP = zeros(Y.da_Cy, 1); 
EnV1 = zeros(Y.da_Cy, Ne);
EnV2 = zeros(Y.da_Cy, Ne);

% DA loop
model.lag = FL;

t  = 1;
XA = zeros(model.Nx, Ne, Y.da_Cy);

while t <= Y.da_Cy
    
    % Propagation
    Xf = Xa; 
    for i = 1:Y.freq_obs_time
        f = F(:, Y.freq_obs_time * (t-1) + i);
        
        for e = 1:Ne
            Xf(:, e) = forecast_model(Xf(:, e), model, f) + Y.mod_err * randn(model.Nx, 1);
        end
    end
    Xfm = mean(Xf, 2);
    
%     % inflate
%     for e = 1:Ne
%         Xf(:, e) = 1.02*(Xf(:, e) - mean(Xf, 2)) + mean(Xf, 2);
%     end

    % Update
    disp( [ 'Update at DA step: ' num2str(t) ] );
    
    Xa   = Xf;
    data = Y.data(:, t);
    for o = 1:Y.num_obs_space
        
        % obs inc.
        y  = Xa(Y.obs_loc_space(o), :);
        ym = sum(y) / Ne;
        yc = y - ym;
        yv = sum(yc.^2) / (Ne - 1);
        
        dy = my_obs_inc(y, ym, yv, data(o), Y.obs_err_var(o,o));
        
        for k = 1:model.Nx
            
            % state inc.
            x   = Xa(k, :);
            xm  = sum(x) / Ne;
            xc  = x - xm;
            xyc = sum(xc .* yc) / (Ne - 1); 
            
            Xa(k, :) = x + xyc/yv * dy;
        end
        
        % Smoothing
        if t > 1 && model.lag > 1
            Sl= max( 1,t-model.lag+1 );
            Sr= t-1 ;
            
            disp( [ 'Smoothing the update between DA steps: #' num2str(Sl) ' and #' num2str(Sr) ' using current observations (i.e., DA step #' num2str(t) ')' ] );
            
            for tS = Sl:Sr
                Xs = XA(:, :, tS);

                for k = 1:model.Nx

                    % state inc.
                    x   = Xs(k, :);
                    xm  = sum(x) / Ne;
                    xc  = x - xm;
                    xyc = sum(xc .* yc) / (Ne - 1); 

                    Xs(k, :) = x + xyc/yv * dy;
                end
                XA(:, :, tS) = Xs;
            end
        end

    end
    
    XA( :,:,t )= Xa;
    
    % Calculations
    RMSF(t)    = 1 / sqrt(model.Nx) * norm(Xfm - XR(:, t * Y.freq_obs_time)); 
    AESP(t)    = sqrt(1 / model.Nx * sum(var(Xf, 0, 2))) ; 
    EnV1(t, :) = Xf(Y.vars_to_diag(1), :);
    EnV2(t, :) = Xf(Y.vars_to_diag(2), :);
    
    t = t + 1 ; 
end

for t = 1:Y.da_Cy
    RMSA(t) = 1 / sqrt(model.Nx) * norm(mean(XA(:, :, t), 2) - XR(:, t * Y.freq_obs_time)); 
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
    title(Lg, [ 'Ensemble size: ' num2str(Ne) ', lag: ' num2str(model.lag)])
    
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
