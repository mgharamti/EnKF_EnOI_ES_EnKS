function [model, obs, RMSE, F, Xa, XP] = da_init(model, XR, F, obs, Ne)

    % perturb the initial state: Mean of the truth
    X = mean(XR, 2);  

    if strcmp(model.type, 'diffusion')
        
        % perturb model parameter
        model.alpha = model.alpha - 0.20 * model.alpha .* randn;       
        
        % perturb forcing
        F = F + 0.15 * F .* randn(model.Nx, model.Cy);                
        

    elseif strcmp(model.type, 'advection')

        % perturbed transport parameters
        model.phi = model.phi + 0.20 * model.phi .* randn;     

        % perturbed source and velocity components
        model.sc = 5.00e-06;   
        model.uc = 0.43e-04;  
        
        F = smooth( repmat((model.sc * sin(1:model.Nx)'), 1, model.Cy) .* randn(model.Nx, model.Cy) ); 
        
        F = reshape(F, model.Nx, model.Cy);
        

    else
        error('Model name is not available!')
        
    end
    
    % perturbed model run (free run; using perturbed params, forcing and IC)
    XP = zeros(model.Nx, model.Cy) ;
    for t = 1:model.Cy

        X = forecast_model(X, model, F(:, t)); 

        XP(:, t) = X;
    end

    % Free-run RMSE
    RMSE = zeros(obs.da_Cy, 1);
    for t = 1:obs.da_Cy

        j = t * obs.freq_obs_time;

        RMSE(t) = 1/sqrt(model.Nx) * norm(XR(:, j) - XP(:, j));  
    end
    
    % Initial Ensemble
    if Ne > 1
        B  = cov(XP');
        Xa = mvnrnd(mean(XP, 2), B, Ne)';  
    else
        % EnOI case
        Xa = mean(XP, 2);
    end
    
    % additive model noise parameter (as a form of inflation)
    obs.mod_err = 0.05;

end