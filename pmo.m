function [model, XR, obs, F] = pmo(model)

    [X, F] = sys_init(model); 

    XR = zeros(model.Nx, model.Cy);
    for t = 1:model.Cy

        X = forecast_model(X, model, F(:, t)); 

        XR(:, t) = X;
    end

    if strcmp(model.type, 'diffusion')
        
        obs.num_obs_space = 3 ;                                         % no. of obs point in the domain
        obs.freq_obs_time = 50;                                         % frequency of obs (assimilate every # cycles)
        obs.obs_loc_space = [15, 25, 35];                               % obs locations
        obs.obs_err_sd    = 0.25 * std(XR, 0, 2);                       % Observation perturbation 
        obs.frwd_operator = zeros(obs.num_obs_space, model.Nx);         % obs matrix
        
        obs.vars_to_diag  = [25, 40];                                   % observed and un-observed variable

        disp( [ 'Consider ' , num2str(obs.num_obs_space) ' observation points spatially every ' , ...
                num2str(obs.freq_obs_time) ' time steps (i.e. ' , num2str(model.dt*obs.freq_obs_time) ' sec)' ] )


    elseif strcmp(model.type, 'advection')
        
        obs.num_obs_space = 9;                                     
        obs.freq_obs_time = 20;     
        obs.obs_loc_space = 10:10:90;              
        obs.obs_err_sd    = 0.40 * std(XR, 0, 2); 
        obs.frwd_operator = zeros(obs.num_obs_space, model.Nx );    
        
        obs.vars_to_diag  = [30, 73];   
        
        disp( [ 'Consider ' , num2str(obs.num_obs_space) ' observation points spatially every ' , ...
                num2str(obs.freq_obs_time) ' time steps (i.e. ' , num2str( model.dt*obs.freq_obs_time/( 24*3600 ) ) ' days)' ] )

    else
        error('Model name is not available!')
        
    end
    
    obs.da_steps = obs.freq_obs_time:obs.freq_obs_time:model.Cy;    % Assimilation time
    obs.da_Cy    = length(obs.da_steps);                            % DA cycles

    for i = 1:obs.num_obs_space
        obs.frwd_operator(i, obs.obs_loc_space(i)) = 1;
    end

    obs.data  = obs.frwd_operator * XR(:, obs.freq_obs_time:obs.freq_obs_time:model.Cy);   % Observations (in space and time)

    for i = 1:obs.da_Cy
        obs.data(:, i) = obs.data(:, i) + (obs.frwd_operator * obs.obs_err_sd) .* randn(obs.num_obs_space, 1);   % Perturbed observations
    end
    obs.obs_err_var = diag((obs.frwd_operator * obs.obs_err_sd).^2);


end