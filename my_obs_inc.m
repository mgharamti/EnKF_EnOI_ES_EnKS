function obs_increments =  my_obs_inc(ensemble, prior_mean, prior_var, observation, obs_error_var)

% Compute the posterior variance
post_var = 1 / (1 / prior_var + 1 / obs_error_var);

% Compute posterior mean
post_mean = post_var * (prior_mean / prior_var + observation / obs_error_var);

% Shift the prior ensemble to have the posterior mean
updated_ensemble = ensemble - prior_mean + post_mean;

% Contract the ensemble to have the posterior_variance
var_ratio = post_var / prior_var;
updated_ensemble = sqrt(var_ratio) * (updated_ensemble - post_mean) + post_mean;

% Compute the increments
obs_increments = updated_ensemble - ensemble; 