clear 
close all
clc

rng('default')

% Test the different linear models 

eqn = 'advection';

model  = sys_conf(eqn);
[X, F] = sys_init(model); 
    
XR = zeros(model.Nx, model.Cy);

for t = 1:model.Cy
    
    X = forecast_model(X, model, F(:, t)); 
    
    XR(:, t) = X;
    
    plot(model.rx, X, '-b'); grid on 
    set(gca, 'FontSize', 16, 'YLim', [0, 10], 'XLim', [model.rx(1), model.rx(model.Nx)])
    xlabel('Length (m)', 'FontSize', 18)
    ylabel('State', 'FontSize', 18)
    title([ 'Model: ', model.type, ', Time: ' sprintf('%10.3f', model.dt*t) ' sec'], 'FontSize', 20)
    
    drawnow
    
end 

B = cov(XR'); contour(B, 50); colorbar

set(gca, 'FontSize', 16)

title([ model.type ': Covariance model'], 'FontSize', 20)