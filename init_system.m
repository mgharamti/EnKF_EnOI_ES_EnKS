function [model, Ur, Up, RMSE, Force, Y, H, R, Cy, pt, sig, p, P] = init_system

%% Spatial Parameters:
model.Nx    = 50;       % total number of variables 
model.x0    = 0;        % length at first point
model.xL    = 10;       % full length 
model.alpha = 0.20;     % diffusion parameter

model.dx = ( model.xL - model.x0 ) / model.Nx;                          % spatial discretization
model.xi = ( model.x0+model.dx/2 : model.dx : model.xL-model.dx/2 )';   % length of intervals

disp( [ '1D linear Heat Equation with ' , num2str(model.Nx) ' variables' ] )

%% Temporal Parameters:
model.dt = 0.04;        % Time step
model.T = 5000;         % Total modeling steps

disp( [ 'Modeling time is ' , num2str( model.dt*model.T ) ' sec with a total ' , ...
        num2str( model.T ) ' time steps' ] )

    
%% INITIALIZATION
U = 5 * ones(model.Nx, 1);          % Initial state 

F = sin(model.xi);                  % determinstic forcing 

pert  = 5 .* randn(1, model.T);     % pertubation to be added to F
Force = F * pert ;                  % different forcing every time step

Ur = zeros( model.Nx,model.T ) ;    % Reference States (Truth)
for t = 1:model.T
    U = HeatModel1D(U, model, Force(:, t)) ;
    
    Ur(:, t) = U ;
end


%% OBSERVATIONS
p  = 3 ;                    % no. of obs point in the domain
pt = 50;                    % frequency of obs (assimilate every pt cycles)
P  = [15, 25, 35];          % obs locations
H  = zeros( p,model.Nx );   % obs matrix

for i = 1:p
    H(i, P(i) ) = 1;
end

sig = 0.25 * std(Ur, 0, 2);         % Observation perturbation == 25%
Y   = H * Ur( :, pt : pt : end);    % Observations (in space and time)

disp( [ 'Consider ' , num2str(p) ' observation points spatially every ' , ...
        num2str(pt) ' time steps (i.e. ' , num2str(model.dt*pt) ' sec)' ] )
    
Cy = size(Y, 2);    % DA cycles
for i = 1:Cy
    Y(:, i) = Y(:, i) + (H * sig) .* randn(p, 1);   % Perturbed observations
end
R = diag( ( H*sig ).^2 ) ;


%% Simulate model error in the forecast model
% perturb the initial state
xa = mean(Ur, 2); 
U  = xa ;
Up = zeros(model.Nx, model.T) ;

% perturb model parameters
model.alpha = model.alpha - 0.20*model.alpha .* randn ;     % Parameters perturbation == 20%
Force =  Force + 0.15*Force .* randn( model.Nx,model.T ) ;  % Forcing perturbation == 15%


%% Free-Run:    
for t = 1:model.T
    U = HeatModel1D(U, model, Force( :,t )) ;
    
    Up(:, t) = U ; 
end

RMSE = zeros( Cy,1 ) ;
for t = 1:Cy
    RMSE(t) = 1/sqrt(model.Nx) * norm( Ur(:, t*pt) - Up(:, t*pt) );  
end
