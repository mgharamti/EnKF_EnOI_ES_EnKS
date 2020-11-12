function [model, Ur, Up, RMSE, Force, Y, H, R, Cy, pt, sig, p, P] = init_system

%% Spacial Parameters:
model.Nx = 50 ; model.x0 = 0 ; model.xL = 10 ; model.alpha = 0.20 ;  
model.dx = ( model.xL - model.x0 ) / model.Nx ; 
model.xi = ( model.x0+model.dx/2 : model.dx : model.xL-model.dx/2 )' ;
disp( [ '1D linear Heat Equation with ' , num2str( model.Nx ) ' variables' ] )
%% Temporal Parameters:
model.dt = 0.04 ; model.T = 5000 ; 
disp( [ 'Modeling time is ' , num2str( model.dt*model.T ) ' sec with a total ' , ...
        num2str( model.T ) ' time steps' ] )

    
%% INITIALIZATION
U = 5*ones( model.Nx,1 ) ;
F = sin( model.xi ) ; pert = 5 .* randn( 1,model.T ) ;
Force = F * pert ;
Ur = zeros( model.Nx,model.T ) ;
for t = 1:model.T
    U = HeatModel1D( U,model,Force( :,t ) ) ;
    Ur( :,t ) = U ;
end


%% OBSERVATIONS
p = 3 ; pt = 50; 
P( 1:p ) = [ 15,25,35 ] ;
P = sort( P,'ascend' ) ; 
H = zeros( p,model.Nx ) ;

for i = 1:p
    H( i,P( i ) ) = 1 ;
end
sig = 0.25 * std( Ur,0,2 ) ;    %Observation perturbation == 25%
Y   = H * Ur( : , pt:pt:end ) ; 
disp( [ 'Consider ' , num2str( p ) ' observation points spatially every ' , ...
        num2str( pt ) ' time steps (i.e. ' , num2str( model.dt*pt ) ' sec)' ] )
    
Cy= size( Y,2 );
for i = 1:Cy
    Y( :,i ) = Y( :,i ) + ( H*sig ) .* randn( p,1 ) ;
end
R = diag( ( H*sig ).^2 ) ;


%% Simulate model error in the forecast model
xa = mean( Ur,2 ) ; U = xa ;
Up = zeros( model.Nx,model.T ) ;


model.alpha = model.alpha - 0.20*model.alpha .* randn ;     % Parameters perturbation == 20%
Force =  Force + 0.15*Force .* randn( model.Nx,model.T ) ;  % Forcing perturbation == 15%


%% Free-Run:    
for t = 1:model.T
    U = HeatModel1D( U,model,Force( :,t ) ) ;
    Up( :,t ) = U ; 
end

RMSE = zeros( Cy,1 ) ;
for t= 1:Cy
    RMSE(t)= 1/sqrt( model.Nx ) * norm( Ur( :,t*pt ) - Up( :,t*pt ) ) ;  
end
