clear 
close all
clc


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
    
    plot(model.xi, U, '-b'); grid on
    set(gca, 'FontSize', 14, 'YLim', [1, 10], 'XLim', [model.xi(1), model.xi(model.Nx)])
    xlabel('Length (m)', 'FontSize', 16)
    ylabel('Temperature', 'FontSize', 16)
    title(['Time: ' num2str(model.dt*t) ' sec'], 'FontSize', 20)
    drawnow
end