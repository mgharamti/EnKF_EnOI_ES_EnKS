function U = HeatModel1D(U, model, F)

    % A simple 1D (linear) heat diffusion model
    
    G = -2 * diag( ones( model.Nx,1 ),0 ) + diag( ones( model.Nx - 1,1 ),1 ) + ...
        diag( ones( model.Nx - 1,1 ) , -1 ) ; 
    
    G( 1,1 )               = -1; 
    G( model.Nx,model.Nx ) = -1;   
    
    % you can spy(G) to see how does it look like!
    
    M = model.alpha * model.dt/model.dx^2 * G + eye( model.Nx ); 
 
    U = max(0, M*U + model.dt*F);
    
end






