function X = forecast_model(X, model, Forcing)

    num_vars  = model.Nx; 
    spacing   = model.dx;
    time_step = model.dt;

    switch model.type
        
        case 'diffusion'
            velocity = model.alpha;
            
            % A simple 1D (linear) heat diffusion model
    
            G = -2 * diag(ones(num_vars    , 1),  0) + ...
                     diag(ones(num_vars - 1, 1),  1) + ...
                     diag(ones(num_vars - 1, 1), -1) ; 

            G(1, 1)               = -1; 
            G(num_vars, num_vars) = -1;   

            % you can spy(G) to see how does it look like!

            M = velocity * time_step / spacing^2 * G + eye(num_vars); 

            X = M * X + time_step * Forcing;
            
            
        case 'advection'
            velocity = smooth(smooth(abs(model.uc .* randn(model.Nx+1, 1))));
            permporo = model.phi;
            bound_L  = model.cL; 
            bound_R  = model.cR; 
            
            for i = 1:num_vars

                if i == 1 
                    X(i) = X(i) + time_step / permporo * ( ...
                            ( ( max(velocity(i)  , 0) * bound_L   ) + ( min(velocity(i)  , 0) * X(i) ) )   / spacing - ...
                            ( ( max(velocity(i+1), 0) * X(i)      ) + ( min(velocity(i+1), 0) * X(i+1) ) ) / spacing + Forcing(i) ) ;

                elseif i == num_vars
                    X(i) = X(i) + time_step / permporo * ( ... 
                            ( ( max(velocity(i)  , 0) * X(i-1) ) + ( min(velocity(i)  , 0) * X(i)      ) ) / spacing - ...
                            ( ( max(velocity(i+1), 0) * X(i)   ) + ( min(velocity(i+1), 0) * bound_R   ) ) / spacing + Forcing(i) ) ;

                else
                    X(i) = X(i) + time_step / permporo * ( ...
                            ( ( max(velocity(i)  , 0) * X(i-1) ) + ( min(velocity(i)  , 0) * X(i) )   ) / spacing - ...
                            ( ( max(velocity(i+1), 0) * X(i)   ) + ( min(velocity(i+1), 0) * X(i+1) ) ) / spacing + Forcing(i) ) ;
                end

            end
             
    end
    
    X = max(0.0, X);

end