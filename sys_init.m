function [X, F] = sys_init(model)

    if strcmp(model.type, 'diffusion')

        X = 5 * ones(model.Nx, 1);        
        F = sin(model.rx); 
        s = 5 * randn(1, model.Cy);
        
        F = F * s;


    elseif strcmp(model.type, 'advection')

        X = sin(model.rx' * 5) + 3; 
        
        F = smooth( repmat((model.sc * sin(1:model.Nx)'), 1, model.Cy) .* randn(model.Nx, model.Cy) ); 
        
        F = reshape(F, model.Nx, model.Cy);
        
    else
        error('Model name is not available!')

    end

end