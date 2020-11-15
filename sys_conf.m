function model = sys_conf(eqn)

model.type = eqn;

if strcmp(model.type, 'diffusion') 
  
    model.Nx    = 50; 
    model.x0    = 0; 
    model.xL    = 10; 
    model.alpha = 0.20;  
    model.dx    = ( model.xL - model.x0 ) / model.Nx; 
    model.rx    = ( model.x0+model.dx/2 : model.dx : model.xL-model.dx/2 )';
    
    disp( [ '1D linear diffusion with ' , num2str(model.Nx) ' variables' ] )

    model.dt = 0.04; 
    model.T  = 200;
    model.Cy = ceil(model.T / model.dt); 
    
    disp( [ 'Modeling time is ', num2str(model.T) ' sec with a total of ', num2str(model.Cy) ' time steps' ] )
    
        
elseif strcmp(model.type, 'advection')
    
    model.Nx = 100 ;                %no. of cells 
    model.Lx = 1000 ;               %1 km
    model.dx = model.Lx/model.Nx ; 
    model.rx = ( model.dx:model.dx:model.Lx );
    
    disp( [ '1D linear advection with ' , num2str(model.Nx) ' variables' ] )

    model.dt = 10*3600 ;                %time step (10 days)     
    model.T  = 2.3*365.25*24*3600 ;     %total simulation time: # year
    model.Cy = ceil(model.T / model.dt); 
    
    disp( [ 'Modeling time is ', num2str(model.T/(24 * 3600)) ' days with a total of ', num2str(model.Cy) ' time steps' ] )

    model.cL = 5 ;                  %left boundary concentration
    model.cR = 0 ;                  %right boundary conc.            
    
    Kd = 0.792;                     %distribution coeff.
    Ph = 0.334;                     %porosity
    Rh = 2.650;                     %retardation factor
    
    model.phi = Ph + (1 - Ph) * Rh * Kd;

    model.sc = 3.00e-06 ;           %Source constant
    model.uc = 1.18e-04 ;           %Velocity constant
 
    
else
    error('Model name is not available!')
    
end