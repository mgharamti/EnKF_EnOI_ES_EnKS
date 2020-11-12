function [Ur, RMSE, RMSF, RMSA, EnV1, EnV2, EnV3] = EnOI(tag, Ne)

if nargin < 1
    tag = 1;
end

rng( 'default' )

% Get the observations, free run, and initialize the model
[model, Ur, Up, RMSE, Force, Y, H, R, Cy, pt] = init_system;

%Prior Statistics:
Xa = mean( Up,2 ) ;  

% Static "historical" ensmeble
t = 1 ;
Xo = Up( :,randperm( model.T,Ne ) ) ;
Ao = ( Xo - repmat( mean( Xo,2 ),1,Ne ) ) / sqrt( Ne-1 ) ;
Co = ( H*Ao )*( H*Ao )' + R ;

[ U1,S1,V1 ]= svd(Co);
Fi = find( cumsum( diag(S1)/sum(S1(:)) ) > 0.99,1 );
Ci = V1( :,1:Fi )*diag( 1./diag( S1(1:Fi,1:Fi) ) )*U1( :,1:Fi )';

Go = Ao * ( ( H*Ao )' * Ci ); 

% Diagnostic Variables:
RMSF = zeros( Cy,1 ) ; 
RMSA = zeros( Cy,1 ) ; 
EnV1 = zeros( Cy,1 ) ;
EnV2 = zeros( Cy,1 ) ;
EnV3 = zeros( Cy,1 ) ;

vb1 = 10 ;
vb2 = 25 ;
vb3 = 40 ;

model.rep= pt;
model.err= 0.0;

while( t <= Cy )
    
    % Propagation
    Xf= Xa; 
    for i= 1:model.rep
        St = Force( :,model.rep*(t-1) + i );
        Xf = HeatModel1D( Xf,model,St ) + model.err*randn( model.Nx,1 ) ;
    end
    D= Y( :,t ) - H*Xf;
    
    % Update
    disp( [ 'Update at observation step: ' num2str( t ) ' (Time = ' num2str( model.dt*t*model.rep ) ' sec)' ] ) ;
    Xa= Xf + Go*D ;
    
    
    % Calculations
    RMSF( t ) = 1/sqrt( model.Nx ) * norm( Xf - Ur( :,t*model.rep ) ) ; 
    RMSA( t ) = 1/sqrt( model.Nx ) * norm( Xa - Ur( :,t*model.rep ) ) ; 
    EnV1( t ) = Xf( vb1 ) ;
    EnV2( t ) = Xf( vb2 ) ;
    EnV3( t ) = Xf( vb3 ) ;
    t = t + 1 ;
    
end


% Plotting:
if tag

    figure( 'uni','pi','pos',[ 300,700,600,500 ] )

    h1= plot( model.dt*(1:20)*model.rep,RMSE,'k','LineWidth',2 ) ; hold on
    h2= plot( model.dt*(1:20)*model.rep,RMSF,'b','LineWidth',2 ) ; grid on
    h3= plot( model.dt*(1:20)*model.rep,RMSA,'r','LineWidth',2 ) ; 
    xlabel( 'Time (sec)' ); ylabel( 'EnOI Statistics' )
    set( gca,'FontSize',16 ); Lg= legend( [ h1,h2,h3 ],'RMS: Free-Run','RMS: Forecast','RMS: Analysis', ...
             'Location','Best' ); set( Lg,'EdgeColor','w' )
    title( [ '$\widehat{RMSE}$ = ' num2str( mean( RMSE ) ) ', $\widehat{RMSA}$ = ' num2str( mean( RMSA ) ) ], ...
             'FontSize',20,'Interpreter','Latex' )


    figure( 'uni','pi','pos',[ 300,700,600,800 ] )

    subplot( 311 )
    hE= plot( model.dt*(1:20)*model.rep,EnV1,'-.r','LineWidth',2 ); hold on; grid on
    hR= plot( model.dt*(1:20)*model.rep,Ur( vb1,pt:pt:end ),'-b','LineWidth',2 );
    xlabel( 'Time (sec)' ); ylabel( [ 'Variable #' num2str(vb1) ] ); set( gca,'FontSize',16 )
    title( 'Ensemble Optimal Interpolation: Evolution of Variables','FontSize',18 )
    legend( [ hR,hE ],'Reference','EnOI Analysis','Location','Best' )

    subplot( 312 )
    hE= plot( model.dt*(1:20)*model.rep,EnV2,'-.r','LineWidth',2 ); hold on; grid on
    hR= plot( model.dt*(1:20)*model.rep,Ur( vb2,pt:pt:end ),'-b','LineWidth',2 );
    xlabel( 'Time (sec)' ); ylabel( [ 'Variable #' num2str(vb2) ] ); set( gca,'FontSize',16 )
    legend( [ hR,hE ],'Reference','EnOI Analysis','Location','Best' )

    subplot( 313 )
    hE= plot( model.dt*(1:20)*model.rep,EnV3,'-.r','LineWidth',2 ); hold on; grid on
    hR= plot( model.dt*(1:20)*model.rep,Ur( vb3,pt:pt:end ),'-b','LineWidth',2 );
    xlabel( 'Time (sec)' ); ylabel( [ 'Variable #' num2str(vb3) ] ); set( gca,'FontSize',16 )
    legend( [ hR,hE ],'Reference','EnOI Analysis','Location','Best' )
    
end