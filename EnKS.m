function [Ur, RMSE, RMSF, RMSA, AESP, EnV1, EnV2, EnV3]= EnKS(tag, FL, Ne)

if nargin < 1
    tag = 1;
    FL = 5;
    Ne = 100;
end

rng( 'default' )

% Get the observations, free run, and initialize the model
[model, Ur, Up, RMSE, Force, Y, H, R, Cy, pt, sig, p] = init_system;

%Prior Statistics:
B  = cov( Up' ) ;
t = 1 ;
Xa = mvnrnd( mean( Up,2 ),B,Ne )' ; 
XA = zeros( model.Nx,Ne,Cy );

% Diagnostic Variables:
RMSF = zeros( Cy,1 ) ; 
RMSA = zeros( Cy,1 ) ; 
AESP = zeros( Cy,1 ) ;
EnV1 = zeros( Cy,Ne ) ;
EnV2 = zeros( Cy,Ne ) ;
EnV3 = zeros( Cy,Ne ) ;

vb1 = 10 ;
vb2 = 25 ;
vb3 = 40 ;

model.rep= pt;
model.err= 0.05;
model.lag= FL;

while( t <= Cy )
    
    % Propagation
    Xf= Xa; 
    for i= 1:model.rep
        St = Force( :,model.rep*(t-1) + i );
        for e = 1:Ne
            Xf( :,e ) = HeatModel1D( Xf( :,e ),model,St ) + model.err*randn( model.Nx,1 ) ;
        end
    end

    S= H * ( Xf - repmat( mean( Xf,2 ),1,Ne ) ); 
    D= Y( :,t )*ones( 1,Ne ) + ( ( H*sig )*ones( 1,Ne ) ).*randn( p,Ne ) - H*Xf;
    C= S*S' + (Ne-1)*R;

    [ U1,S1,V1 ]= svd(C);
    Fi = find( cumsum( diag(S1)/sum(S1(:)) ) > 0.99,1 );
    Ci = V1( :,1:Fi )*diag( 1./diag( S1(1:Fi,1:Fi) ) )*U1( :,1:Fi )';
    
    % Update
    disp( [ 'Update at observation step: ' num2str( t ) ' (Time = ' num2str( model.dt*t*model.rep ) ' sec)' ] ) ;
    X5= ( eye(Ne)+S'*Ci*D ); %Evensen-Style
    Xa= Xf*X5;
    
    XA( :,:,t )= Xa;

    if t > 1
        Sl= max( 1,t-model.lag+1 );
        Sr= t-1 ;
        disp( [ 'Smoothing the update between times: ' num2str( model.dt*Sl*model.rep ) ' and ' num2str( model.dt*Sr*model.rep ) ...
                ' sec using current observations' ] ) ;
        for tS = Sl:Sr
            XA( :,:,tS )= XA( :,:,tS )*X5;
        end
    end
    
    % Calculations
    RMSF( t ) = 1/sqrt( model.Nx ) * norm( mean( Xf,2 ) - Ur( :,t*model.rep ) ) ;  
    AESP( t ) = sqrt( 1/model.Nx * sum( var( Xf,0,2 ) ) ) ; 
    EnV1( t,: ) = Xf( vb1,: ) ;
    EnV2( t,: ) = Xf( vb2,: ) ;
    EnV3( t,: ) = Xf( vb3,: ) ;
    t = t + 1 ;
    
end

for t= 1:Cy
    RMSA( t ) = 1/sqrt( model.Nx ) * norm( mean( XA( :,:,t ),2 ) - Ur( :,t*model.rep ) ) ; 
end


% Plotting:
if tag

    figure( 'uni','pi','pos',[ 300,700,600,500 ] )

    h1= plot( model.dt*(1:Cy)*model.rep,RMSE,'k','LineWidth',2 ) ; hold on
    h2= plot( model.dt*(1:Cy)*model.rep,RMSF,'b','LineWidth',2 ) ; grid on
    h3= plot( model.dt*(1:Cy)*model.rep,RMSA,'r','LineWidth',2 ) ; xlabel( 'Time (sec)' )
    h4= plot( model.dt*(1:Cy)*model.rep,AESP,'g','LineWidth',2 ) ; ylabel( 'EnKS Statistics' )
    set( gca,'FontSize',16 ); Lg= legend( [ h1,h2,h3,h4 ],'RMS: Free-Run','RMS: Forecast','RMS: Analysis', ...
             'AES: Spread','Location','Best' ); set( Lg,'EdgeColor','w' )
    title( [ '$\widehat{RMSE}$ = ' num2str( mean( RMSE ) ) ', $\widehat{RMSA}$ = ' num2str( mean( RMSA ) ) ], ...
             'FontSize',20,'Interpreter','Latex' )


    figure( 'uni','pi','pos',[ 300,700,600,800 ] )

    subplot( 311 )
    hE= plot( model.dt*(1:Cy)*model.rep,EnV1,'--','Color',[ 0.5,0.5,0.5 ] ); hold on; grid on
    hM= plot( model.dt*(1:Cy)*model.rep,mean( EnV1,2 ),'-.r','LineWidth',2 );
    hR= plot( model.dt*(1:Cy)*model.rep,Ur( vb1,pt:pt:end ),'-b','LineWidth',2 );
    xlabel( 'Time (sec)' ); ylabel( [ 'Variable #' num2str(vb1) ] ); set( gca,'FontSize',16 )
    title( 'Ensemble Kalman Smoother: Evolution of Variables','FontSize',18 )
    legend( [ hR,hE(1),hM ],'Reference','EnKS Ensemble Members','EnKS Ensemble Mean','Location','Best' )

    subplot( 312 )
    hE= plot( model.dt*(1:Cy)*model.rep,EnV2,'--','Color',[ 0.5,0.5,0.5 ] ); hold on; grid on
    hM= plot( model.dt*(1:Cy)*model.rep,mean( EnV2,2 ),'-.r','LineWidth',2 );
    hR= plot( model.dt*(1:Cy)*model.rep,Ur( vb2,pt:pt:end ),'-b','LineWidth',2 );
    xlabel( 'Time (sec)' ); ylabel( [ 'Variable #' num2str(vb2) ] ); set( gca,'FontSize',16 )
    legend( [ hR,hE(1),hM ],'Reference','EnKS Ensemble Members','EnKS Ensemble Mean','Location','Best' )

    subplot( 313 )
    hE= plot( model.dt*(1:Cy)*model.rep,EnV3,'--','Color',[ 0.5,0.5,0.5 ] ); hold on; grid on
    hM= plot( model.dt*(1:Cy)*model.rep,mean( EnV3,2 ),'-.r','LineWidth',2 );
    hR= plot( model.dt*(1:Cy)*model.rep,Ur( vb3,pt:pt:end ),'-b','LineWidth',2 );
    xlabel( 'Time (sec)' ); ylabel( [ 'Variable #' num2str(vb3) ] ); set( gca,'FontSize',16 )
    legend( [ hR,hE(1),hM ],'Reference','EnKS Ensemble Members','EnKS Ensemble Mean','Location','Best' )
    
end