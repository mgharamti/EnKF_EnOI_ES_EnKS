function [Ur, RMSE, RMSF, RMSA, AESP, EnV1, EnV2, EnV3] = ES(tag, Ne)

if nargin < 1
    tag = 1;
end

rng( 'default' )

% Get the observations, free run, and initialize the model
[model, Ur, Up, RMSE, Force, Y, H, R, Cy, pt, sig, p] = init_system;

%Prior Statistics:
B  = cov( Up' ) ;
t = 1 ;
Xa = mvnrnd( mean( Up,2 ),B,Ne )' ; 
Xf = Xa ; 

X_tel = zeros( model.Nx*Cy,Ne );
S_tel = zeros( p*Cy,Ne );
R_tel = zeros( p*Cy,p*Cy );
D_tel = zeros( p*Cy,Ne );

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

while( t <= Cy )
    
    % Propagation
    for i= 1:model.rep
        St = Force( :,model.rep*(t-1) + i );
        for e = 1:Ne
            Xf( :,e ) = HeatModel1D( Xf( :,e ),model,St ) + model.err*randn( model.Nx,1 ) ;
        end
    end
    
    S= H * ( Xf - repmat( mean( Xf,2 ),1,Ne ) ); 
    D= Y( :,t )*ones( 1,Ne ) + ( ( H*sig )*ones( 1,Ne ) ).*randn( p,Ne ) - H*Xf;
    
    X_tel( (t-1)*model.Nx+1:t*model.Nx,: )= Xf; 
    S_tel( (t-1)*p+1:t*p,: )= S;
    R_tel( (t-1)*p+1:t*p,(t-1)*p+1:t*p )= R;
    D_tel( (t-1)*p+1:t*p,: )= D;
    
    t = t + 1 ;
end

% Update
disp( [ 'Update, Time = ' num2str( model.dt*(t-1)*model.rep ) ' sec' ] ) ;

C_tel= S_tel*S_tel' + (Ne-1)*R_tel;

[ U1,S1,V1 ]= svd(C_tel);
Fi = find( cumsum( diag(S1)/sum(S1(:)) ) > 0.99,1 );
Ci_tel = V1( :,1:Fi )*diag( 1./diag( S1(1:Fi,1:Fi) ) )*U1( :,1:Fi )';

X5= ( eye(Ne)+S_tel'*Ci_tel*D_tel ); %Evensen-Style
X= X_tel*X5;
    
    
% Calculations
for t= 1:Cy
    RMSF( t ) = 1/sqrt( model.Nx ) * norm( mean( X_tel( (t-1)*model.Nx+1:t*model.Nx,: ),2 ) - Ur( :,t*model.rep ) ) ; 
    RMSA( t ) = 1/sqrt( model.Nx ) * norm( mean( X( (t-1)*model.Nx+1:t*model.Nx,: ),2 ) - Ur( :,t*model.rep ) ) ; 
    AESP( t ) = sqrt( 1/model.Nx * sum( var( X_tel( (t-1)*model.Nx+1:t*model.Nx,: ),0,2 ) ) ) ; 
    EnV1( t,: ) = X_tel( (t-1)*model.Nx+vb1,: ) ;
    EnV2( t,: ) = X_tel( (t-1)*model.Nx+vb2,: ) ;
    EnV3( t,: ) = X_tel( (t-1)*model.Nx+vb3,: ) ;
end
    




% Plotting:
if tag

    figure( 'uni','pi','pos',[ 300,700,600,500 ] )

    h1= plot( model.dt*(1:20)*model.rep,RMSE,'k','LineWidth',2 ) ; hold on
    h2= plot( model.dt*(1:20)*model.rep,RMSF,'b','LineWidth',2 ) ; grid on
    h3= plot( model.dt*(1:20)*model.rep,RMSA,'r','LineWidth',2 ) ; xlabel( 'Time (sec)' )
    h4= plot( model.dt*(1:20)*model.rep,AESP,'g','LineWidth',2 ) ; ylabel( 'ES Statistics' )
    set( gca,'FontSize',16 ); Lg= legend( [ h1,h2,h3,h4 ],'RMS: Free-Run','RMS: Forecast','RMS: Analysis', ...
             'AES: Spread','Location','Best' ); set( Lg,'EdgeColor','w' )
    title( [ '$\widehat{RMSE}$ = ' num2str( mean( RMSE ) ) ', $\widehat{RMSA}$ = ' num2str( mean( RMSA ) ) ], ...
             'FontSize',20,'Interpreter','Latex' )


    figure( 'uni','pi','pos',[ 300,700,600,800 ] )

    subplot( 311 )
    hE= plot( model.dt*(1:20)*model.rep,EnV1,'--','Color',[ 0.5,0.5,0.5 ] ); hold on; grid on
    hM= plot( model.dt*(1:20)*model.rep,mean( EnV1,2 ),'-.r','LineWidth',2 );
    hR= plot( model.dt*(1:20)*model.rep,Ur( vb1,pt:pt:end ),'-b','LineWidth',2 );
    xlabel( 'Time (sec)' ); ylabel( [ 'Variable #' num2str(vb1) ] ); set( gca,'FontSize',16 )
    title( 'Ensemble Smoother: Evolution of Variables','FontSize',18 )
    legend( [ hR,hE(1),hM ],'Reference','ES Ensemble Members','ES Ensemble Mean','Location','Best' )

    subplot( 312 )
    hE= plot( model.dt*(1:20)*model.rep,EnV2,'--','Color',[ 0.5,0.5,0.5 ] ); hold on; grid on
    hM= plot( model.dt*(1:20)*model.rep,mean( EnV2,2 ),'-.r','LineWidth',2 );
    hR= plot( model.dt*(1:20)*model.rep,Ur( vb2,pt:pt:end ),'-b','LineWidth',2 );
    xlabel( 'Time (sec)' ); ylabel( [ 'Variable #' num2str(vb2) ] ); set( gca,'FontSize',16 )
    legend( [ hR,hE(1),hM ],'Reference','ES Ensemble Members','ES Ensemble Mean','Location','Best' )

    subplot( 313 )
    hE= plot( model.dt*(1:20)*model.rep,EnV3,'--','Color',[ 0.5,0.5,0.5 ] ); hold on; grid on
    hM= plot( model.dt*(1:20)*model.rep,mean( EnV3,2 ),'-.r','LineWidth',2 );
    hR= plot( model.dt*(1:20)*model.rep,Ur( vb3,pt:pt:end ),'-b','LineWidth',2 );
    xlabel( 'Time (sec)' ); ylabel( [ 'Variable #' num2str(vb3) ] ); set( gca,'FontSize',16 )
    legend( [ hR,hE(1),hM ],'Reference','ES Ensemble Members','ES Ensemble Mean','Location','Best' )
    
end