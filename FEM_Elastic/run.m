close all;
clear all;

NX = 100;
NY = 100;
DH = 6;
Nperiod = 2;
DefCoef = 2;
[X, Y] = GenerateGrid( NX, NY, DH, Nperiod, DefCoef );
[P, T] = ConstructTriangle( X, Y, NX, NY );
%plotTriangle( P, T );
%网格生成
% load homo_allpml.mat
% P = p;
% T = t;

NT = 500;





src = zeros( NT, 1 );


Nodes = size( P, 1 );

%{
Unew  = zeros( Nodes, 1 );
Uold  = zeros( Nodes, 1 );
U     = zeros( Nodes, 1 );
%}

Unewx  = zeros( Nodes, 1 );
Uoldx  = zeros( Nodes, 1 );
Ux     = zeros( Nodes, 1 );

Unewy  = zeros( Nodes, 1 );
Uoldy  = zeros( Nodes, 1 );
Uy     = zeros( Nodes, 1 );

display( "=======产生高斯点=======" );
n = 8;
[r, s, w] = genGaussPoint( n );
Enodes = size( T, 2 );


Me = zeros( Enodes, Enodes );
Me = CalMe( r, s, w, T );

display( "=======雅可比系数求解=======" );
JField = CalJacobian( P, T );
display( "=======质量矩阵组装=======" );
M = AssembleMassMatrix ( Me, JField, P, T );
display( "=======刚度矩阵组装=======" );
[Sxx, Sxy, Syx, Syy] = AssembleStiffMatrix( JField, P, T, r, s, w );


CFL = 0.3;

DT = 0.001;

PPW = 12;

Vp = CFL*DH/DT;
Vs = 0.5774*Vp;
rho= 1500.0;

f0 = Vs/(PPW*DH);
mu = rho * ( Vs * Vs );
lam = rho * ( Vp * Vp - 2.0 * Vs * Vs );
lam2mu = lam + 2 * mu;

display( "=======矩阵求逆=======" );

Minv = DT ^ 2 * inv( M ) / rho;
srcX = NX/2;
srcY = NY/2;


tic

F1 = zeros( Nodes, 1 );
F2 = zeros( Nodes, 1 );
for it = 1 : NT
    display( it );
    src( it ) = ( 1 - 2 * ( pi * f0 * ( it * DT - 1 / f0) ) ^ 2 ) * exp ( - ( pi * f0 * ( it * DT - 1 / f0 ) ) ^ 2 );
    
    F2( srcX + srcY * NX, 1 ) = src( it );
    Unewx = Minv * ( F1 - ( lam2mu * Sxx * Ux + lam * Syx * Uy + mu * Sxy * Uy + mu * Syy * Ux ) ) +  2 * Ux - Uoldx; 
    Unewy = Minv * ( F2 - ( mu * Sxx * Uy + mu * Syx * Ux + lam * Sxy * Ux + lam2mu * Syy * Uy ) ) +  2 * Uy - Uoldy; 
    

    %Unew = Minv / rho * ( F - ( lam * ( Sxx + Syy ) * U ) ) +  2 * U - Uold; 

    Uoldx = Ux;
    Ux = Unewx;

    Uoldy = Uy;
    Uy = Unewy;


    %Uold = U;
    %U = Unew;
    
    subplot( 1, 2, 1 );
    pcolor( X, Y, reshape( Ux, [NX, NY] ) );
    colorbar;
    vm = max( Ux );
    if vm ~= 0
        caxis( [-vm, vm] );
    end
    axis equal;
    subplot( 1, 2, 2 );
    pcolor( X, Y, reshape( Uy, [NX, NY] ) );
    %trisurf( T, P(:, 1), P(:, 2), zeros( Enodes, 1 ), U );
    %shading interp;
    %view(2);
    colorbar;
    vm = max( Uy );
    caxis( [-vm, vm] );
    axis equal;
    
    drawnow;
end
toc
%{
%}