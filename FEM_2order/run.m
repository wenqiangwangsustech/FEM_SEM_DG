close all;
clear all;

xp = [0, 1, 0, 0.5, 0.5, 0];
yp = [0, 0, 1, 0, 0.5, 0.5];

order = 2;
coef  = basisTriFunc( xp, yp, order );
NP = (order + 1) * (order + 2) / 2;

NX = 100;
NY = 100;
DH = 200;
Nperiod = 2;
DefCoef = 2;
[X, Y] = GenerateGrid( NX, NY, DH, Nperiod, DefCoef );
[P, T] = ConstructTriangle( X, Y, NX, NY, order );

nx = order * (NX - 1) + 1;
ny = order * (NY - 1) + 1;

%plotTriangle( P, T );
%网格生成
% load homo_allpml.mat
% P = p;
% T = t;


NT = 500;

src = zeros( NT, 1 );


Nodes = size( P, 1 );
Nodes

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
Indexes = WaveIndex( nx, ny, NX, NY );

Me = zeros( Enodes, Enodes );
Me = CalMe( coef, r, s, w, T );

display( "=======雅可比系数求解=======" );
JField = CalJacobian( P, T );
display( "=======质量矩阵组装=======" );
M = AssembleMassMatrix ( Me, JField, Nodes, T );
display( "=======刚度矩阵组装=======" );
[Sxx, Sxy, Syx, Syy] = AssembleStiffMatrix( coef, JField, Nodes, P, T, r, s, w);

DT = 0.005;

Vp = 5000;
Vs = 0.5774*Vp;

f0 = 2;
mu = Vs * Vs;
lam = Vp * Vp - 2.0 * Vs * Vs;
lam2mu = Vp * Vp;

display( "=======矩阵求逆=======" );

[L, U, P, Q, D] = lu(M);

%Minv = DT ^ 2 * inv( M );
srcX = NX;
srcY = NY;


tic

F1 = zeros( Nodes, 1 );
F2 = zeros( Nodes, 1 );
for it = 1 : NT
    display( it );
    
    src( it ) = ( 1 - 2 * ( pi * f0 * ( it * DT - 1 / f0) ) ^ 2 ) * exp ( - ( pi * f0 * ( it * DT - 1 / f0 ) ) ^ 2 );
    
    F2( srcX + srcY * nx, 1 ) = src( it );
    %{
    Unewx = Minv * ( F1 - ( lam2mu * Sxx * Ux + lam * Syx * Uy + mu * Sxy * Uy + mu * Syy * Ux ) ) +  2 * Ux - Uoldx; 
    Unewy = Minv * ( F2 - ( mu * Sxx * Uy + mu * Syx * Ux + lam * Sxy * Ux + lam2mu * Syy * Uy ) ) +  2 * Uy - Uoldy;
    %}
    Tmp1 = F1 - ( lam2mu * Sxx * Ux + lam * Syx * Uy + mu * Sxy * Uy + mu * Syy * Ux );
    Tmp2 = F2 - ( mu * Sxx * Uy + mu * Syx * Ux + lam * Sxy * Ux + lam2mu * Syy * Uy );
    Unewx = Q * ( U \ ( L \ ( P * ( D \ Tmp1 ) ) ) ) * DT^2 +  2 * Ux - Uoldx; 
    Unewy = Q * ( U \ ( L \ ( P * ( D \ Tmp2 ) ) ) ) * DT^2 +  2 * Uy - Uoldy; 
    
    
    Uoldx = Ux;
    Ux = Unewx;

    Uoldy = Uy;
    Uy = Unewy;

    subplot( 1, 2, 1 );
    pcolor( X, Y, reshape( Ux(Indexes), [NX, NY] ) );
    shading interp;
    view(2);
    colorbar;
    vm = max( Ux );
    if vm ~= 0
        caxis( [-vm, vm] );
    end
    axis equal;
    subplot( 1, 2, 2 );
    pcolor( X, Y, reshape( Uy(Indexes), [NX, NY] ) );
    %trisurf( T, P(:, 1), P(:, 2), zeros( Enodes, 1 ), U );
    shading interp;
    view(2);
    colorbar;
    vm = max( Uy );
    caxis( [-vm, vm] );
    axis equal;
    
    
%{
    Unew = Minv * ( F2 - Vp ^ 2 * ( Sxx + Syy ) * U ) +  2 * U - Uold; 
    Uold = U;
    U = Unew;
    
    pcolor( reshape( U, [nx, ny] ) );
    %trisurf( T, P(:, 1), P(:, 2), zeros( Enodes, 1 ), U );
    %shading interp;
    %view(2);
    colorbar;
    vm = max( U );
    if vm ~= 0
        caxis( [-vm, vm] );
    end
    axis equal;
%}    
    drawnow;
end
toc
%{%}