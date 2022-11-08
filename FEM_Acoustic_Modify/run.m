close all;
clear all;

NX = 100;
NY = 100;
DH = 100;
Nperiod = 4;
DefCoef = 0;
[X, Y] = GenerateGrid( NX, NY, DH, Nperiod, DefCoef );
[P, T] = ConstructTriangle( X, Y, NX, NY );
%plotTriangle( P, T );
%网格生成
% load homo_allpml.mat
% P = p;
% T = t;

NT = 200;

DT = 0.02;
f0 = 2;


v = 3000;

src = zeros( NT, 1 );


Nodes = size( P, 1 );
Unew  = zeros( Nodes, 1 );
Uold  = zeros( Nodes, 1 );
U     = zeros( Nodes, 1 );


display( "=======产生高斯点=======" );
n = 3;
[r, s, w] = genGaussPoint( n );
Enodes = size( T, 2 );


Me = zeros( Enodes, Enodes );
Me = CalMe( r, s, w, T );

display( "=======雅可比系数求解=======" );
JField = CalJacobian( P, T );
display( "=======质量矩阵组装=======" );
M = AssembleMassMatrix ( Me, JField, P, T );
display( "=======刚度矩阵组装=======" );
S= AssembleStiffMatrix( JField, P, T, r, s, w );
%[S1,S2,M] = FEM_MA(P,T);
%S = S1 + S2;



display( "=======矩阵求逆=======" );

Minv = DT ^ 2 * inv( M );
S = v^2 * S;

srcX = NX/2;
srcY = NY/2;


tic

F = zeros( Nodes, 1 );
for it = 1 : NT
    display( it );
    src( it ) = ( 1 - 2 * ( pi * f0 * ( it * DT - 1 / f0) ) ^ 2 ) * exp ( - ( pi * f0 * ( it * DT - 1 / f0 ) ) ^ 2 );
    F( srcX + srcY * NX,1 ) =  src( it );
    Unew = Minv * ( F - S * U ) +  2 * U - Uold; 
    
    Uold = U;
    U = Unew;
    
    pcolor( X, Y, reshape( U, [NX, NY] ) );
    %trisurf( T, P(:, 1), P(:, 2), zeros( Enodes, 1 ), U );
    shading interp;
    view(2);
    colorbar;
    caxis( [-1e-9, 5e-9] );
    axis equal;
    
    drawnow;
end
toc
%{
%}