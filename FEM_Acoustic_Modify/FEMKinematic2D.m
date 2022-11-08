clear all;
clc;

Nt = 1000;
Nx =100;
Nz =100;
grid = Nx * Nz;
v = 3000;
srcLocationX  = 50;
srcLocationZ  = 50;
src = zeros( Nt, 1 );
density = 2500;

V = v * ones(grid, 1);
h = 10;

deltaT = 0.0005;
f0 = 10;

Unew = zeros( grid, 1);
Uold = zeros( grid, 1);
U = zeros( grid, 1);

F = zeros(grid, 1);
M = zeros( grid, grid );
K = zeros( grid, grid );
Minv_K = zeros( grid, grid );


%%                                                                                     
for i = 1 : ( Nx*( Nz - 1 ) - 1 )
    if mod( i, Nx )  ~= 0 %&& mod( i - 1, Nx )  ~= 0
        M( i,  i)  =M( i,  i) +  h^2 / 9;
        M( i,  i+1)  = M( i,  i+1) + h^2 / 18;
        M( i,  i+Nx)  =  M( i,  i+Nx) + h^2 / 18;
        M( i,  i+1+Nx)  = M( i,  i+1+Nx) + h^2 / 36;
        
        M( i+1,  i)  = M( i+1,  i)  + h^2 / 18;
        M( i+1,  i+1)  = M( i+1,  i+1) + h^2 / 9;
        M( i+1,  i+Nx)  = M( i+1,  i+Nx) + h^2 / 36;
        M( i+1,  i+1+Nx)  = M( i+1,  i+1+Nx) + h^2 / 18;
                
        M( i+Nx,  i)  =M( i+Nx,  i) +  h^2 / 18;
        M( i+Nx,  i+1)  = M( i+Nx,  i+1) + h^2 / 36;
        M( i+Nx,  i+Nx)  =M( i+Nx,  i+Nx)  + h^2 / 9;
        M( i+Nx,  i+1+Nx)  = M( i+Nx,  i+1+Nx)  + h^2 / 18;
        
        M( i+1+Nx,  i)  = M( i+1+Nx,  i) + h^2 / 36;
        M( i+1+Nx,  i+1)  = M( i+1+Nx,  i+1) + h^2 / 18;
        M( i+1+Nx,  i+Nx)  = M( i+1+Nx,  i+Nx) + h^2 / 18;
        M( i+1+Nx,  i+1+Nx)  =M( i+1+Nx,  i+1+Nx) +  h^2 / 9;       
        
        K( i,  i)  =K( i,  i) +  2 / 3;
        K( i,  i+1)  = K( i,  i+1) - 1 / 6;
        K( i,  i+Nx)  =  K( i,  i+Nx)  - 1 / 6;
        K( i,  i+1+Nx)  = K( i,  i+1+Nx) - 1 / 3;
        
        K( i+1,  i)  = K( i+1,  i)   - 1 / 6;
        K( i+1,  i+1)  = K( i+1,  i+1) + 2 / 3;
        K( i+1,  i+Nx)  = K( i+1,  i+Nx) - 1 / 3;
        K( i+1,  i+1+Nx)  = K( i+1,  i+1+Nx)  - 1 / 6;
                
        K( i+Nx,  i)  =K( i+Nx,  i)  - 1 / 6;
        K( i+Nx,  i+1)  = K( i+Nx,  i+1) - 1 / 3;
        K( i+Nx,  i+Nx)  =K( i+Nx,  i+Nx)  + 2 / 3;
        K( i+Nx,  i+1+Nx)  = K( i+Nx,  i+1+Nx)   - 1 / 6;
        
        K( i+1+Nx,  i)  = K( i+1+Nx,  i) - 1 / 3;
        K( i+1+Nx,  i+1)  = K( i+1+Nx,  i+1)  - 1 / 6;
        K( i+1+Nx,  i+Nx)  = K( i+1+Nx,  i+Nx)  - 1 / 6;
        K( i+1+Nx,  i+1+Nx)  =K( i+1+Nx,  i+1+Nx) +  2 / 3;
    end
end
%%


%%
Minv =deltaT^2 * inv( density * M );
K = v^2 *density *K;
%Minv_K = Minv* K;
%%
%%


tic
for it = 1 : Nt
    display( it );
    src( it ) = ( 1 - 2 * ( pi * f0 * ( it * deltaT - 1 / f0) ) ^ 2 ) * exp ( - ( pi * f0 * ( it * deltaT - 1 / f0 ) ) ^ 2 );
 %%    
%     U( srcLocationX + srcLocationZ * Nx ) = U( srcLocationX + srcLocationZ * Nx )+ deltaT^2 * src( it );
%     Unew = - Minv_K * U  +  2 * U - Uold;
    
%%     
    F( srcLocationX+ srcLocationZ * Nx ) = h^4 * src( it );
    Unew = Minv * ( F - K* U ) +  2 * U - Uold;

%%
    %Unew = Minv_F - Minv_K* U +  2 * U - Uold;%%This is error, because F  varies with time.
%%  
    
%%
    Uold = U;
    U = Unew;
    
    drawGif( reshape( U,Nx , Nz ), it, deltaT, 'FEM' );
    drawnow;
end
toc