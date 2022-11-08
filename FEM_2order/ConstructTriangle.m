function [P, T] = ConstructTriangle( X, Y, NX, NY, order )
    N = ( NX - 1 ) * ( NY - 1 ) * 2;
    NP = ( order + 1 ) * ( order + 2 ) / 2;
    nx = order * (NX - 1) + 1;
    ny = order * (NY - 1) + 1;
    P = zeros( nx * ny, 2 );
    T = zeros( N, NP ) + 1;
    for j = 1 : NY - 1
        for i = 1 : NX - 1
            J = j - 1;

            e = ( i + J*(NX - 1)) * 2;
            
            T(e - 1, 1) = 2 * i - 1 + 2 * J * nx;
            T(e - 1, 2) = 2 * i + 1 + 2 * J * nx;
            T(e - 1, 3) = 2 * i - 1 + 2 * (J + 1) * nx;
            
            T(e - 1, 4) = 2 * i + 2 * J * nx;
            T(e - 1, 5) = 2 * i + ( 2 * J + 1 ) * nx;
            T(e - 1, 6) = 2 * i - 1 + ( 2 * J + 1 ) * nx;
            
            
            P( T(e - 1, 1), 1 ) = X(i,j);
            P( T(e - 1, 2), 1 ) = X(i+1,j);
            P( T(e - 1, 3), 1 ) = X(i,j+1);

            P( T(e - 1, 1), 2 ) = Y(i,j);
            P( T(e - 1, 2), 2 ) = Y(i+1,j);
            P( T(e - 1, 3), 2 ) = Y(i,j+1);
            
            T(e, 1) = 2 * i + 1 + 2 * (J + 1) * nx;
            T(e, 2) = 2 * i - 1 + 2 * (J + 1) * nx;
            T(e, 3) = 2 * i + 1 + 2 * J * nx;

            T(e, 4) = 2 * i + 2 * (J + 1) * nx;
            T(e, 5) = 2 * i + ( 2 * J + 1 ) * nx;
            T(e, 6) = 2 * i + 1 + ( 2 * J + 1 ) * nx;

            P( T(e, 1), 1 ) = X(i+1, j+1);
            P( T(e, 2), 1 ) = X(i,   j+1);
            P( T(e, 3), 1 ) = X(i+1, j  );

            P( T(e, 1), 2 ) = Y(i+1, j+1);
            P( T(e, 2), 2 ) = Y(i,   j+1);
            P( T(e, 3), 2 ) = Y(i+1, j  );

        end
    end

end