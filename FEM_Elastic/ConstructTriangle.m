function [P, T] = ConstructTriangle( X, Y, NX, NY )
    N = ( NX - 1 ) * ( NY - 1 ) * 2;
    P = zeros( NX * NY, 2 );
    T = zeros( N, 3 );
    for j = 1 : NY - 1
        for i = 1 : NX - 1
            J = j - 1;

            e = ( i + J*(NX - 1)) * 2;
            
            T(e - 1, 1) = i    +J    *NX;
            T(e - 1, 2) = (i+1)+J    *NX;
            T(e - 1, 3) = i    +(J+1)*NX;
            
            P( T(e - 1, 1), 1 ) = X(i,j);
            P( T(e - 1, 2), 1 ) = X(i+1,j);
            P( T(e - 1, 3), 1 ) = X(i,j+1);
            

            P( T(e - 1, 1), 2 ) = Y(i,j);
            P( T(e - 1, 2), 2 ) = Y(i+1,j);
            P( T(e - 1, 3), 2 ) = Y(i,j+1);
            
            T(e, 1) = (i+1)+(J+1)*NX;
            T(e, 2) = i    +(J+1)*NX;
            T(e, 3) = (i+1)+J    *NX;

            P( T(e, 1), 1 ) = X(i+1, j+1);
            P( T(e, 2), 1 ) = X(i,   j+1);
            P( T(e, 3), 1 ) = X(i+1, j  );

            P( T(e, 1), 2 ) = Y(i+1, j+1);
            P( T(e, 2), 2 ) = Y(i,   j+1);
            P( T(e, 3), 2 ) = Y(i+1, j  );

        end
    end

end