function JField = CalJacobian( P, T )
    N = size( T, 1 );
    JField = zeros( N, 1 );
    for i = 1 : N
        x1 = P( T(i, 1), 1 ); y1 = P( T(i, 1), 2 );
        x2 = P( T(i, 2), 1 ); y2 = P( T(i, 2), 2 );
        x3 = P( T(i, 3), 1 ); y3 = P( T(i, 3), 2 );

        J = (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1);
        JField( i ) = J;
    end
end