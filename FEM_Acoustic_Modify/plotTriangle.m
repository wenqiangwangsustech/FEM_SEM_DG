function plotTriangle( P, T )
    Ne = size( T );
    x = zeros( 4, 1 );
    y = zeros( 4, 1 );
    for i = 1 : Ne
        for j = 1 : 3
            x(j) = P( T( i, j ), 1 );
            y(j) = P( T( i, j ), 2 );
        end
        x(4) = P( T( i, 1 ), 1 );
        y(4) = P( T( i, 1 ), 2 );
        line(x, y);
        hold on;
    end
    axis equal;
end