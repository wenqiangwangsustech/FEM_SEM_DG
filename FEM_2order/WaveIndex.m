function Indexes = WaveIndex( nx, ny, NX, NY )
    Indexes = zeros( NX * NY, 1 );
    n = 1;
    for j = 1 : 2 : ny
        for i = 1 : 2 : nx
            Indexes(n) = i + ( j - 1 ) * nx;
            n = n + 1;
        end
    end
end