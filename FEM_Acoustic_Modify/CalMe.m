function Me = CalMe( r, s, w, T )
    Enodes = size( T, 2);
    Me = zeros( Enodes, Enodes );
    
    for j = 1 : Enodes
        for i = 1 : Enodes
            f = phi( r, s, i) .* phi( r, s, j);
            Me( i, j ) = GaussianQuadrature2D( f, r, s, w );
        end
    end
    
    %{Me = [1/12,1/24,1/24;1/24,1/12,1/24;1/24,1/24,1/12];%}

end