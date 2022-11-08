function M = AssembleMassMatrix( Me, JField, P, T )
    Nodes = size( P, 1 );
    Enodes = size( T, 2 );
    M = zeros( Nodes, Nodes );
    Elements = size( T, 1 );
    for i = 1 : Elements
        J = JField( i );
        for ln1 = 1 : Enodes
            for ln2 = 1 : Enodes
                gn1 = T( i, ln1 );
                gn2 = T( i, ln2 );
                M( gn1, gn2 ) = M( gn1, gn2 ) + abs( J  ) * Me( ln1, ln2 );
            end
        end
    end
    
end