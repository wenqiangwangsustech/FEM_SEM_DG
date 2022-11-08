function S = AssembleStiffMatrix( JField, P, T, r, s, w )
    Nodes = length( P );
    Enodes = size( T, 2 );
    S = sparse( Nodes, Nodes );
    %S1 = zeros( Nodes, Nodes );
    %S2 = zeros( Nodes, Nodes );
    Elements = size( T, 1 );
    for i = 1 : Elements
        J = JField( i );
        x1 = P( T(i, 1), 1 ); y1 = P( T(i, 1), 2 );
        x2 = P( T(i, 2), 1 ); y2 = P( T(i, 2), 2 );
        x3 = P( T(i, 3), 1 ); y3 = P( T(i, 3), 2 );
               
        x_r = x2 - x1;
        x_s = x3 - x1;
        y_r = y2 - y1;
        y_s = y3 - y1;

        r_x = y_s / J;
        s_x = - y_r / J;
        r_y = - x_s / J;
        s_y = x_r / J;

        [Se, Se1, Se2] = CalSe( r, s, w, r_x, r_y, s_x, s_y );
        for ln1 = 1 : Enodes
            for ln2 = 1 : Enodes
                gn1 = T( i, ln1 );
                gn2 = T( i, ln2 );
                S( gn1, gn2 ) = S( gn1, gn2 ) + abs( J ) * Se( ln1, ln2 );
                %S1( gn1, gn2 ) = S1( gn1, gn2 ) + J * Se1( ln1, ln2 );
                %S2( gn1, gn2 ) = S2( gn1, gn2 ) + J * Se2( ln1, ln2 );
            end
        end
    end
   %S = S1 + S2;
    
end