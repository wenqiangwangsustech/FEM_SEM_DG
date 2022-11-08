function [ Sxx, Sxy, Syx, Syy ] = AssembleStiffMatrix( JField, P, T, r, s, w )
    Nodes = length( P );
    Enodes = size( T, 2 );
    Sxx= zeros( Nodes, Nodes );
    Sxy= zeros( Nodes, Nodes );
    Syx= zeros( Nodes, Nodes );
    Syy= zeros( Nodes, Nodes );
    Elements = size( T, 1 );
    
    Ngp = length( r );
    phiR = zeros( Ngp, Enodes );
    phiS = zeros( Ngp, Enodes );
    
    for in = 1 : Enodes
        phiR( :, in ) = derivPhi( r, s, in, 'R' );
        phiS( :, in ) = derivPhi( r, s, in, 'S' );
    end

    %Triangle make r_x r_y s_x s_y constant.
    %Therefor, the Jacobian Inverse Matrix can move out of the integral.

    Serr = zeros( Enodes, Enodes );
    Sers = zeros( Enodes, Enodes );
    Sesr = zeros( Enodes, Enodes );
    Sess = zeros( Enodes, Enodes );

    for j = 1 : Enodes
        for i = 1 : Enodes
            Frr = phiR( :, i ) .* phiR( :, j );
            Frs = phiR( :, i ) .* phiS( :, j );
            Fsr = phiS( :, i ) .* phiR( :, j );
            Fss = phiS( :, i ) .* phiS( :, j );
            Serr( i, j ) = GaussianQuadrature2D( Frr, r, s, w );
            Sers( i, j ) = GaussianQuadrature2D( Frs, r, s, w );
            Sesr( i, j ) = GaussianQuadrature2D( Fsr, r, s, w );
            Sess( i, j ) = GaussianQuadrature2D( Fss, r, s, w );
        end
    end


    Sexx = zeros( Enodes, Enodes );
    Sexy = zeros( Enodes, Enodes );
    Seyx = zeros( Enodes, Enodes );
    Seyy = zeros( Enodes, Enodes );

    for i = 1 : Elements
        J = JField( i );
        x1 = P( T(i, 1), 1 ); y1 = P( T(i, 1), 2 );
        x2 = P( T(i, 2), 1 ); y2 = P( T(i, 2), 2 );
        x3 = P( T(i, 3), 1 ); y3 = P( T(i, 3), 2 );
               
        x_r = x2 - x1;
        x_s = x3 - x1;
        y_r = y2 - y1;
        y_s = y3 - y1;

        r_x =   y_s / J;
        s_x = - y_r / J;
        r_y = - x_s / J;
        s_y =   x_r / J;

        Sexx = r_x * r_x * Serr + r_x * s_x * Sers + s_x * r_x * Sesr + s_x * s_x * Sess; 
        Sexy = r_x * r_y * Serr + r_x * s_y * Sers + s_x * r_y * Sesr + s_x * s_y * Sess; 
        Seyx = r_y * r_x * Serr + r_y * s_x * Sers + s_y * r_x * Sesr + s_y * s_x * Sess; 
        Seyy = r_y * r_y * Serr + r_y * s_y * Sers + s_y * r_y * Sesr + s_y * s_y * Sess;        


        for ln1 = 1 : Enodes
            for ln2 = 1 : Enodes
                gn1 = T( i, ln1 );
                gn2 = T( i, ln2 );
                %S( gn1, gn2 ) = S( gn1, gn2 ) + abs( J ) * Se( ln1, ln2 );
                Sxx( gn1, gn2 ) = Sxx( gn1, gn2 ) + abs( J ) * Sexx( ln1, ln2 );
                Sxy( gn1, gn2 ) = Sxy( gn1, gn2 ) + abs( J ) * Sexy( ln1, ln2 );
                Syx( gn1, gn2 ) = Syx( gn1, gn2 ) + abs( J ) * Seyx( ln1, ln2 );
                Syy( gn1, gn2 ) = Syy( gn1, gn2 ) + abs( J ) * Seyy( ln1, ln2 );
            end
        end
    end
    
end