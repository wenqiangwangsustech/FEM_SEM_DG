function [Se, Se1, Se2] = CalSe( r, s, w, r_x, r_y, s_x, s_y )
    Np = length(r);    
    Se = zeros( 3, 3 );
    Se1 = zeros( 3, 3 );
    Se2 = zeros( 3, 3 );

    f11 = ( derivPhi( r, s, 1, 'R' ) * r_x + derivPhi( r, s, 1, 'S' ) * s_x ) .* ( derivPhi( r, s, 1, 'R' ) * r_x + derivPhi( r, s, 1, 'S' ) * s_x ); 
    f12 = ( derivPhi( r, s, 1, 'R' ) * r_x + derivPhi( r, s, 1, 'S' ) * s_x ) .* ( derivPhi( r, s, 2, 'R' ) * r_x + derivPhi( r, s, 2, 'S' ) * s_x ); 
    f13 = ( derivPhi( r, s, 1, 'R' ) * r_x + derivPhi( r, s, 1, 'S' ) * s_x ) .* ( derivPhi( r, s, 3, 'R' ) * r_x + derivPhi( r, s, 3, 'S' ) * s_x ); 
    f22 = ( derivPhi( r, s, 2, 'R' ) * r_x + derivPhi( r, s, 2, 'S' ) * s_x ) .* ( derivPhi( r, s, 2, 'R' ) * r_x + derivPhi( r, s, 2, 'S' ) * s_x ); 
    f23 = ( derivPhi( r, s, 2, 'R' ) * r_x + derivPhi( r, s, 2, 'S' ) * s_x ) .* ( derivPhi( r, s, 3, 'R' ) * r_x + derivPhi( r, s, 3, 'S' ) * s_x ); 
    f33 = ( derivPhi( r, s, 3, 'R' ) * r_x + derivPhi( r, s, 3, 'S' ) * s_x ) .* ( derivPhi( r, s, 3, 'R' ) * r_x + derivPhi( r, s, 3, 'S' ) * s_x ); 

    Se1( 1, 1 ) = GaussianQuadrature2D( f11, r, s, w );
    Se1( 1, 2 ) = GaussianQuadrature2D( f12, r, s, w );
    Se1( 1, 3 ) = GaussianQuadrature2D( f13, r, s, w );
    Se1( 2, 2 ) = GaussianQuadrature2D( f22, r, s, w );
    Se1( 2, 3 ) = GaussianQuadrature2D( f23, r, s, w );
    Se1( 3, 3 ) = GaussianQuadrature2D( f33, r, s, w );

    Se1( 2, 1 ) = Se1( 1 , 2 );
    Se1( 3, 1 ) = Se1( 1 , 3 );
    Se1( 3, 2 ) = Se1( 2 , 3 );

    f11 = ( derivPhi( r, s, 1, 'R' ) * r_y + derivPhi( r, s, 1, 'S' ) * s_y ) .* ( derivPhi( r, s, 1, 'R' ) * r_y + derivPhi( r, s, 1, 'S' ) * s_y ); 
    f12 = ( derivPhi( r, s, 1, 'R' ) * r_y + derivPhi( r, s, 1, 'S' ) * s_y ) .* ( derivPhi( r, s, 2, 'R' ) * r_y + derivPhi( r, s, 2, 'S' ) * s_y ); 
    f13 = ( derivPhi( r, s, 1, 'R' ) * r_y + derivPhi( r, s, 1, 'S' ) * s_y ) .* ( derivPhi( r, s, 3, 'R' ) * r_y + derivPhi( r, s, 3, 'S' ) * s_y ); 
    f22 = ( derivPhi( r, s, 2, 'R' ) * r_y + derivPhi( r, s, 2, 'S' ) * s_y ) .* ( derivPhi( r, s, 2, 'R' ) * r_y + derivPhi( r, s, 2, 'S' ) * s_y ); 
    f23 = ( derivPhi( r, s, 2, 'R' ) * r_y + derivPhi( r, s, 2, 'S' ) * s_y ) .* ( derivPhi( r, s, 3, 'R' ) * r_y + derivPhi( r, s, 3, 'S' ) * s_y ); 
    f33 = ( derivPhi( r, s, 3, 'R' ) * r_y + derivPhi( r, s, 3, 'S' ) * s_y ) .* ( derivPhi( r, s, 3, 'R' ) * r_y + derivPhi( r, s, 3, 'S' ) * s_y ); 

    Se2( 1, 1 ) = GaussianQuadrature2D( f11, r, s, w );
    Se2( 1, 2 ) = GaussianQuadrature2D( f12, r, s, w );
    Se2( 1, 3 ) = GaussianQuadrature2D( f13, r, s, w );
    Se2( 2, 2 ) = GaussianQuadrature2D( f22, r, s, w );
    Se2( 2, 3 ) = GaussianQuadrature2D( f23, r, s, w );
    Se2( 3, 3 ) = GaussianQuadrature2D( f33, r, s, w );

    Se2( 2, 1 ) = Se2( 1 , 2 );
    Se2( 3, 1 ) = Se2( 1 , 3 );
    Se2( 3, 2 ) = Se2( 2 , 3 );

    Se = Se1 + Se2;

end