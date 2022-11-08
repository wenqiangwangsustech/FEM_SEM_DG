function v = derivPhi( coef, r, s, flag, Direction )
    v = zeros( length( r ), 1 );
    if Direction == 'R'
        v = coef( 2, flag) + ...
            2 * coef( 4, flag) * r + ...
            coef( 5, flag) * s;
    end
    if Direction == 'S'
        v = coef( 3, flag) + ...
            coef( 5, flag) * r + ...
            2 * coef( 6, flag) * s;
    end
    %v
end 