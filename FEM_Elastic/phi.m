function v = phi( r, s, flag )
    v = zeros( length( r ), 1 );
    if flag == 1
        v = - r - s + 1;
    end
        
    if flag == 2
        v = r;
    end
     
    if flag == 3
        v = s;
    end
end 