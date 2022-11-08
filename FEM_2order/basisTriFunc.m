function  x = basisTriFunc( r, s, order)
    NP = length( r );
    A = zeros( NP, NP );
    for j = 1 : NP 
        A(j, 1) = 1; 
        A(j, 2) = r(j); 
        A(j, 3) = s(j); 
        A(j, 4) = r(j)*r(j);
        A(j, 5) = r(j)*s(j);
        A(j, 6) = s(j)*s(j);
    end

    B = eye( NP );
    x = inv( A ) * eye;

    
end


