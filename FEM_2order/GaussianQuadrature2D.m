function GaussIntegral = GaussianQuadrature2D( funcV, r, s, w )
    Np = length(r);
    sum0 = 0;
    for i = 1 : Np
        sum0 = sum0 + funcV(i) * w( i );
    end
    GaussIntegral = sum0;
end