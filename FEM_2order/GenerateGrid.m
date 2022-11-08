function [X, Y] = GenerateGrid( NX, NY, DH, Nperiod, DefCoef )
    X = zeros( NX, NY );
    omiga = 2 * pi * Nperiod / ( NX * DH );
    for j = 1 : NY
        for i = 1 : NX
            X( i, j ) = DH * i;
            Y( i, j ) = DefCoef * DH * sin( i * DH * omiga ) + DH * j;
            %Y( i, j ) = DH * j;
        end
    end
end