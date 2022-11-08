function plotGrid( X, Y, sample, NX, NY )
    for j = 1 : sample : NY
        plot( X( :, j ), Y( :, j ), 'k' );
        hold on;
    end
    for i = 1 : sample : NX
        plot( X( i, : ), Y( i, : ), 'k' );
        hold on;
    end
end
