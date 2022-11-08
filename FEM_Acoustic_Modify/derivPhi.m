function v = derivPhi( r, s, flag, Direction )
    tmp = ones( size(r) );
    if Direction == 'R'
        if flag == 1
            v = - 1 .* tmp;
        end
            
        if flag == 2
            v = 1 .* tmp;
        end
         
        if flag == 3
            v = 0 .* tmp;
        end
    end
    if Direction == 'S'
        if flag == 1
            v = - 1 .* tmp;
        end

        if flag == 2
            v = 0 .* tmp;
        end

        if flag == 3
            v = 1 .* tmp;
        end
    end
    %v
end 