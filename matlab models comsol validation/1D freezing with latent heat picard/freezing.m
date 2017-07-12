% freezing function

function rhophrate_i = freezing(T,phi_i(ie))

        if phi_i(ie) == 0 && T > T_m ;
            rhophrate_i = 0; % when there is no ice, no more melting
        elseif phi - phi_i(ie) == 0 && T < T_m ;
            rhophrate_i = 0; % when there is no water, no more freezing
        else 
            rhophrate_i = alpha*(T - T_m) ;
        end