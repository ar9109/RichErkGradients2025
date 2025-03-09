function [value, isterminal, direction] = odeEvent_stop_at_boundary(T, Y, mask)
    value      = mask(round(Y(2)),round(Y(1)));
    isterminal = 1;   % Stop the integration
    direction  = 0;
end

