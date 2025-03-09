function dydt = ode_match_eigvec(t,Y,e2)
%MATCH_EIGVEC Summary of this function goes here
%   Detailed explanation goes here 
    x = Y(1); y = Y(2);
    if any(round(Y)<=1) || any(round(Y)>=size(e2,[2,1])')
        dx = 0;
        dy = 0;
    else
%         dx = bilinear_interp(x,y,e2(:,:,2)',floor(x),ceil(x),floor(y),ceil(y));
%         dy = bilinear_interp(x,y,e2(:,:,1)',floor(x),ceil(x),floor(y),ceil(y));
        dx = e2(round(y),round(x),2);
        dy = e2(round(y),round(x),1);
    end
    dydt = [dx;dy];
end

