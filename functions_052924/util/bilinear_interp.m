function [fxy] = bilinear_interp(x,y,F,x1,x2,y1,y2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fxy = 1./(x2-x1)./(y2-y1)*[x2-x x-x1]*[F(x1,y1) F(x1,y2); F(x2,y1) F(x2,y2)]*[y2-y;y-y1];
end

