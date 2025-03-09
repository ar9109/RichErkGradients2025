function [arr] = nan2zero(arr)
%NAN2ZERO Summary of this function goes here
%   Detailed explanation goes here
arr(isnan(arr)) = 0;

end

