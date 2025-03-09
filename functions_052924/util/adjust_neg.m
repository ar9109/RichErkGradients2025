function [arr] = adjust_neg(arr)
%ADJUST_NEG Summary of this function goes here
%   Detailed explanation goes here
arr(arr<0) = 0;
end

