function [arr_out,notnanidx] = remove_nan_inf(arr,col)
%REMOVE_NAN_INF Summary of this function goes here
%   Detailed explanation goes here
arguments
    arr
    col=0; % remove for all columns
end

% x1 = x(:);
% x_out = x(~isnan(x)&~isinf(x));

if isvector(arr)
    notnanidx = ~isnan(arr)&~isinf(arr);
    arr_out = arr(notnanidx);
elseif col == 0
    notnanidx = all(~isnan(arr)&~isinf(arr),2);
    arr_out = arr(notnanidx,:);
else
    notnanidx = ~isnan(arr(:,col))&~isinf(arr(:,col));
    arr_out = arr(notnanidx,:);
end

end


