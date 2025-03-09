function [arr_out,notnanidx] = remove_nan(arr,col)
arguments
    arr
    col=0; % remove for all columns
end
%REMOVE_NAN Summary of this function goes here
%   Detailed explanation goes here
if isvector(arr)
    notnanidx = ~isnan(arr);
    arr_out = arr(notnanidx);
elseif col == 0
    notnanidx = all(~isnan(arr),2);
    arr_out = arr(notnanidx,:);
else
    notnanidx = ~isnan(arr(:,col));
    arr_out = arr(notnanidx,:);
end

end

