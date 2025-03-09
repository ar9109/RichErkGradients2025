function [varargout] = expand_cell(cellin)
%EXPAND_CELL Summary of this function goes here
%   Detailed explanation goes here
varargout = {cellin{:}};
if nargout~=numel(varargout)
    error('size doesn''t match');
end
end

