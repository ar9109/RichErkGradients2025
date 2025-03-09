function [idx] = find_from_plot(f,line_option)
arguments
    f
    line_option = [];
end
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
datacursormode on
dcm_obj = datacursormode(f);
cursor_info = getCursorInfo(dcm_obj);

if isempty(line_option)
    idx = find(all(cursor_info.Position==[cursor_info.Target.XData' cursor_info.Target.YData'],2));
elseif strcmp(line_option, 'line')
    idx = numel(cursor_info.Target.Parent.Children)-find(cursor_info.Target.Parent.Children==cursor_info.Target)+1;
else
    error('line_option has to be either empty or ''line''')
end
end

