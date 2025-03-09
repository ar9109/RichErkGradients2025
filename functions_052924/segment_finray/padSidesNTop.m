function [stack1,stack2] = padSidesNTop(stack1,stack2)
%PADSIDESNTOP Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(stack1) && ~isempty(stack2)
    sizediff = size(stack1,[1,2]) - size(stack2,[1,2]);
    % pad rows
    if sizediff(1) < 0
        stack1 = padarray(stack1,abs(sizediff(1)),0,'pre');
    elseif sizediff(1) > 0
        stack2 = padarray(stack2,abs(sizediff(1)),0,'pre');
    end
    % pad columns
    if sizediff(2) < 0
        stack1 = padarray(stack1,[0,floor(abs(sizediff(2))./2)],0,'both');
        if mod(abs(sizediff(2)),2)
            stack1 = padarray(stack1,[0,1],0,'post');
        end
    elseif sizediff(2) > 0
        stack2 = padarray(stack2,[0,floor(abs(sizediff(2))./2)],0,'both');
        if mod(abs(sizediff(2)),2)
            stack2 = padarray(stack2,[0,1],0,'post');
        end
    end
else
%     error("one of the input stacks is empty")

end

% size(stack1)
% size(stack2)
