function [der,f] = find_derivative(xdata,ydata,xinput,smoothParam)
arguments
    xdata;
    ydata;
    xinput = xdata;
    smoothParam = [];
end
%FIND_DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
if numel(xdata)>1
    if isempty(smoothParam)
        [f] = fit(xdata(:),ydata(:),'smoothingspline');
    else
        [f] = fit(xdata(:),ydata(:),'smoothingspline','SmoothingParam',smoothParam);
    end
    der = ppval(fnder(f.p),xinput);
else
    der = nan;
end

