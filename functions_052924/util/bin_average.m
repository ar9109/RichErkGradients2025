function [yMean,xMean,ySEM,yN] = bin_average(xData,yData,nbins,min_x,max_x)
% min_x and max_x must contain all datapoints !!
arguments
    xData
    yData
    nbins
    min_x = min(xData);
    max_x = max(xData);


end
%BINNING Summary of this function goes here
%   Detailed explanation goes here
removeNaN=~isnan(xData)&~isnan(yData);
xData = xData(removeNaN);xData = xData(:);
yData = yData(removeNaN);yData = yData(:);
xBins = linspace(min_x,max_x,nbins+1);xBins = xBins(:);
idxBins=discretize(xData,xBins);
idxBins = idxBins(~isnan(idxBins));
xMean = [(xBins(1:end-1)+xBins(2:end))/2];
yMean = accumarray(idxBins,yData,size(xMean),@nanmean,NaN);
yStd = accumarray(idxBins,yData,size(xMean),@nanstd,NaN);
yN = accumarray(idxBins,yData,size(xMean),@(x)sum(~isnan(x)),NaN);
ySEM = yStd./sqrt(yN);
end

