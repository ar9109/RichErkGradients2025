function [single_exp] = generate_oneray(single_exp,ktr,gem,ccrot,ampPointrot,endPointrot,dbins, gem_thres)
%GENERATE_ONERAY Summary of this function goes here
%   Detailed explanation goes here



% binning
%     dbins = 50;
%     bins = linspace(0,(s.endPointrot(1)-s.ampPointrot(1)),21);
bins = 0:dbins:(endPointrot(1)-ampPointrot(1));
%     bins = flip((endPointrot(1)-ampPointrot(1)):-dbins:0);

idxBin = discretize(ccrot(:,1),bins); % assigns each column to a particular bin
% Initialise
% KTR
averageKTR = nan(1,numel(bins)-1);
semKTR = nan(1,numel(bins)-1);
averageKTR = [];
semKTR = [];

% GEM
fractionGEM = nan(1,numel(bins)-1);
nucleiGEM = nan(1,numel(bins)-1);
nucleiALL = nan(1,numel(bins)-1);
fractionGEM = [];
nucleiGEM = [];
nucleiALL = [];
averageGEM = [];
medianGEM = [];

% binvalue
binvalue = nan(1,numel(bins)-1);
binvalue = [];

for bin = 1:(numel(bins)-1)
    % KTR
    if ~isempty(ktr)
        averageKTR(bin) = mean(ktr(idxBin==bin));
        semKTR(bin) = std(ktr(idxBin==bin))/sqrt(sum(idxBin==bin));
    end
    % GEM
    if ~isempty(gem)
        nucleiGEM(bin) = sum(idxBin==bin&gem>gem_thres);
        nucleiALL(bin) = sum(idxBin==bin);
        fractionGEM(bin) = nucleiGEM(bin)./nucleiALL(bin);
        averageGEM(bin) = mean(gem(idxBin==bin));
        medianGEM(bin) = median(gem(idxBin==bin));
    end
    % binvalue
    binvalue(bin) = (bins(bin)+bins(bin+1))/2;
end

% trim at amp points and end points
trim_logical = ccrot(:,1)<=(endPointrot(1)-ampPointrot(1)) & ccrot(:,1)>=0; 
 

% add to structure
single_exp.ktr = ktr;
single_exp.averageKTR = averageKTR;
single_exp.semKTR = semKTR;


single_exp.gem = gem;
single_exp.nucleiGEM = nucleiGEM;
single_exp.nucleiALL = nucleiALL;
single_exp.fractionGEM = fractionGEM;
single_exp.averageGEM = averageGEM;
single_exp.medianGEM = medianGEM;



single_exp.ccrot = ccrot;
single_exp.ampPointrot = ampPointrot;
single_exp.endPointrot = endPointrot;
single_exp.L_reg = endPointrot(1)-ampPointrot(1);


single_exp.trim_logical = trim_logical;

single_exp.binvalue = binvalue;


end

