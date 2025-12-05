function [stack,roi,roiTight] = selectROIstack_fast(stack,roi,opts)
%  Copyright (C) 2020  Alessandro De Simone
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%
%
adj=opts.adj;

% if(isempty(roi))
%     roi=ones(size(stack),'logical');
% end
roi = [];


stack_YX = max(stack,[],3);
stack_YZ = squeeze(max(stack,[],2));
stack_ZX = squeeze(max(stack,[],1))';

% stack_YX = stack.YX;
% stack_YZ = stack.YZ;
% stack_ZX = stack.ZX;

disp('Please select the scale of interest in the three projections (Press Enter to move to the next view)');
 

% whether to use BF, edited 110922
if(isfield(opts,'BF')) && ~isempty(opts.BF)
        
    stack_BF = loadtiff(opts.BFpath);
    h=figure;
    h.Name = [opts.basename ' BF']; % edited 032624
    imshow(imadjust(stack_BF));
    h.WindowState = 'maximized';
    [mask_BF,xi2]=roipoly;
    close(h);
    if ~(isempty(mask_BF) || numel(xi2)<=4)
        stack_YX(~mask_BF)=0;
        stack_YZ(~any(mask_BF,2),:)=0;
        stack_ZX(:,~any(mask_BF,1))=0;
    else
        mask_BF = ones(size(stack_YX));
    end
    % save roi.BF
    roi.BF = mask_BF;

else
    mask_BF = ones(size(stack_YX));
end


% Y-X
h=figure;
h.Name = [opts.basename ' XY']; % edited 032624
imshow(adapthisteq(stack_YX,NumTiles=round([size(stack_YX)./100]))); % edited 032724
% imshow(imadjust(stack_YX,adj,[0 1]));
h.WindowState = 'maximized';
[mask,xi2]=roipoly;
close(h);  % edited 080422
if ~(isempty(mask) || numel(xi2)<=4)
    stack_YX(~mask)=0;
    stack_YZ(~any(mask,2),:)=0;
    stack_ZX(:,~any(mask,1))=0;
else
    mask = ones(size(stack_YX));
end
roi_YX = mask & mask_BF;
clear mask_BF
clear mask;

% Y-Z

% close(h);
h=figure;
h.Name = [opts.basename ' YZ']; % edited 032624
imshow(imadjust(stack_YZ,adj,[0 1]));
h.WindowState = 'maximized';
[mask,xi2]=roipoly;
close(h); % edited 080422
if ~(isempty(mask) || numel(xi2)<=4)
    stack_YZ(~mask)=0;
    stack_YX(~any(mask,2),:)=0;
    stack_ZX(~any(mask,1),:)=0;
else
    mask = ones(size(stack_YZ));
end
roi_YZ = mask;
clear mask;


% Z-X
h=figure;
h.Name = [opts.basename ' ZX']; % edited 032624
imshow(imadjust(stack_ZX,adj,[0 1]));
h.WindowState = 'maximized';
[mask,xi2]=roipoly;
close(h);
if ~(isempty(mask) || numel(xi2)<=4)
%     stack_ZX(~mask)=0;
%     stack_YZ(:,~any(mask,2))=0;
%     stack_YX(:,~any(mask,1))=0;
else
    mask = ones(size(stack_ZX));
end
roi_ZX = mask;
clear mask;

% edited 071822
roiTight = [];
if isfield(opts,'roiTight') && ~isempty(opts.roiTight)
    for c = opts.roiTight
        stack_tight = loadtiff(opts.roiTightpath{c});
        h=figure;
        h.Name = [opts.basename ' Tight!!!']; % edited 032624
        imshow(imadjust(stack_tight));
        h.WindowState = 'maximized';
        [mask_tight,xi2]=roipoly;
        close(h);
        if (isempty(mask_tight) || numel(xi2)<=4)
        mask_tight = ones(size(stack_YX));
        end
        roiTight_here.channel = c;
        roiTight_here.mask = mask_tight;
        roiTight = [roiTight roiTight_here];
    end
end

roi.YX = roi_YX;
roi.YZ = roi_YZ;
roi.ZX = roi_ZX;
end


