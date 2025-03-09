function [stack,roi,roiTight] = selectROIstack(stack,roi,opts)
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

if(isempty(roi))
    roi=ones(size(stack),'logical');
end

disp('Please select the scale of interest in the three projections (Press Enter to move to the next view)');
 

% whether to use BF, edited 110922
if(isfield(opts,'BF')) 
    if ~isempty(opts.BF)
        
        stack_BF = loadtiff(opts.BFpath);
        h=figure;
        imshow(imadjust(stack_BF));
        [mask,xi2]=roipoly;
        close(h);
        if ~(isempty(mask) || numel(xi2)<=4)
            for z=1:size(stack,3)
                img=roi(:,:,z);
                img(~mask)=0;
                roi(:,:,z)=img;
            end
            stack(~roi)=0;
        end

    end
end


% Y-X
proj=max(stack,[],3);
h=figure;
imshow(imadjust(proj,adj,[0 1]));
[mask,xi2]=roipoly;
close(h);  % edited 080422
if ~(isempty(mask) || numel(xi2)<=4)
    for z=1:size(stack,3)
        img=roi(:,:,z);
        img(~mask)=0;
        roi(:,:,z)=img;
    end
    stack(~roi)=0;
end


% Y-Z
stack=permute(stack,[1 3 2]);
roi  =permute(roi  ,[1 3 2]);

% close(h);
proj=max(stack,[],3);
h=figure;
imshow(imadjust(proj,adj,[0 1]));
[mask,xi2]=roipoly;
close(h); % edited 080422
if ~(isempty(mask) || numel(xi2)<=4)
    for z=1:size(stack,3)
        img=roi(:,:,z);
        img(~mask)=0;
        roi(:,:,z)=img;
    end
    stack(~roi)=0;
end
% stack=permute(stack,[1 3 2]);
% roi  =permute(roi  ,[1 3 2]);


% Z-X
stack=permute(stack,[2 3 1]);
roi  =permute(roi  ,[2 3 1]);
proj=max(stack,[],3);
% close(h);
h=figure;
imshow(imadjust(proj,adj,[0 1]));
[mask,xi2]=roipoly;
close(h);
if ~(isempty(mask) || numel(xi2)<=4)
    for z=1:size(stack,3)
        img=roi(:,:,z);
        img(~mask)=0;
        roi(:,:,z)=img;
    end
    stack(~roi)=0;
end
stack=permute(stack,[3 2 1]);
roi  =permute(roi  ,[3 2 1]);

% edited 071822
if (isfield(opts,'roiTight')) && opts.roiTight
        roiTight = roi;
        proj=max(stack,[],3);
        h=figure;
        imshow(imadjust(proj,adj,[0 1]));
        mask=roipoly;
        if(isempty(mask))
            mask=ones(size(proj));
        end
        
        for z=1:size(stack,3)
            img=roiTight(:,:,z);
            img(~mask)=0;
            roiTight(:,:,z)=img;
        end
else
    roiTight = [];
end

end


