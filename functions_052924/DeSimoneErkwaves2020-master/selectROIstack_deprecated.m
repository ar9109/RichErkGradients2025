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

% extract bounding box directly from polygons edited 071822
rects = nan(4,4); 
rects(4,:) = 0;


adj=opts.adj;

if(isempty(roi))
    roi=ones(size(stack),'logical');
end

disp('Please select the scale of interest in the three projections (Press Enter to move to the next view)');
 
% X-Y
proj=max(stack,[],3);
h=figure;
imshow(imadjust(proj,adj,[0 1]));
% edited 071822
[mask,xi2,yi2] =roipoly;
minx = max(min(xi2)-30,1);
miny = max(min(yi2)-30,1);
maxx = min(max(xi2)+30,size(proj,2));
maxy = min(max(yi2)+30,size(proj,1));
rects(1,:) = [minx,miny,maxx-minx,maxy-miny];

% edited 071822
if (isfield(opts,'roiTight'))
    if opts.roiTight
        [mask,xi2,yi2] =roipoly;
        roiTight = [xi2,yi2];
    end
else
    roiTight = [];
end

% old
if(isempty(mask))
    mask=ones(size(proj));
end

for z=1:size(stack,3)
    img=roi(:,:,z);
    img(~mask)=0;
    roi(:,:,z)=img;
end

stack(~roi)=0;

% X-Z
stack=permute(stack,[1 3 2]);
roi  =permute(roi  ,[1 3 2]);

close(h);
proj=max(stack,[],3);


h=figure;
imshow(imadjust(proj,adj,[0 1]));

% edited 071822
[mask,xi2,yi2] =roipoly;
minx = max(min(xi2)-5,1);
miny = max(min(yi2)-30,1);
maxx = min(max(xi2)+5,size(proj,2));
maxy = min(max(yi2)+30,size(proj,1));
rects(3,:) = [minx,miny,maxx-minx,maxy-miny];

if(isempty(mask))
    mask=ones(size(proj));
end


for y=1:size(stack,3)
    img=roi(:,:,y);
    img(~mask)=0;
    roi(:,:,y)=img;
end

stack(~roi)=0;

stack=permute(stack,[1 3 2]);
roi  =permute(roi  ,[1 3 2]);

% Y Z
stack=permute(stack,[2 3 1]);
roi  =permute(roi  ,[2 3 1]);

proj=max(stack,[],3);
close(h);
h=figure;
imshow(imadjust(proj,adj,[0 1]));

% edited 071822
[mask,xi2,yi2] =roipoly;
minx = max(min(xi2)-5,1);
miny = max(min(yi2)-30,1);
maxx = min(max(xi2)+5,size(proj,2));
maxy = min(max(yi2)+30,size(proj,1));
rects(2,:) = [minx,miny,maxx-minx,maxy-miny];

close(h);
if(isempty(mask))
    mask=ones(size(proj));
end

for x=1:size(stack,3)
    img=roi(:,:,x);
    img(~mask)=0;
    roi(:,:,x)=img;
end

stack(~roi)=0;

stack=permute(stack,[3 1 2]);
roi  =permute(roi  ,[3 1 2]);

end


