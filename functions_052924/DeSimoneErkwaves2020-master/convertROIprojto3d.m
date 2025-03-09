function [roi] = convertROIprojto3d(myScale,paths, suffix)
arguments
    myScale
    paths
    suffix = 'ROIprojYX';
end
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
path_YX = [paths.fastFolder filesep myScale.basename '_' suffix];
roi_YX=logical(loadtiff([path_YX]));
path_YZ = [paths.fastFolder filesep myScale.basename '_ROIprojYZ.tif'];
roi_YZ=logical(loadtiff([path_YZ]));
path_ZX = [paths.fastFolder filesep myScale.basename '_ROIprojZX.tif'];
roi_ZX=logical(loadtiff([path_ZX]));

roi=ones([size(roi_YX),size(roi_YZ,2)],'logical');

% Y-X
mask = roi_YX;
for z=1:size(roi,3)
    img=roi(:,:,z);
    img(~mask)=0;
    roi(:,:,z)=img;
end



% Y-Z
roi  =permute(roi  ,[1 3 2]);

mask = roi_YZ;
for z=1:size(roi,3)
    img=roi(:,:,z);
    img(~mask)=0;
    roi(:,:,z)=img;
end


% Z-X
roi  =permute(roi  ,[2 3 1]);

mask = roi_ZX;
for z=1:size(roi,3)
    img=roi(:,:,z);
    img(~mask)=0;
    roi(:,:,z)=img;
end

roi  =permute(roi  ,[3 2 1]);

% sancheck
% stack_YX = max(roi,[],3);
% stack_YZ = squeeze(max(roi,[],2));
% stack_ZX = squeeze(max(roi,[],1))';
% sum(stack_YX~=roi_YX,'all')
% sum(stack_YZ~=roi_YZ,'all')
% sum(stack_ZX~=roi_ZX,'all')

end

