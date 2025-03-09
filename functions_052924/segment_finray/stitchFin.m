function [imHere, tform] = stitchFin(imPos1, imPos2, varargin)
%STITCHFIN Summary of this function goes here
%   Detailed explanation goes here
if nargin == 3 && (isa(varargin{1}, 'simtform2d') || isa(varargin{1}, 'rigidtform2d'))
    tform = varargin{1};
elseif nargin == 2
    % calculate transformation
%     tic
%     imPos1_half = imPos1;
%     imPos1_half(1:floor(size(imPos1_half,1)/2),:) = 0;
%     tform = imregcorr(imPos1_half,imPos2);
    tform = matlabSIFTregist(imPos1,imPos2);
%     toc
else
    error("only 3 input arguments required and the third one needs to be a tform object");
end

% transform Pos1 and its mask
% R_imPos1 =
% affineOutputView(size(imPos1),tform,"BoundsStyle","FollowOutput"); % more cumbersome way to get the R_imPos1
% imPos1_T = imwarp(imPos1,tform,"OutputView",R_imPos1);

[imPos1_T, R_imPos1] = imwarp(imPos1,tform);
R_imPos2 = imref2d(size(imPos2));

maskPos1_T = imwarp(true(size(imPos1)-1),tform,"OutputView",R_imPos1);

% overlay transformed images to new positions in the canvas
[imPos1_reg,imPos2_reg] = calculateOverlayImages(imPos1_T,R_imPos1, imPos2,R_imPos2);
[maskPos1_reg] = calculateOverlayImages(maskPos1_T,R_imPos1, imPos2,R_imPos2);

% merge
imHere = imPos2_reg;
maskPos1_erode = imerode(maskPos1_reg, strel("diamond",1));
imHere(maskPos1_erode) = imPos1_reg(maskPos1_erode);

%%
% figure;
% imshow(imHere)
% imshowpair(imPos1_reg,imPos2_reg, 'falsecolor');
% imshowpair(imPos1_T,R_imPos1, imPos2,R_imPos2, "falsecolor");
end

