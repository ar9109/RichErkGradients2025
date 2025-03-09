function [PB_bw] = ilastikPB2BW(PB,options)
arguments
    PB
    options.thres = 0.5;
end
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% old using ilastik mask
PB_bw = PB>options.thres;
% figure; imshow(PB_bw,[])
% morphological


% PB_bw = imopen(PB_bw,strel('disk',1));
% PB_bw = imfill(PB_bw,'holes');

PB_bw = bwmorph(PB_bw,'clean');
PB_bw = bwmorph(PB_bw,'fill');
PB_bw = bwareaopen(PB_bw,5);

% PB_bw = imfill(PB_bw,'holes');

% PB_bw = imdilate(PB_bw,strel('disk',1));
% PB_bw = imfill(PB_bw,'holes');
% PB_bw = imerode(PB_bw,strel('disk',1));

PB_bw = bwareaopen(PB_bw,100);


% figure; imshow(PB_bw);
end

