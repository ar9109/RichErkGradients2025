function [fin_bw] = finmask(PB_bw,PB_BF)
arguments
    PB_bw
    PB_BF = [];
end
%FINMASK Summary of this function goes here
BW2 = PB_bw;

%% use bright field ilastik(doesnt work)
if ~isempty(PB_BF)
fin_bw1 = PB_BF>0.5;
% largest
CC = bwconncomp(fin_bw1);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,maxidx] = max(numPixels);
fin_bw1 = false(size(fin_bw1));
fin_bw1(CC.PixelIdxList{maxidx}) = 1;
fin_bw1 = imfill(fin_bw1,'holes');

% figure;imshow(fin_bw1)
BW2(~fin_bw1) = 0;

end
%%
% figure;imshow(BW2)











%% thresold
% % figure; imshow(IM,[])
% % PB2 = -imfilter(im2double(PB),fspecial('log',100,50));
% IM2 = imgaussfilt(IM_clean,1)-imgaussfilt(IM_clean,20);
% % PB2 = PB-imgaussfilt(PB,10);
% % IM2 = imgaussfilt(IM,1);
% % figure; imshow(IM2,[]);
% 
% thres_otsu = graythresh(IM2);
% % thres_adp = adaptthresh(PB2);
% BW = imbinarize(IM2,thres_otsu);
% % PB_bw = PB>100;
% figure; imshow(BW,[])

%% use edge detection
% IM2 = imgaussfilt(IM_clean,1)-imgaussfilt(IM_clean,100);
% figure; imshow(IM2,[]);
% 
% hx = [-1 -1 -1;2 2 2;-1 -1 -1];
% hxy = [2 -1 -1;-1 2 -1;-1 -1 2];
% hy = [-1 2 -1;-1 2 -1;-1 2 -1];
% hyx = [-1 -1 2;-1 2 -1;2 -1 -1];
% IM3 = sqrt(imfilter(im2double(IM2),hx).^2 ...
%     +imfilter(im2double(IM2),hxy).^2+ ...
%     imfilter(im2double(IM2),hy).^2+ ...
%     imfilter(im2double(IM2),hyx).^2);
% figure; imshow(IM3,[])
% 
% thres_otsu = graythresh(IM3);
% % thres_adp = adaptthresh(IM3);
% BW = imbinarize(IM3,thres_otsu);
% % PB_bw = PB3>0.03;
% figure; imshow(BW,[])
%% morphological
% PB_bw = imopen(PB_bw,strel("cube",3));
% BW2 = BW;
BW2 = bwmorph(BW2,'clean');
BW2 = bwareaopen(BW2,1000);
BW2 = bwmorph(BW2,'fill');
BW2 = imfill(BW2,'holes');


% take the largest region
BW2_closed = imclose(BW2,strel('disk',50));
% figure; imshow(BW2_closed,[])
CC = bwconncomp(BW2_closed);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,maxidx] = max(numPixels);
BW2_closed = false(size(BW2_closed));
BW2_closed(CC.PixelIdxList{maxidx}) = 1;
BW2(~BW2_closed) = 0;

% close
BW2 = imclose(BW2,strel('disk',100,8));
BW2 = imfill(BW2,'holes');

% figure; imshow(BW2,[])

%% merge
if ~isempty(PB_BF)
    fin_bw = fin_bw1&BW2;
else
    fin_bw = BW2;
end

end

