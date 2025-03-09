function [length_mat,seeds,fiducial,labeled_seeds,fiducial_fullimage,imQCBFHere,imQCosxHere] = segment_initial_seeds_STvec(length_mat,fish,s,BF,osx,PB_bw,fin_bw,dpa_ls,options)
arguments
    length_mat
    fish
    s
    BF
    osx
    PB_bw
    fin_bw
    dpa_ls
    options.verbose = [];
end
%SEGMENT_INITIAL_SEEDS Summary of this function goes here
%   Detailed explanation goes here
%% outer bound
% thres_otsu = graythresh(PB);
% PB_bw = imbinarize(PB,thres_otsu);
% PB_bw = PB>0.5;
% figure;
% imshow(PB_bw,[])

% morphological


% PB_bw = bwmorph(PB_bw,'clean');
% PB_bw = bwmorph(PB_bw,'fill');
% tic;
mask_out = PB_bw;
mask_out(~fin_bw) = 0;

mask_out = bwareaopen(mask_out,100);
mask_out_closed = fin_bw;



% mask_hull = bwconvhull(mask_out);

% figure; imshow(mask_out);
% figure; imshow(mask_out_closed);
% toc
%% take upper region of the outer bound
% mask_out2 = false(size(mask_hull));
% top_line = false(size(mask_out_closed));
top_line = [];
for c = 1:size(mask_out_closed,2)
    idx_1 = find(mask_out_closed(:,c),1);
%     mask_out2(idx_1:idx_1+200,c) = true;
%     top_line(idx_1,c) = true;
    if ~isempty(idx_1)
        top_line = [top_line;[idx_1,c]];
    end
end
left_line = [];
right_line = [];
for c = 1:size(mask_out_closed,1)
    idx_1 = find(mask_out_closed(c,:),1,'first');
    idx_2 = find(mask_out_closed(c,:),1,'last');
%     mask_out2(idx_1:idx_1+200,c) = true;
    if ~isempty(idx_1)
        left_line= [left_line;[c,idx_1]];
    end
    if ~isempty(idx_2)
        right_line = [right_line;[c,idx_2]];
    end
end
% figure; imshow(top_line)

% old method - 200 pixels down from the top line
% mask_fiducial = mask_out.*mask_out2;
% im_fiducial = PB_bw.*mask_out2;
% im_fiducial = bwareaopen(im_fiducial,2000);
% 
% % figure; imshow(im_fiducial);
% CC = bwconncomp(im_fiducial);
%% Find the two spikes on the sides and use as markers for initial seeds
[~,maxidx] = max(top_line(:,2));
[~,minidx] = min(top_line(:,2));
top_line_row_left = top_line(minidx,1);
top_line_row_right = top_line(maxidx,1);

% find the longer edge
if top_line_row_left<top_line_row_right 
    row = left_line(:,1);
    col = size(mask_out_closed,2)-left_line(:,2);
    trunc = row>top_line_row_left;
    row = row(trunc);
    col = col(trunc);
else
    row = right_line(:,1);
    col = right_line(:,2);
    trunc = row>top_line_row_right;
    row = row(trunc);
    col = col(trunc);
end

% figure; imshow(mask_out);
% old method - use sum of pixels in row
    % [row,col] = find(top_line);
%     row = [1:size(mask_out_closed,1)]';
%     col = sum(mask_out_closed,2);
%     use = find(sum(top_line,2),1,'last')+10:size(mask_out_closed,1); % the lowest point of top_line
%     x_use = row(use);
%     y_use = col(use);
    x_use = row;
    y_use = col;
    fit_obj = fit(x_use,y_use,'smoothingspline', SmoothingParam = 1E-5);
    
    cut_idx = x_use(find(diff(fit_obj(x_use))>0,1));
    
%     figure;plot(x_use,(fit_obj(x_use)),'-')
%     figure;plot(x_use(1:end-1),diff(fit_obj(x_use)),'-')
if cut_idx-max(top_line_row_left,top_line_row_right)>100
    line_fiducial = cut_idx;
elseif cut_idx-max(top_line_row_left,top_line_row_right)>50
    line_fiducial = cut_idx+50;
else
    line_fiducial = cut_idx+150;
end

im_fiducial = mask_out;
im_fiducial([1:line_fiducial-200,line_fiducial:end],:) = 0;
im_fiducial = bwareaopen(im_fiducial,2000);

% figure; imshow(im_fiducial);
CC = bwconncomp(im_fiducial);

im_rays = mask_out;
im_rays([line_fiducial:end],:) = 0;
im_rays = bwareaopen(im_rays,2000);
labeled_rays = bwlabel(im_rays);
% figure; imshow(im_rays,[]);
%% iterate to segment_initial_seeds
% initialize
seeds = cell([1,numel(CC.PixelIdxList)]);
fiducial = cell([1,numel(CC.PixelIdxList)]);
QC_arcBW = uint8(zeros(size(PB_bw)));

labeled_seeds = uint8(zeros(size(PB_bw)));
fiducial_fullimage = cell(size(fiducial));

for raynum = 1:numel(CC.PixelIdxList)
    % raynum
    %% find one ray
    % raynum = 3;
    BW_oneray_seed = false(size(im_fiducial));
    BW_oneray_seed(CC.PixelIdxList{raynum}) = true;
    BW_oneray = labeled_rays==mode(labeled_rays(BW_oneray_seed),'all'); % groupcounts % find the largest piece

    % figure; imshow(BW_oneray);

    %% save seeds
    labeled_seeds(BW_oneray_seed) = uint8(raynum); % record for QC
        % figure; imshow(labeled_seeds,[]);
    %% crop canvas and measure length of one ray
    [~,~,end_bot,BW_oneray_small,~,row_excerpt,col_excerpt] = measure_one_ray(BW_oneray_seed,verbose = options.verbose,initial_seeds = 1);
    
    
    %% save length, seed and end_bot
    length_matHere = struct('fish',fish,'ray',raynum,'dpa',dpa_ls,'arcL',nan(size(dpa_ls)),'fiducial',nan([numel(dpa_ls),2]));
%     length_matHere.arcL(s) = arcLHere;
    length_matHere.fiducial(s,:) = end_bot+[row_excerpt(1) col_excerpt(1)]-1;
    
    if isempty(length_mat)||~any([length_mat.fish]==fish&[length_mat.ray]==raynum,"all")
        length_mat = [length_mat;length_matHere];
    else
        length_mat([length_mat.fish]==fish&[length_mat.ray]==raynum) = length_matHere;
    end
    % seeds
    BW_oneray_small_seed = imdilate(BW_oneray_small,strel('disk',100,8));
    BF_oneray_small = BF(row_excerpt,col_excerpt);
    BF_oneray_small(~BW_oneray_small_seed) = 0;
    seeds{raynum} = cat(3,BF_oneray_small,im2uint8(BW_oneray_small_seed),im2uint8(BW_oneray_small));
    % figure; imshow(BF_oneray_small);
    % fiducial
    fiducial{raynum} = end_bot;
    fiducial_fullimage{raynum} = end_bot+[row_excerpt(1),col_excerpt(1)]-1;

    
    % close all;
end
%% plot QC
% BF_rgb = repmat(BF,[1,1,3]);
% QC_rgb = repmat(QC_arcBW,[1,1,3]);
% cm = [0,0,0;lines(30)]; % colormap
% imQCHere = im2uint8(ind2rgb(QC_arcBW,cm));
% imQCHere(QC_rgb==0) = BF_rgb(QC_rgb==0);

% figure;imshow(imQCHere);
% imQCHere = imfuse(BF,QC_arcBW,'ColorChannels',[2 1 2]);

imQCBFHere = labeloverlay(BF,QC_arcBW,Colormap=lines(30),Transparency=0);
imQCosxHere = labeloverlay(osx,QC_arcBW,Colormap=lines(30),Transparency=0);

%% Debug Boundary
% BD_mask = false(size(BW_oneray_small_closed));
% BD_mask(sub2ind(size(BD_mask),BD(:,1),BD(:,2))) = true;
% figure;
% imshow(BD_mask)



%% Debug2 Bright field
% BF2 = BF;
% BF2(~im_fiducial) = 0;
% % figure; imshowpair(BF,im_fiducial);
% figure; imshow(BF2);

end

