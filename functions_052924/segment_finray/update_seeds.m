function [seeds_new,fiducial_new,labeled_seeds,fiducial_fullimage] = ...
    update_seeds(seeds,fiducial,BF,PB_bw,options)
arguments
    seeds;
    fiducial;
    BF;
    PB_bw;
    options.BW_oneray_cell = [];
    options.verbose = 0;
end
%UPDATE_SEEDS Summary of this function goes here
%   Detailed explanation goes here
%% outer bound
mask_out = PB_bw;

mask_out = bwareaopen(mask_out,500);
% PB_bw = imclose(PB_bw,strel('disk',10));

mask_out = imdilate(mask_out,strel('disk',100,8));
mask_out = imfill(mask_out,'holes');
mask_out = imerode(mask_out,strel('disk',100,8));

mask_out = imdilate(mask_out,strel('disk',20,8));

% mask_out = bwconvhull(mask_out);

% figure; imshow(mask_out);

% crop BF
BF_fin = BF;
BF_fin(~mask_out) =  0;
% figure; imshow(BF_fin);

%% loop over the rays to get seeds_new
seeds_new = cell(size(seeds));
fiducial_new = cell(size(fiducial));
labeled_seeds = uint8(zeros(size(PB_bw)));
% labeled_seeds_watershed = uint8(zeros(size(PB_bw))); % alternative fullray watershed
fiducial_fullimage = cell(size(fiducial));

for raynum = 1:numel(seeds)
    %% find overlap
    % close all;
    % raynum = 9;
    
%     raynum
    seedHere = seeds{raynum}; 
    fiducialHere = fiducial{raynum};
    
    % take a look
%     figure;
%     imshow(seedHere); hold on;
%     plot(fiducialHere(2),fiducialHere(1),'r*')
% %     plot(end_bot(1),end_bot(2),'r*')
%     hold off;
    
    % start registration
    imPos1 = seedHere(:,:,1);
    maskPos1 = logical(seedHere(:,:,2));
    seedPos1 = logical(seedHere(:,:,3));
    imPos2 = BF_fin;
    
    %% method1 cross corr
%     corrhere = normxcorr2(imPos1,imPos2);
%     [ypeak, xpeak] = find(corrhere==max(corrhere(:)),1,'first');
%     figure;
%     surf((corrhere))
%     shading flat
    
    % %% method3 SIFT
%     tic;
    [tform,match_numSIFT1] = matlabSIFTregist(imPos1,imPos2,...
        selectStrongest=[],...
        MatchThreshold=100, ...
        MaxDistance=5,Confidence=99.9,MaxNumTrials=10000,...
        verbose = options.verbose);
    if match_numSIFT1<10
%         disp('1st SIFT returns too few matching features, increase features for SIFT')
%         [tform,match_numSIFT2] = matlabSIFTregist(imPos1,imPos2,...
%             selectStrongest=[],sigma = 1.6,EdgeThreshold=50,ContrastThreshold=0.003,NumLayersInOctave=4,...
%             MatchThreshold=100, ...
%             MaxDistance=5,Confidence=99.9,MaxNumTrials=10000,...
%             verbose = options.verbose);
%         if match_numSIFT2<10 % use kaze feature if SIFT doesnt work
            disp('SIFT returns too few matching features, use KAZE')
            [tform] = matlabSIFTregist(imPos1,imPos2,...
                selectStrongest=[], ...
                kaze = 1, ...
                MatchThreshold=100, ...
                MaxDistance=5,Confidence=99.9,MaxNumTrials=10000,...
                verbose = options.verbose);
%         end
    end
%     toc
    % bu


    % transform Pos1 and its mask to new positions in Pos2 canvas
    
    R_imPos2 = imref2d(size(imPos2));
    [imPos1_reg, R_imPos1] = imwarp(imPos1,tform,"OutputView",R_imPos2);

    maskPos1_reg = imwarp(maskPos1>0,tform,"OutputView",R_imPos2);
    seedPos1_reg = imwarp(seedPos1,tform,"OutputView",R_imPos2);
    
    % overlay transformed images to new positions in the canvas !!!!!!!!! calculateOverlayImages deform images
%     [imPos1_reg,imPos2_reg] = calculateOverlayImages(imPos1_T,R_imPos1, imPos2,R_imPos2);
%     [maskPos1_reg] = calculateOverlayImages(maskPos1_T,R_imPos1, imPos2,R_imPos2);
%     [seedPos1_reg] = calculateOverlayImages(seedPos1_T,R_imPos1, imPos2,R_imPos2);
    
    %% alternative fullray watershed
    % transform oneray (full length of a ray for better seed for watershed)
%     BW_oneray_here = BW_oneray_cell{1,raynum};
%     BW_oneray_fiducial_here = BW_oneray_cell{2,raynum};
%     offset_here = fiducial{raynum}-BW_oneray_fiducial_here;
%     xWorldLimits = 0.5+[0,size(BW_oneray_here,2)]+offset_here(2);
%     yWorldLimits = 0.5+[0,size(BW_oneray_here,1)]+offset_here(1);
%     R_oneray = imref2d(size(BW_oneray_here),xWorldLimits,yWorldLimits);
%     BW_oneray_reg = imwarp(BW_oneray_here,R_oneray,tform,"OutputView",R_imPos2);
% 
%     % take only the largest object
%     BW_oneray_watershed = BW_oneray_reg&PB_bw;
%     CC = bwconncomp(BW_oneray_watershed);
%     numPixels = cellfun(@numel,CC.PixelIdxList);
%     [~,maxidx] = max(numPixels);
%     BW_oneray_watershed(:) = 0;
%     BW_oneray_watershed(CC.PixelIdxList{maxidx}) = 1;
    %% calculate fiducial
    v_fullimage = tform.A*[fiducialHere(2) fiducialHere(1) 1]';
    v_fullimage = double(round(v_fullimage([2 1])));
    
    %% plot and view the registration
    if options.verbose
        figure;
        imshowpair(imPos1_reg,imPos2)
        hold on;
        plot(v_fullimage(2),v_fullimage(1),'y*')
        hold off;
        figure;
        imshow(maskPos1_reg)
        figure;
        imshowpair(seedPos1_reg,PB_bw)
        
    
        figure;
        imshowpair(BF,PB_bw)
    end
    %% add to new seeds and fiducial
    [row,col] = find(maskPos1_reg);
    row_excerpt = min(row):max(row);
    col_excerpt = min(col):max(col);
    
    % seed for registration
    BF_seed_new = imPos2;
    BF_seed_new(~maskPos1_reg)=0;
    BF_seed_new_small = BF_seed_new(row_excerpt,col_excerpt);
    maskPos1_reg_small = maskPos1_reg(row_excerpt,col_excerpt);
    % new binary seed for watershed
    seedPos1_reg_overlap = imclose(seedPos1_reg,strel('disk',20))&PB_bw; % take overelap between transformed seed and new mask
    seedPos1_reg_overlap = imdilate(seedPos1_reg_overlap,strel('line',20,0))&PB_bw; % dilate to prevent shrinking

    if any(seedPos1_reg_overlap,"all") % should have overlap with mask
        % the part between tic and toc is to make sure that only minima larger than 30% area is selected during watershed for each seed
%     tic;
        CC = bwconncomp(seedPos1_reg_overlap);
        numPixels = cellfun(@numel,CC.PixelIdxList);
%         retain_idx = numPixels>=sum(numPixels,"all").*0.3;
        [~,maxidx] = max(numPixels);
        seedPos1_reg_overlap(:) = 0;
        seedPos1_reg_overlap(cat(1,CC.PixelIdxList{maxidx})) = 1;
    %     toc
        seedPos1_reg_small = seedPos1_reg_overlap(row_excerpt,col_excerpt);
        seedHere_new = cat(3,BF_seed_new_small,im2uint8(maskPos1_reg_small),im2uint8(seedPos1_reg_small));
        
        fiducialHere_new = v_fullimage(:)'-[min(row),min(col)]+1;
    else
        seedHere_new = seedHere;
        fiducialHere_new = fiducialHere;
    end

%     toc
    %% take a look at the new seeds and fiducial
%     figure;imshowpair(seedPos1_reg_overlap,PB_bw)
    % figure;
    % imshow(seedHere_new); hold on;
    % plot(fiducialHere_new(2),fiducialHere(1),'r*')
    % hold off;
    % 
    %% add to cell list
    seeds_new{raynum} = seedHere_new;
    fiducial_new{raynum} = fiducialHere_new;
    
    %% generate labeled seed image for watershed
    previous_seeds = labeled_seeds>0;
    labeled_seeds(seedPos1_reg_overlap) = uint8(raynum);
    labeled_seeds(imdilate(previous_seeds,strel('square',3))&seedPos1_reg_overlap) = 0;
%     labeled_seeds_watershed(BW_oneray_watershed) = labeled_seeds_watershed(BW_oneray_watershed)+raynum; % alternative fullray watershed
    fiducial_fullimage{raynum} = v_fullimage(:)';

    %% debug

end
%     labeled_seeds(labeled_seeds>numel(seeds)) = 0; % make overlapping region 0
%     labeled_seeds_watershed(labeled_seeds_watershed>numel(seeds)) = 0; % alternative fullray watershed
%% remove outliers
coord = cell2mat(fiducial_fullimage');

too_far = hypot(coord(1:end-1,1)-coord(2:end,1),coord(1:end-1,2)-coord(2:end,2))>300;

idx_remove = find([0;too_far]&[too_far;0]);
if too_far(1)==1&&too_far(2)==0
    idx_remove = [1;idx_remove];
end
if too_far(end)==1&&too_far(end-1)==0
    idx_remove = [idx_remove;size(coord,1)];
end

[seeds_new{idx_remove}]= seeds{idx_remove};
[fiducial_new{idx_remove}] = fiducial{idx_remove};
if ~isempty(idx_remove)
    for id = 1:numel(idx_remove)
        fiducial_fullimage{idx_remove(id)} = [];
    end
end

labeled_seeds(ismember(labeled_seeds,idx_remove)) = 0;
% figure;imshow(labeled_seeds,[])
% figure;
% plot(coord(:,1),coord(:,2),'k*','MarkerSize',5);
%%
end

