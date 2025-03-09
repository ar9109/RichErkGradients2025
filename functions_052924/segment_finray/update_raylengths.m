function [length_mat,imQCBFHere,imQCosxHere,L_display] = update_raylengths(length_mat,labeled_seeds,fiducial_fullimage,BF,osx,PB_bw,fish,s,options)
arguments
    length_mat;
    labeled_seeds;
    fiducial_fullimage;
    BF;
    osx;
    PB_bw;
    fish;
    s;
    options.verbose = [];
end
%UPDATE_RAYLENGTHS Summary of this function goes here
%   Detailed explanation goes here
%% truncate mask below seeds
% [rowSeeds,colSeeds] = find(labeled_seeds);
% cutoff = max(rowSeeds(:))+100;

closed_seeds = imclose(labeled_seeds>0,strel('disk',200));
[~,maxidxs] = max(cumsum(closed_seeds),[],1);
maxidxs(maxidxs==1) = inf;
[~,YY] = meshgrid(1:size(closed_seeds,2),1:size(closed_seeds,1));
truncation_mask = YY<=maxidxs;

% cutoff = max(cellfun(@(c)c(1),fiducial_fullimage),[],'all')+100;
PB_bw_trunc = PB_bw;
PB_bw_trunc(~truncation_mask) = 0;
PB_bw_trunc_dilated = imdilate(PB_bw_trunc,strel('disk',3));
% figure; imshow(PB_bw_trunc_dilated)
% figure; imshow(closed_seeds)
%% find upper bound
mask_out = PB_bw;

mask_out = bwareaopen(mask_out,500);
% PB_bw = imclose(PB_bw,strel('disk',10));

mask_out = imdilate(mask_out,strel('disk',50,8));
mask_out = imfill(mask_out,'holes');
mask_out = imerode(mask_out,strel('disk',50,8));

% make sure only the largest object preserved
CC = bwconncomp(mask_out);
[~,maxidx] = max(cellfun(@numel,CC.PixelIdxList));
mask_out = false(size(mask_out));
mask_out(CC.PixelIdxList{maxidx}) = 1;

% mask_out = imdilate(mask_out,strel('disk',20,8));
% figure; imshow(mask_out)

% upper line of mask_out
% find the two higest points
BD_y = zeros(1,(size(mask_out,2)));
for c = 1:size(mask_out,2)
    idx_1 = find(mask_out(:,c),1);
    if ~isempty(idx_1)
        BD_y(c) = size(mask_out,1)-idx_1+1;
    end
end
x_anchor_left = find(BD_y>0,1,"first");
x_anchor_right = find(BD_y>0,1,"last");

BD_cell = bwboundaries(imdilate(mask_out,strel('disk',1))); 
[~,maxidx] = max(cellfun(@numel,BD_cell));
BD_coor = BD_cell{maxidx};

dist_anchor_left = hypot(BD_coor(:,1)-0,BD_coor(:,2)-x_anchor_left);
dist_anchor_left(abs(BD_coor(:,2)-x_anchor_left)>500) = max(dist_anchor_left); % make sure the scan doesnt go beyond 200px (horizontal) of the anchor
dist_anchor_right = hypot(BD_coor(:,1)-0,BD_coor(:,2)-x_anchor_right);
dist_anchor_right(abs(BD_coor(:,2)-x_anchor_right)>500) = max(dist_anchor_right); % make sure the scan doesnt go beyond 200px (horizontal) of the anchor

[~,minidx_left] = min(dist_anchor_left);
[~,minidx_right] = min(dist_anchor_right);

if minidx_left<=minidx_right
    upper_bound_coor = BD_coor(minidx_left:minidx_right,:);
else
    upper_bound_coor = BD_coor([minidx_right:end,1:minidx_left],:);
end
BD_up = false(size(mask_out));
BD_up(sub2ind(size(mask_out),upper_bound_coor(:,1),upper_bound_coor(:,2))) = 1;

% BD = imdilate(mask_out,strel('disk',1))-mask_out;
% figure;imshow(BD);
% figure;imshow(BD_up);
% figure;imshowpair(BD_up,BD);



%% perform watershed
% define minima from seeds
watershed_minima = labeled_seeds>0;
% watershed_minima(~mask_out) = max(watershed_minima(:))+1;
% watershed_minima = imerode(watershed_minima>0,strel('disk',1));
watershed_minima(~PB_bw_trunc_dilated) = 1; % define minima for background

% figure; imshowpair(PB_bw_trunc,watershed_minima);

% construct propagation matrix with distance transform and upwards queue

% !!method1
% use absolute numbers of vertical increment to control watershed
% [~,YY] = meshgrid(1:size(PB_bw_trunc,2),1:size(PB_bw_trunc,1)); 
% D = flipud(YY); 
% D2 = -bwdist(BD_up);
% % distance transform from the upper bound
% D2_mask = D2>-200;
% D2(~D2_mask) = 0;
% % merge two matrices
% D(D2_mask) = D2(D2_mask)+size(D,1)+abs(min(D2(:)))+100;

% !!method2
D1 = bwdist(closed_seeds);
D2 = bwdist(BD_up);
D = D1./(D1+D2);

% !! method3 find potential from the structure tensor % seem to perform worse than method2
% tic;
% end_bot = fiducial_fullimage{15};
% mask_in = imclose(PB_bw_trunc_dilated,strel('disk',50));
% [potential] = get_watershed_potential(osx,mask_in,end_bot);
% toc
% D = potential-min(potential(:))+1;

% !! look at isosurface
% a = 0.16;
% lowbound = min(D(:))+a.*(max(D(:))-min(D(:)));
% highbound = min(D(:))+(a+0.05).*(max(D(:))-min(D(:)));
% figure;imshow(D>lowbound&D<highbound);
% figure;imshowpair(D>lowbound&D<highbound,closed_seeds);

% !!continue
% queue the background after foreground
D(~PB_bw_trunc) = max(D,[],'all').*1.1;
% queue the boundary last
D(logical(imdilate(PB_bw_trunc,strel('disk',1,8))-PB_bw_trunc)) = max(D,[],'all').*1.1; %!!! change to outside boundary

% figure; imshow(D,[])
% figure; imshowpair(D,mask_out2);

% impose minima for seeds
J = imimposemin(D,watershed_minima);
% figure; imshow(J,[])

% watershed
L = watershed(J);
% L = watershed_mod(D,watershed_minima); % the modified watershed function
L(~PB_bw) = 0;
% figure; imshow(L,[])

% look at watershed
L_display = labeloverlay(im2uint8(PB_bw_trunc),L,Colormap=colorcube(double(max(L(:))+10)),Transparency=0);
% figure; imshow(L_display)


%% assign watershed objects to rays and measure lengths
QC_arcBW = uint8(zeros(size(labeled_seeds)));
% if ~istable(length_mat) % better to use structure in this code
%     length_mat = struct2table(length_mat);
% end

% alternative fullray watershed
% BW_oneray_cell = cell(2,numel(fiducial_fullimage)); 


for raynum = unique(labeled_seeds(labeled_seeds>0))'
%     close all;
    %% find corresponding ray in watershed L and truncate with fiducial
%     raynum
    BW_oneray_seed = labeled_seeds==raynum;
    % figure; imshow(labeled_seeds)
    
    BW_oneray = L==mode(L(BW_oneray_seed),'all'); % find the largest piece
%     [B,BG] = groupcounts(L(BW_oneray_seed&L>0&~isinf(L))); % groupcounts
%     BW_oneray = ismember(L,BG(B>0.3*sum(B,"all"))); % find all pieces that match with the seeds, but larger than 30% total area
    BW_oneray(fiducial_fullimage{raynum}(1)+1:end,:) = 0;
%     figure; imshow(BW_oneray)
    if ~any(BW_oneray,'all')
        continue
    end
    %% measure arcLength
    [arcBW,arcLHere,end_bot,BW_oneray_small,BW_oneray_small_closed,row_excerpt,col_excerpt]...
        = measure_one_ray(BW_oneray,verbose = options.verbose);

    % alternative fullray watershed
%     BW_oneray_cell{1,raynum} = BW_oneray_small;
%     BW_oneray_cell{2,raynum} = fiducial_fullimage{raynum} - [row_excerpt(1),col_excerpt(1)] +1;
    %% get QC overlaid image
    if ~isempty(options.verbose)&&options.verbose
        figure;
        imshowpair(BW_oneray_small_closed,arcBW)
    end
    [arc_row,arc_col] = find(arcBW);
    original_coor_row = arc_row+row_excerpt(1)-1;
    original_coor_col = arc_col+col_excerpt(1)-1;
    % paint
    QC_arcBW(sub2ind(size(QC_arcBW),original_coor_row,original_coor_col)) = raynum;
    %% save length
    thisray = [length_mat.fish]==fish&[length_mat.ray]==raynum; % better to use structure in this code
    length_mat(thisray).arcL(s) = arcLHere;
    length_mat(thisray).fiducial(s,:) = fiducial_fullimage{raynum};
end

% length_mat = table2struct(length_mat); % better to use structure in this code
% figure; imshow(QC_arcBW,[]);

imQCBFHere = labeloverlay(BF,QC_arcBW,Colormap=lines(30),Transparency=0);
imQCosxHere = labeloverlay(osx,QC_arcBW,Colormap=lines(30),Transparency=0);
end

