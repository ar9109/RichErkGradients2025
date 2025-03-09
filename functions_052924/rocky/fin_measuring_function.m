%%fin_measuring_function

pixelSize = [3.632 3.632];

if exist('path_fin','var')
    [file,path_fin] = uigetfile({'*.tif'},'Select Fin Image',path_fin);
    fin = double(importdata(fullfile(path_fin,file)));
else
    [file,path_fin] = uigetfile({'*.tif'},'Select Fin Image');
    fin = double(importdata(fullfile(path_fin,file)));
end

num_rays = input("How many rays in the image? ");
masks = zeros(size(fin));
for i = 1:num_rays
    [file,path] = uigetfile({'*.tif'},strcat('Select Mask for Ray ',num2str(i)),path_fin);
    if isequal(file,0)
        disp('User selected Cancel');
    else
        disp(['User selected ', fullfile(path,file)]);
    end
    mask = fullfile(path,file);
    %tempray = double(importdata(strcat('C:\Users\Rocky\Desktop\mask_',num2str(i),'.tif')));
    tempray = double(importdata(mask));
    tempray = tempray.*fin;
    tempray = tempray./max(tempray(:));
    masks(:,:,i) = tempray;
end

% get skeletons of rays
skel = zeros(size(masks));
full_masks = zeros(size(masks));
ray_length = zeros(num_rays,1);
branches = cell(num_rays,1);
x_pos = cell(num_rays,1);
y_pos = cell(num_rays,1);
terminal_pts = cell(num_rays,1);
branching_pts = cell(num_rays,1);
root_pts = zeros(size(branches,1),1);
for i = 1:num_rays
    ray = masks(:,:,i);
    num_thresholds = 40;
    binary = zeros(size(ray,1), size(ray,2), num_thresholds);

    for j = 1:num_thresholds
        binary(:,:,j) = imbinarize(ray, 0.01*j);
    end
    figure; imagesc3D(binary);

    thresh = input("Which slide is the best contrast? ");
    close all;

    thresh = thresh/100;
    fin_ray = ray > thresh;

    full_fin_ray = bwconncomp(fin_ray);
    numOfPixels = cellfun(@numel,full_fin_ray.PixelIdxList);
    [~,indexOfMax] = max(numOfPixels);
    full_ray = zeros(size(fin_ray));
    full_ray(full_fin_ray.PixelIdxList{indexOfMax}) = 1;
    full_ray_out = full_ray;

    full_ray = imfill(full_ray,'holes');
    full_ray = imclose(full_ray,strel('disk',2));
    full_ray = imfill(full_ray,'holes');
    full_ray = imopen(full_ray,strel('disk',2));
    full_ray = imfill(full_ray,'holes');
    new_skeleton = bwskel(logical(full_ray));

    %Alternative splitting method to 'branchpoint'
    %Use convolution to identify points with more than 2 neighboring pixels
    filter = [1 1 1; 1 0 1; 1 1 1];

    new_skeleton_disconnect = new_skeleton & ~(new_skeleton & conv2(double(new_skeleton), filter, 'same')>2);

    cc = bwconncomp(new_skeleton_disconnect);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [sorted_px, ind] = sort(numPixels);

    %Remove components shorter than threshold
    threshold  = 40;
    for ii=ind(sorted_px<threshold)
        cur_comp = cc.PixelIdxList{ii};
        new_skeleton(cur_comp) = 0;

        %Before removing component, check whether image is still connected
        full_cc = bwconncomp(new_skeleton);
        if full_cc.NumObjects>1
            new_skeleton(cur_comp) = 1;
        end
    end
    new_skeleton = bwmorph(new_skeleton, 'spur');
    new_skeleton_disconnect = new_skeleton & ~(new_skeleton & conv2(double(new_skeleton), filter, 'same')>2);

    cc = bwconncomp(new_skeleton_disconnect);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [sorted_px, ind] = sort(numPixels);

    %Remove components shorter than threshold
    for ii=ind(sorted_px<threshold)
        cur_comp = cc.PixelIdxList{ii};
        new_skeleton(cur_comp) = 0;

        %Before removing component, check whether image is still connected
        full_cc = bwconncomp(new_skeleton);
        if full_cc.NumObjects>1
            new_skeleton(cur_comp) = 1;
        end
    end

    %Clean up left over spurs
    new_skeleton = bwmorph(new_skeleton,'skel',Inf);
    new_skeleton = bwmorph(new_skeleton, 'thin',Inf);
    new_skeleton = bwmorph(new_skeleton, 'spur');
    disp_new_skeleton = imdilate(new_skeleton,strel('disk',1));

    B = bwmorph(new_skeleton,'branchpoints');
    E = bwmorph(new_skeleton,'endpoints');

    [row, ~] = find(new_skeleton);
    [~, col] = find(new_skeleton(min(row),:));
    root_pixel = [min(row), min(col)];
    skel_mask = false(size(new_skeleton));
    skel_mask(root_pixel(1),root_pixel(2)) = true;
    D = bwdistgeodesic(logical(new_skeleton),skel_mask,'quasi-euclidean');
    D(isnan(D)) = 0;
    terminal_values = D.*E;
    branch_values = D.*B;
    branch_lengths = unique(branch_values).*pixelSize(1);
    branch_lengths = branch_lengths(2:end);
    ray_lengths = unique(terminal_values).*pixelSize(1);
    ray_lengths = ray_lengths(2:end);

    [x_endpts, y_endpts] = find(E);
    endpt_locs = [x_endpts, y_endpts];
    [x_branchpts, y_branchpts] = find(B);
    branchpt_locs = [x_branchpts, y_branchpts];
    distance = mean(ray_lengths);
    ray_lengths = ray_lengths(ray_lengths > 0.1*distance);
    distance = median(ray_lengths);

    new_skeleton_disconnect = new_skeleton & ~(new_skeleton & conv2(double(new_skeleton), filter, 'same')>2);

    all_points = B+E;
    all_points = imdilate(all_points,strel('disk',3));
    label_pts = bwlabeln(all_points);

    root_endpt = unique(label_pts.*skel_mask);
    root_endpt = root_endpt(end);

    term_endpts = unique(E.*label_pts);
    term_endpts = term_endpts(2:end);

    branch_pts = unique(B.*label_pts);
    branch_pts = branch_pts(2:end);

    if root_endpt == 0
        root_endpt = max(term_endpts);
    end

    branch_adj =  (label_pts.*new_skeleton_disconnect)+new_skeleton_disconnect;

    segment_props = regionprops(new_skeleton_disconnect,'Area','BoundingBox','PixelIdxList');

    adj_matrix = zeros(max(label_pts(:)));
    lengths = zeros(1,3);

    for j = 1:size(segment_props,1)

        segment = zeros(size(new_skeleton_disconnect));
        pixel_list = segment_props(j).PixelIdxList;
        segment(pixel_list) = 1;
        if size(pixel_list,1) > 100

            [row, ~] = find(segment);
            [~, col] = find(segment(min(row),:));
            new_root = [min(row), min(col)];
            seg_mask = false(size(segment));
            seg_mask(new_root(1),new_root(2)) = true;
            length = bwdistgeodesic(logical(segment),seg_mask,'quasi-euclidean');
            branch_length = max(length(:)).*pixelSize(1);

            adj_pair = unique(segment.*branch_adj);
            adj_pair = adj_pair(3:end)-1;

            lengths(j,:) = [adj_pair(1) adj_pair(2) branch_length];
            lengths = lengths(any(lengths,2),:);

        end
    end

    disp_new_skeleton = imdilate(disp_new_skeleton,strel('disk',1));
    skel(:,:,i) = disp_new_skeleton;
    full_masks(:,:,i) = full_ray_out;
    ray_length(i) = distance;
    branches{i} = lengths;

    index = [term_endpts; branch_pts]';
    perm_ind = zeros(1,size(index,2));
    for k = 1:size(index,2)
        perm_ind(k) = find(index == k);
    end

    x_vals = -[x_endpts; x_branchpts]';
    x_vals = x_vals(perm_ind);
    y_vals = [y_endpts; y_branchpts]';
    y_vals = y_vals(perm_ind);

    x_pos{i} = x_vals;
    y_pos{i} = y_vals;
    terminal_pts{i} = term_endpts;
    branching_pts{i} = branch_pts;
    root_pts(i) = root_endpt;
end

full_fin_masks = sum(full_masks,3);
full_skeleton = sum(skel,3);

max_nodes = [];
for i = 1:size(branches,1)
    max_nodes(i) = max(terminal_pts{i});
end

testx = [];
for i = 1:size(branches,1)
    testx = [testx x_pos{i}];
end

testy = [];
for i = 1:size(branches,1)
    testy = [testy y_pos{i}];
end

for i = 2:size(branches,1)
    branches{i}(:,1) = branches{i}(:,1) + sum(max_nodes(1:i-1));
    branches{i}(:,2) = branches{i}(:,2) + sum(max_nodes(1:i-1));
    terminal_pts{i} = terminal_pts{i} + sum(max_nodes(1:i-1));
    branching_pts{i} = branching_pts{i} + sum(max_nodes(1:i-1));
    root_pts(i) = root_pts(i) + sum(max_nodes(1:i-1));
end

test_term = [];
for i = 1:size(branches,1)
    test_term = [test_term; terminal_pts{i}];
end

test = [];
for i = 1:size(branches,1)
    test = [test; branches{i}];
end

figure;
subplot(1,3,1);
imshowpair(fin, full_fin_masks);

subplot(1,3,2);
imshowpair(fin, full_skeleton);

subplot(1,3,3);
ray_graph = graph(test(:,1), test(:,2), test(:,3));
%line_widths = 5*ray_graph.Edges.Weight./max(ray_graph.Edges.Weight);
%h = plot(ray_graph, 'EdgeLabel', ray_graph.Edges.Weight, 'LineWidth',line_widths,'layout','force');
h = plot(ray_graph, 'EdgeLabel', ray_graph.Edges.Weight);
highlight(h, test_term, 'NodeColor','r');
highlight(h, root_pts, 'NodeColor','g');
h.EdgeFontSize = 11;

h.XData = testy';
h.YData = testx';