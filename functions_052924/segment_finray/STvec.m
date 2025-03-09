function [e2_out] = STvec(I_in,mask_in,options)
arguments
    I_in;
    mask_in;
    options.bd = 6; % smoothening bandwidth for gradient
    options.scaling = 0.1;
end
% Structure Tensor vector field
%   Detailed explanation goes here
%%
% initialize
I = imresize(I_in,options.scaling);
% mask = imresize(mask_in,options.scaling);
% figure; imshow(I)
n = size(I,1);
m = size(I,2);

sigma = 1;
cutoff = ceil(3.5*sigma);
rw = options.bd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gradients
h = fspecial('gaussian',[1,2*cutoff+1],sigma);
dh = h .* (-cutoff:cutoff) / (-sigma^2);
gx = conv2(dh,h,I,'same');

h = fspecial('gaussian',[2*cutoff+1,1],sigma);
dh = h .* (-cutoff:cutoff)' / (-sigma^2);
gy= conv2(h,dh,I,'same');

% figure; imshow(gx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute tensors

gx2 = gx.^2;
gy2 = gy.^2;
gxy = gx.*gy;

% smooth
gx2_sm = imgaussfilt(gx2,rw); %rw/sqrt(2*log(2))
gy2_sm = imgaussfilt(gy2,rw);
gxy_sm = imgaussfilt(gxy,rw);
H = zeros(n,m,2,2);
H(:,:,1,1) = gx2_sm; 
H(:,:,2,2) = gy2_sm; 
H(:,:,1,2) = gxy_sm; 
H(:,:,2,1) = gxy_sm; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigen decomposition
l1 = zeros(n,m);
l2 = zeros(n,m);
e1 = zeros(n,m,2);
e2 = zeros(n,m,2);
for i = 1:n
    for j = 1:m
        Hmat = zeros(2);
        Hmat(:,:) = H(i,j,:,:);
        [V,D] = eigs(Hmat);
        D = abs(D);
        l1(i,j) = D(1,1); % eigen values
        l2(i,j) = D(2,2); 
        e1(i,j,:) = V(:,1); % eigen vectors
        e2(i,j,:) = V(:,2); 
    end
end

%%
% figure;imshow(I)
% figure;imshow(abs(l1./l2)>12,[]);
%% correct directions
fy = e2(:,:,1);
fx = e2(:,:,2);
flip_logical = fy>0;
fx(flip_logical) = -fx(flip_logical);
fy(flip_logical) = -fy(flip_logical);
e2 = cat(3,fy,fx);
%% mask
e2_large = imresize(e2,size(I_in));

fx = e2_large(:,:,2);
fx(~mask_in) = 0;
fy = e2_large(:,:,1);
fy(~mask_in) = 0;

e2_out = cat(3,fy,fx);

    %% view vector field
%     figure; imshow(I,[])
%     hold on;
%     xx = 1:1:size(I,2);
%     yy = 1:1:size(I,1);
%     u = e2(yy,xx,2);
%     v = e2(yy,xx,1);
%     % quiver(xx,yy,u,v,LineWidth=1,AutoScaleFactor=2);
%     quiver(xx,yy,u,v);
end

