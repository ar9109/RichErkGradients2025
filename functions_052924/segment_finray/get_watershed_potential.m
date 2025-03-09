function [potential] = get_watershed_potential(I_in,mask_in,end_bot)
%GET_WATERSHED_POTENTIAL Summary of this function goes here
%   Detailed explanation goes here
%% dip_image, not so useful
% run('C:\Program Files\DIPimage 2.9\dipstart.m')
% %%
% sigma1 = 3;
% rw = 3;
% I2 = dip_image(I);
% g = gradient(I2, sigma1);
% H = gaussf(g*g.', rw);
% [e,l] = eig(H);
% % Equivalences with your outputs:
% l1 = l{2};
% l2 = l{1};
% e1 = e{2,:};
% e2 = e{1,:};
% 
% tt = dip_array(l1);
% figure; imshow(tt,[]);
% hold on;
% quiver(1:10:size(I,2),1:10:size(I,1),e1)
%%
% initialize
I = imresize(I_in,0.1);
mask = imresize(mask_in,0.1);
% figure; imshow(I)
n = size(I,1);
m = size(I,2);

sigma = 1;
cutoff = ceil(3.5*sigma);
rw = 3;

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
%% view vector field
figure; imshow(I,[])
hold on;
xx = 1:size(I,2);
yy = 1:size(I,1);
u = e2(yy,xx,2);
v = e2(yy,xx,1);
quiver(xx,yy,u,v);
% quiver(x,y,fx_corrected,fy_corrected);
plot(end_bot(2)./10,end_bot(1)./10, 'ro')
plot(y(:,1),y(:,2), 'r-')

hold off;
% 
% figure;imshow(mask)
%% solve trajectory
x0 = end_bot(2).*0.1; y0 = end_bot(1).*0.1;

odefun = @(t,y)match_eigvec(t,y,e2);
[t,y] = ode45(odefun,[0,1000],[x0;y0],odeset('MaxStep',1e-1));




%%

%% invererse gradient
fx = e2(:,:,2);
fx(~mask) = 0;
fy = e2(:,:,1);
fy(~mask) = 0;
% fy_corrected = fy;
% fy_corrected(:) = -1;
% fx_corrected = fy_corrected./fy.*fx;
% tic;
fhat = intgrad2(fx,fy,1,1,0);
% toc
%   Time required was 4.332 seconds

figure;surf(fhat)
potential = imresize(fhat,size(I_in),'bilinear');
end

