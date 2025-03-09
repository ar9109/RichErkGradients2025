function [length_mat,imQCBFHere,imQCosxHere,tt] = update_raylengths_STvec(length_mat,fiducial_fullimage,BF,osx,fin_bw,fish,s,options)
arguments
    length_mat;
    fiducial_fullimage;
    BF;
    osx;
    fin_bw;
    fish;
    s;
    options.verbose = [];
    options.rays_to_loop = [];
end
%UPDATE_RAYLENGTHS Summary of this function goes here
%   Detailed explanation goes here

%% calculate vector field

% tic;
% % mask_in = imclose(PB_bw_trunc_dilated,strel('disk',50));
% mask_in = imclose(PB_bw,strel('disk',100,8));
% % mask_in = imdilate(mask_in,strel('disk',50));
mask_in = fin_bw;
% I = imadjust(ITT,[0.03,0.07]);
I = osx;
[e2] = STvec(I,mask_in,scaling = 0.2, bd = 12);
% toc

%% initialize functions
odefun = @(t,y)ode_match_eigvec(t,y,e2);
myEvent = @(T,Y)odeEvent_stop_at_boundary(T,Y,mask_in);


%% iterate through rays and solve trajectory
QC_arcBW = uint8(zeros(size(osx)));
tt = [];
if isempty(options.rays_to_loop)
    rays_to_loop = 1:numel(fiducial_fullimage);
else
    rays_to_loop = options.rays_to_loop;
end
for raynum = rays_to_loop
    end_bot = fiducial_fullimage{raynum};
    if isempty(end_bot)
        continue
    end
    x0 = end_bot(2); y0 = end_bot(1);
    [t,y] = ode45(odefun,[0,100000],[x0;y0],odeset('MaxStep',1,'Events', myEvent));


    %% view vector field
%     figure; imshow(I,[])
%     hold on;
%     xx = 1:5:size(I,2);
%     yy = 1:5:size(I,1);
%     u = e2(yy,xx,2);
%     v = e2(yy,xx,1);
%     % quiver(xx,yy,u,v,LineWidth=1,AutoScaleFactor=2);
%     quiver(xx,yy,u,v);
%     
%     % quiver(x,y,fx_corrected,fy_corrected);
%     plot(end_bot(2),end_bot(1), 'ro')
%     plot(y(:,1),y(:,2), 'r-',LineWidth=2)
%     
%     hold off;
    % 
    % figure;imshow(mask)
    
    % %% try look at bright field
    % figure;
    % sliceViewer(e2)
    % %%
    % ITT = loadtiff('E:\21Jul23_osx_H2AmCherry_calcineurin_dissecting\raw\0dpa\fish1_12x_0dpa_post-Image Export-28_c2.tif');
    % ITT = ITT(:,:,1);
    
    %% measure length
    arcLHere = sum(hypot(diff(y(:,1)),diff(y(:,2))));
    
    % paint 
    arcIdx = unique(round(y),'rows');
    arcBW = false(size(osx));
    arcBW(sub2ind(size(osx),round(arcIdx(:,2)),round(arcIdx(:,1)))) = true;
    arcBW = imdilate(arcBW,strel('disk',5,8));
    QC_arcBW(arcBW) = raynum; % raynum

    %% save length
    thisray = [length_mat.fish]==fish&[length_mat.ray]==raynum; % better to use structure in this code
    length_mat(thisray).arcL(s) = arcLHere;
    length_mat(thisray).fiducial(s,:) = fiducial_fullimage{raynum};
    tt = [tt arcLHere];

end

imQCosxHere = labeloverlay(osx,QC_arcBW,Colormap=lines(30),Transparency=0);
if isempty(BF)
    imQCBFHere = [];
else
    imQCBFHere = labeloverlay(BF,QC_arcBW,Colormap=lines(30),Transparency=0);
end
end

