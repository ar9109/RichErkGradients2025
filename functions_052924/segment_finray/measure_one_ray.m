function [arcBW,arcLHere,end_bot,BW_oneray_small,BW_oneray_small_closed,row_excerpt,col_excerpt] = measure_one_ray(BW_oneray,options)
arguments
    BW_oneray;
    options.verbose = [];
    options.initial_seeds = 0;
end
%MEASURE_ONE_RAY Summary of this function goes here
%   Detailed explanation goes here

%% crop canvas
% outer box
% [row,col] = ind2sub(size(BW_oneray),CC.PixelIdxList{raynum});
[row,col] = find(BW_oneray);

row_excerpt = max(min(row)-10,1):min(max(row)+100,size(BW_oneray,1));
col_excerpt = max(min(col)-100,1):min(max(col)+100,size(BW_oneray,2));

BW_oneray_small = BW_oneray(row_excerpt,col_excerpt);
if options.initial_seeds
    cut_off = find(diff(sum(BW_oneray_small,2))>0,1,'last');
    BW_oneray_small(cut_off:end,:)=0;
end
% figure; imshow(BW_oneray_small);


%% measure length for one ray
% fill the gap between branches
BW_oneray_small_closed = imclose(BW_oneray_small,strel('disk',100,8));

% fit with a spline
[sprow,spcol] = find(BW_oneray_small_closed);
fit_obj = fit(sprow,spcol,'smoothingspline', SmoothingParam = 1E-7);

% sum of pixels in each row - for determining unwanted bend
sum_row = sum(BW_oneray_small_closed,2);
x_sum = find(sum_row);
y_sum = sum_row(x_sum);
fit_obj_sum = fit(x_sum,y_sum,'smoothingspline', SmoothingParam = 1E-5);



% figure; imshow(BW_oneray_small_closed);
% hold on; plot(end_top(2),end_top(1),'ro'); hold off;
%% plot 2nd derivative


BD = cell2mat(bwboundaries(BW_oneray_small_closed));
x_fit = min(BD(:,1)):max(BD(:,1));
y_fit = fit_obj(x_fit)';
y_fit_sum = fit_obj_sum(x_fit)';
der_2nd = -(del2(y_fit_sum));
der_2nd_bend = abs(del2(y_fit));
[peak_val,peak_idx] = max(der_2nd(1:min(numel(der_2nd(1:end-150)),400))); % take max of laplacian only on the first 400 pixels, and not the last 150 pixels
if isempty(peak_val) || options.initial_seeds
    peak_idx = 1;
else 
    peak_val_bend = max(der_2nd_bend(max(peak_idx-50,1):min(peak_idx+50,numel(der_2nd_bend))));
    if peak_val<20/1e4 || peak_val_bend<10/1e4 || (peak_idx>300&&peak_val<50/1e4)% only highly bended curve is corrected
        peak_idx = 1;
    end
end
if peak_idx+100<=numel(x_fit)
    y_extrap = interp1(x_fit(peak_idx+40:end),y_fit(peak_idx+40:end),x_fit(1:peak_idx+39),'linear','extrap'); % shift to the right(or bottom) 40 pixels and extrapolate
    y_fit2 = [y_extrap, y_fit(peak_idx+40:end)];
    y_fit2(y_fit2<=0) = 1; % usually when this happens this ray is wrong
    y_fit2(y_fit2>size(BW_oneray_small_closed,2)) = size(BW_oneray_small_closed,2);
else
    y_fit2 = y_fit;
end
%% debug look at the ray arc
if ~isempty(options.verbose)&&options.verbose
    figure;
    plot(x_fit,abs(del2(y_fit)).*1e4); 
    hold on;
    plot(x_fit,y_fit); 
    plot(x_fit(peak_idx+40),y_fit(peak_idx+40),'r*');hold off;
    figure;
    plot(x_fit,-(del2(y_fit_sum)).*1e4); 
    hold on;
    plot(x_fit,y_fit_sum); 
    plot(x_fit(peak_idx+40),y_fit_sum(peak_idx+40),'r*');hold off;
end
%% find intersection between spline and boundary of the mask
BD_midpoint = (min(BD(:,1))+max(BD(:,1)))/2;
[~,yidx] = ismember(BD(:,1),x_fit);
BD_spline_dist = abs(reshape(y_fit2(yidx(:)),[],1)-BD(:,2));
dist_top = BD_spline_dist;
dist_top(BD(:,1)>BD_midpoint)=nan;
dist_bot = BD_spline_dist;
dist_bot(BD(:,1)<=BD_midpoint)=nan;
[~,idx_top] = min(dist_top);
[~,idx_bot] = min(dist_bot);
end_top = BD(idx_top,:);
end_bot = BD(idx_bot,:);

x_fit_final = end_top(1):end_bot(1);
[~,yidx] = ismember(x_fit_final,x_fit);
y_fit_final = y_fit2(yidx(:));
%% calculate arc length
% x_fit = linspace(end_top(1),end_bot(1),50);
% x_fit = end_top(1):end_bot(1);
% y_fit = fit_obj(x_fit)';
% arcLHere = trapz(x_fit,sqrt(1+(ppval(fnder(fit_obj.p),x_fit)).^2)); % method1
arcLHere = sum(hypot(diff(x_fit_final(:)),diff(y_fit_final(:)))); % method 2

%% the overlaid image
arcBW = false(size(BW_oneray_small_closed));
arcBW(sub2ind(size(arcBW),round(x_fit_final),round(y_fit_final))) = true;
arcBW = imdilate(arcBW,strel('disk',5,8)); %.*BW_oneray_small_closed;
end

