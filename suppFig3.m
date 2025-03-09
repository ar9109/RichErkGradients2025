%% Supp. Fig. 3B -left

%load "analysis_mat_WTdataSet.mat" from "SuppFig3" folder

% Plot binned average ERK v. position - Used 1 May 24 - color Blue



ERKall = [];
for fish = 4
%for fish = 1:55
    fishNow = fish;
    
    ERKcollect = [];
    xValuesCollect = [];
    xValuesMean = [];
    ERKcollectMean = [];

    %for ray = 2:10
    for ray = 2:2
    %for ray = 7:7
        rayNow = ray;
        
        %for hpp = 72:336
            %hppHere = hpp;
            hppMin = 72;
            %figure;
            
            for i = 1:size(analysis_mat,2)
                if analysis_mat(i).fish == fishNow && analysis_mat(i).ray == rayNow &&  analysis_mat(i).hpa >= hppMin

                    fishHere = analysis_mat(i).fish;
                    rayHere = analysis_mat(i).ray;
                    hppHere = analysis_mat(i).hpa;
                    ktrHere = analysis_mat(i).ktr-0.8;
                    ktrHere = ktrHere.*(ktrHere>0);
                    ccHerepre = analysis_mat(i).ccrot(:,:);
                    endPointHere = 0;
                    LregHere = analysis_mat(i).L_reg;
                    ampPointHere = LregHere;
                    LampHere = analysis_mat(i).L_amp;

                    idxs=~isnan(ktrHere) & ~isinf(ktrHere); %& ccHerepre(:,1)>ampPointHere% &ccHerepre<endPointHere; %All the FUCCI values that are not NaN
                    ccHere=ccHerepre(idxs,:);%
                    ktrHere=ktrHere(idxs);%

                    ccIdx = ccHere(:,1) < ampPointHere;

                    ccHereClean = ccHere(ccIdx,:);
                    ktrHereClean = ktrHere(ccIdx);

                    xValues = LregHere - ccHereClean(:,1);

                    %bin data
                    xData = xValues;
                    yData = ktrHereClean;
                    removeNaN=~isnan(xData)&~isnan(yData);
                    xData = xData(removeNaN);xData = xData(:);
                    yData = yData(removeNaN);yData = yData(:);
                    xBins = 0:250:max(xValues)+250;

                    %xBins = linspace(min_x,max_x,nbins+1);xBins = xBins(:);
                    idxBins=discretize(xData,xBins);
                    idxBins = idxBins(~isnan(idxBins));
                    xMean = [(xBins(1:end-1)+xBins(2:end))/2];
                    yMean = accumarray(idxBins,yData,size(xMean'),@nanmean,NaN);
                    yStd = accumarray(idxBins,yData,size(xMean'),@nanstd,NaN);
                    yN = accumarray(idxBins,yData,size(xMean'),@(x) sum(~isnan(x)),NaN);
                    ySEM = yStd./sqrt(yN);
                    %end bin data
                    
                    BaseName='Fish';
                    BaseName2='Ray';
                    BaseName3='hpa';
                    FileName1=[BaseName,num2str(fishHere)];
                    FileName2=[BaseName2,num2str(rayHere)];
                    FileName3=[num2str(hppHere),BaseName3];

%Set up colors
    tLevels = 10;
    color_map=colormap(cbrewer2('Blues',tLevels));
    %colormap(jet(tLevels));
    %colors = jet(tLevels);
    tmin = 108;
    tmax = 192;
    
%     if rayHere < 4
        LS = '-';
%     else
%         LS = '--';
%     end

    idxColor = round((tLevels-1)*(hppHere-tmin)/(tmax-tmin)+1);
%    plot(xValues,averageKTR./avgKTRray2,'.-','color',colors(idxColor,:),'LineWidth',3,'MarkerSize',30); hold on;
    plot(xMean,yMean','.-','color',color_map(idxColor,:),'LineWidth',3,'MarkerSize',30,'LineStyle',LS); hold on;
    set(gca, 'fontsize', 20);
    ylim([0 0.6]);
    xlim([0 2000]);
    xlabel('Position (\mum)');
    ylabel('ERK Activity (A.U.)');
    %title([FileName1 ' ' FileName2]);
    
    c= colorbar('Ticks', 108:12:192,'TickLabels',{'','120','132','144','156','168','180',''});
    caxis([tmin tmax])
    c.FontSize = 16;
    c.Label.String = 'Time (hpa)';
    set( c, 'YDir', 'reverse' );
    set(gca, 'fontsize', 16);
    %set(gca, 'fontsize', 20, 'XDir','reverse');
    %print('-djpeg',[paths.plotFolder FileName1 '_' FileName2 '_' FileName3 '_ERKvPosition_wOffset_CollectAverageVposition']);    
                end
            end

    end
end

%% Supp. Fig. 3B -right

%load "analysis_mat_WTdataSet.mat" from "SuppFig3" folder

% Plot binned average ERK v. position - Used 1 May 24 - color Blue



ERKall = [];
for fish = 4
%for fish = 1:55
    fishNow = fish;
    
    ERKcollect = [];
    xValuesCollect = [];
    xValuesMean = [];
    ERKcollectMean = [];

    %for ray = 2:10
    %for ray = 2:2
    for ray = 7:7
        rayNow = ray;
        
        %for hpp = 72:336
            %hppHere = hpp;
            hppMin = 72;
            %figure;
            
            for i = 1:size(analysis_mat,2)
                if analysis_mat(i).fish == fishNow && analysis_mat(i).ray == rayNow &&  analysis_mat(i).hpa >= hppMin

                    fishHere = analysis_mat(i).fish;
                    rayHere = analysis_mat(i).ray;
                    hppHere = analysis_mat(i).hpa;
                    ktrHere = analysis_mat(i).ktr-0.8;
                    ktrHere = ktrHere.*(ktrHere>0);
                    ccHerepre = analysis_mat(i).ccrot(:,:);
                    endPointHere = 0;
                    LregHere = analysis_mat(i).L_reg;
                    ampPointHere = LregHere;
                    LampHere = analysis_mat(i).L_amp;

                    idxs=~isnan(ktrHere) & ~isinf(ktrHere); %& ccHerepre(:,1)>ampPointHere% &ccHerepre<endPointHere; %All the FUCCI values that are not NaN
                    ccHere=ccHerepre(idxs,:);%
                    ktrHere=ktrHere(idxs);%

                    ccIdx = ccHere(:,1) < ampPointHere;

                    ccHereClean = ccHere(ccIdx,:);
                    ktrHereClean = ktrHere(ccIdx);

                    xValues = LregHere - ccHereClean(:,1);

                    %bin data
                    xData = xValues;
                    yData = ktrHereClean;
                    removeNaN=~isnan(xData)&~isnan(yData);
                    xData = xData(removeNaN);xData = xData(:);
                    yData = yData(removeNaN);yData = yData(:);
                    xBins = 0:250:max(xValues)+250;

                    %xBins = linspace(min_x,max_x,nbins+1);xBins = xBins(:);
                    idxBins=discretize(xData,xBins);
                    idxBins = idxBins(~isnan(idxBins));
                    xMean = [(xBins(1:end-1)+xBins(2:end))/2];
                    yMean = accumarray(idxBins,yData,size(xMean'),@nanmean,NaN);
                    yStd = accumarray(idxBins,yData,size(xMean'),@nanstd,NaN);
                    yN = accumarray(idxBins,yData,size(xMean'),@(x) sum(~isnan(x)),NaN);
                    ySEM = yStd./sqrt(yN);
                    %end bin data
                    
                    BaseName='Fish';
                    BaseName2='Ray';
                    BaseName3='hpa';
                    FileName1=[BaseName,num2str(fishHere)];
                    FileName2=[BaseName2,num2str(rayHere)];
                    FileName3=[num2str(hppHere),BaseName3];

%Set up colors
    tLevels = 10;
    color_map=colormap(cbrewer2('Blues',tLevels));
    %colormap(jet(tLevels));
    %colors = jet(tLevels);
    tmin = 108;
    tmax = 192;
    
%     if rayHere < 4
        LS = '-';
%     else
%         LS = '--';
%     end

    idxColor = round((tLevels-1)*(hppHere-tmin)/(tmax-tmin)+1);
%    plot(xValues,averageKTR./avgKTRray2,'.-','color',colors(idxColor,:),'LineWidth',3,'MarkerSize',30); hold on;
    plot(xMean,yMean','.-','color',color_map(idxColor,:),'LineWidth',3,'MarkerSize',30,'LineStyle',LS); hold on;
    set(gca, 'fontsize', 20);
    ylim([0 0.6]);
    xlim([0 2000]);
    xlabel('Position (\mum)');
    ylabel('ERK Activity (A.U.)');
    %title([FileName1 ' ' FileName2]);
    
    c= colorbar('Ticks', 108:12:192,'TickLabels',{'','120','132','144','156','168','180',''});
    caxis([tmin tmax])
    c.FontSize = 16;
    c.Label.String = 'Time (hpa)';
    set( c, 'YDir', 'reverse' );
    set(gca, 'fontsize', 16);
    %set(gca, 'fontsize', 20, 'XDir','reverse');
    %print('-djpeg',[paths.plotFolder FileName1 '_' FileName2 '_' FileName3 '_ERKvPosition_wOffset_CollectAverageVposition']);    
                end
            end

    end
end


%% Supp. Fig. 3C -left
% Plot Normalized binned average ERK v. NORMALIZED position - Used 26 Apr 23

%load "analysis_mat_WTdataSet.mat" from "SuppFig3" folder

close all;


ERKall = [];
for fish = 4
%for fish = 1:55
    fishNow = fish;
    
    ERKcollect = [];
    xValuesCollect = [];
    xValuesMean = [];
    ERKcollectMean = [];

    %for ray = 2:10
    for ray = 2:2
    %for ray = 7:7
        rayNow = ray;
        
        %for hpp = 72:336
            %hppHere = hpp;
            hppMin = 72;
            %figure;
            
            for i = 1:size(analysis_mat,2)
                if analysis_mat(i).fish == fishNow && analysis_mat(i).ray == rayNow &&  analysis_mat(i).hpa >= hppMin

                    fishHere = analysis_mat(i).fish;
                    rayHere = analysis_mat(i).ray;
                    hppHere = analysis_mat(i).hpa;
                    ktrHere = analysis_mat(i).ktr-0.8;
                    ktrHere = ktrHere.*(ktrHere>0);
                    ccHerepre = analysis_mat(i).ccrot(:,:);
                    endPointHere = 0;
                    LregHere = analysis_mat(i).L_reg;
                    ampPointHere = LregHere;
                    LampHere = analysis_mat(i).L_amp;


                    idxs=~isnan(ktrHere) & ~isinf(ktrHere); %& ccHerepre(:,1)>ampPointHere% &ccHerepre<endPointHere; %All the FUCCI values that are not NaN
                    ccHere=ccHerepre(idxs,:);%
                    ktrHere=ktrHere(idxs);%

                    ccIdx = ccHere(:,1) < ampPointHere;

                    ccHereClean = ccHere(ccIdx,:);
                    ktrHereClean = ktrHere(ccIdx);

                    xValues = LregHere - ccHereClean(:,1);
                    xValues = xValues./max(xValues);

                    %bin data
                    xData = xValues;
                    yData = ktrHereClean;
                    removeNaN=~isnan(xData)&~isnan(yData);
                    xData = xData(removeNaN);xData = xData(:);
                    yData = yData(removeNaN);yData = yData(:);
                    %xBins = 0:250:max(xValues)+250;
                    xBins = 0:0.1:1;

                    %xBins = linspace(min_x,max_x,nbins+1);xBins = xBins(:);
                    idxBins=discretize(xData,xBins);
                    idxBins = idxBins(~isnan(idxBins));
                    xMean = [(xBins(1:end-1)+xBins(2:end))/2];
                    yMean = accumarray(idxBins,yData,size(xMean'),@nanmean,NaN);
                    yStd = accumarray(idxBins,yData,size(xMean'),@nanstd,NaN);
                    yN = accumarray(idxBins,yData,size(xMean'),@(x) sum(~isnan(x)),NaN);
                    ySEM = yStd./sqrt(yN);
                    %end bin data
                    
                    BaseName='Fish';
                    BaseName2='Ray';
                    BaseName3='hpa';
                    FileName1=[BaseName,num2str(fishHere)];
                    FileName2=[BaseName2,num2str(rayHere)];
                    FileName3=[num2str(hppHere),BaseName3];

%Set up colors
    tLevels = 10;
    color_map=colormap(cbrewer2('Blues',tLevels));
    %colormap(jet(tLevels));
    %colors = jet(tLevels);
    tmin = 108;
    tmax = 192;
    
%     if rayHere < 4
        LS = '-';
%     else
%         LS = '--';
%     end

    idxColor = round((tLevels-1)*(hppHere-tmin)/(tmax-tmin)+1);
%    plot(xValues,averageKTR./avgKTRray2,'.-','color',colors(idxColor,:),'LineWidth',3,'MarkerSize',30); hold on;
    plot(xMean,yMean./nanmean(yMean)','.-','color',color_map(idxColor,:),'LineWidth',3,'MarkerSize',30,'LineStyle',LS); hold on;
    set(gca, 'fontsize', 16);
    ylim([0 2.5]);
    xlim([0 1]);
    xlabel('Normalized Position ({\itu})');
    ylabel('ERK Activity / Average ERK Activity');
    %title([FileName1 ' ' FileName2]);
    
    c= colorbar('Ticks', 108:12:192,'TickLabels',{'','120','132','144','156','168','180',''});
    caxis([tmin tmax])
    c.FontSize = 16;
    c.Label.String = 'Time (hpa)';
    set( c, 'YDir', 'reverse' );
    set(gca, 'fontsize', 16);
    %set(gca, 'fontsize', 20, 'XDir','reverse');
    %print('-djpeg',[paths.plotFolder FileName1 '_' FileName2 '_' FileName3 '_ERKvPosition_wOffset_CollectAverageVposition']);    
                end
            end

    end
end

%% Supp. Fig. 3C -right
% Plot Normalized binned average ERK v. NORMALIZED position - Used 26 Apr 23

%load "analysis_mat_WTdataSet.mat" from "SuppFig3" folder

close all;


ERKall = [];
for fish = 4
%for fish = 1:55
    fishNow = fish;
    
    ERKcollect = [];
    xValuesCollect = [];
    xValuesMean = [];
    ERKcollectMean = [];

    %for ray = 2:10
    %for ray = 2:2
    for ray = 7:7
        rayNow = ray;
        
        %for hpp = 72:336
            %hppHere = hpp;
            hppMin = 72;
            %figure;
            
            for i = 1:size(analysis_mat,2)
                if analysis_mat(i).fish == fishNow && analysis_mat(i).ray == rayNow &&  analysis_mat(i).hpa >= hppMin

                    fishHere = analysis_mat(i).fish;
                    rayHere = analysis_mat(i).ray;
                    hppHere = analysis_mat(i).hpa;
                    ktrHere = analysis_mat(i).ktr-0.8;
                    ktrHere = ktrHere.*(ktrHere>0);
                    ccHerepre = analysis_mat(i).ccrot(:,:);
                    endPointHere = 0;
                    LregHere = analysis_mat(i).L_reg;
                    ampPointHere = LregHere;
                    LampHere = analysis_mat(i).L_amp;


                    idxs=~isnan(ktrHere) & ~isinf(ktrHere); %& ccHerepre(:,1)>ampPointHere% &ccHerepre<endPointHere; %All the FUCCI values that are not NaN
                    ccHere=ccHerepre(idxs,:);%
                    ktrHere=ktrHere(idxs);%

                    ccIdx = ccHere(:,1) < ampPointHere;

                    ccHereClean = ccHere(ccIdx,:);
                    ktrHereClean = ktrHere(ccIdx);

                    xValues = LregHere - ccHereClean(:,1);
                    xValues = xValues./max(xValues);

                    %bin data
                    xData = xValues;
                    yData = ktrHereClean;
                    removeNaN=~isnan(xData)&~isnan(yData);
                    xData = xData(removeNaN);xData = xData(:);
                    yData = yData(removeNaN);yData = yData(:);
                    %xBins = 0:250:max(xValues)+250;
                    xBins = 0:0.1:1;

                    %xBins = linspace(min_x,max_x,nbins+1);xBins = xBins(:);
                    idxBins=discretize(xData,xBins);
                    idxBins = idxBins(~isnan(idxBins));
                    xMean = [(xBins(1:end-1)+xBins(2:end))/2];
                    yMean = accumarray(idxBins,yData,size(xMean'),@nanmean,NaN);
                    yStd = accumarray(idxBins,yData,size(xMean'),@nanstd,NaN);
                    yN = accumarray(idxBins,yData,size(xMean'),@(x) sum(~isnan(x)),NaN);
                    ySEM = yStd./sqrt(yN);
                    %end bin data
                    
                    BaseName='Fish';
                    BaseName2='Ray';
                    BaseName3='hpa';
                    FileName1=[BaseName,num2str(fishHere)];
                    FileName2=[BaseName2,num2str(rayHere)];
                    FileName3=[num2str(hppHere),BaseName3];

%Set up colors
    tLevels = 10;
    color_map=colormap(cbrewer2('Blues',tLevels));
    %colormap(jet(tLevels));
    %colors = jet(tLevels);
    tmin = 108;
    tmax = 192;
    
%     if rayHere < 4
        LS = '-';
%     else
%         LS = '--';
%     end

    idxColor = round((tLevels-1)*(hppHere-tmin)/(tmax-tmin)+1);
%    plot(xValues,averageKTR./avgKTRray2,'.-','color',colors(idxColor,:),'LineWidth',3,'MarkerSize',30); hold on;
    plot(xMean,yMean./nanmean(yMean)','.-','color',color_map(idxColor,:),'LineWidth',3,'MarkerSize',30,'LineStyle',LS); hold on;
    set(gca, 'fontsize', 16);
    ylim([0 2.5]);
    xlim([0 1]);
    xlabel('Normalized Position ({\itu})');
    ylabel('ERK Activity / Average ERK Activity');
    %title([FileName1 ' ' FileName2]);
    
    c= colorbar('Ticks', 108:12:192,'TickLabels',{'','120','132','144','156','168','180',''});
    caxis([tmin tmax])
    c.FontSize = 16;
    c.Label.String = 'Time (hpa)';
    set( c, 'YDir', 'reverse' );
    set(gca, 'fontsize', 16);
    %set(gca, 'fontsize', 20, 'XDir','reverse');
    %print('-djpeg',[paths.plotFolder FileName1 '_' FileName2 '_' FileName3 '_ERKvPosition_wOffset_CollectAverageVposition']);    
                end
            end

    end
end

%%  Supp. Fig. 3D
% f(u) on timeave matrix, aggregate single nuclei first - subtract later - color by Lamp

%load "analysis_mat_timeaverage_042624.mat" from "SuppFig3" folder
E0 = 0.8;


L_reg = [analysis_mat_timeaverage.L_reg_timeave];
L_amp = [analysis_mat_timeaverage.L_amp];
phi = L_reg./L_amp;



% bin average for single ray-time
averageKTR_u_cell = arrayfun(@(s)...
    adjust_neg(bin_average(s.u,...
    s.ktr,10,0,1)-E0),... %changed from 10 to 20
    analysis_mat_timeaverage,'UniformOutput',false);

mean_ktr_bin10_offsetE0 = cellfun(@nanmean,averageKTR_u_cell);



% normalize, e = Erk(u), m = mean(Erk)
norm_averageKTR_u_cell = arrayfun(@(e,m)(e{:})./(m),averageKTR_u_cell,mean_ktr_bin10_offsetE0,'UniformOutput',false);

norm_averageKTR_u = cell2mat(norm_averageKTR_u_cell);

%
% fit not-timeaveraged data with a line
% bin average
data.x = repmat([0.05:0.1:0.95]',[1,size(norm_averageKTR_u,2)]);
data.y = norm_averageKTR_u;
[binned.y,binned.x,binned.sem] = bin_average(data.x(:),data.y(:),10,0,1);
fit_x = linspace(0,1,100);


% fit
xdata = data.x(~isnan(data.y));
ydata = data.y(~isnan(data.y));
ft = fittype({'-x','1'});
foKTR = fitoptions('Method','LinearLeastSquares',...
    'Lower',[-Inf,0],...
    'Weights', []);
[fit_obj,fit_gof] = fit(xdata(:),ydata(:),ft,foKTR);

% plot
tLevels = 101;
color_map=colormap(cbrewer2('Reds',tLevels));

tt=[];
color_res = 100;
color_code = phi;
color_code = L_amp;

c_digit = -2;
cmax = round(max(color_code),c_digit);
cmin = round(min(color_code),c_digit);
color_unit = (cmax-cmin)./color_res; % config color

f = figure('Visible','on');hold on;

cm = colormap(flipud(cool(round(color_res)+1))); % config color


hold on;
for i = numel(analysis_mat_timeaverage):-1:1
%     s = analysis_mat_timeaverage(i);
%     sz = min([numel(s.averageKTR),numel(s.fractionGEM)]);

    cl_idx = round((color_code(i)-cmin)./color_unit)+1;
    if cl_idx <= 0
        cl_idx = 1;
    elseif cl_idx > round(color_res)+1
        cl_idx = round(color_res)+1;
    end
    cl = cm(cl_idx,:); % config color

    x_plot = data.x(:,i);
    y_plot = data.y(:,i);

    x_plot = x_plot(~isnan(y_plot));
    y_plot = y_plot(~isnan(y_plot));
%     if numel(y_plot)<10
%         tt = [tt,i];
%     end
    colorCode = round(((L_amp(i)-700)./3353).*100+1);
    plot(x_plot,y_plot,'-','color',[color_map(colorCode,:),0.5],'markersize',5,'LineWidth',1.5);
    

end
hold off;
colormap(color_map);
%c = colorbar('Ticks',0:0.5:1,'TickLabels', round((linspace(cmin,cmax,3)),c_digit));
c = colorbar;
set(c,'TickLabels',{'700','1380','2060','2740','3420','4100'});
c.Label.String = 'Length Amputated (um)';

hold on;
%plot([0.05:0.1:0.95]',mean(data.y,2),'-','color','k','markersize',5,'LineWidth',1.5);
%hold off;

fit_gof.rsquare = round(fit_gof.rsquare,1);

 plot_fit(fit_obj,fit_x,fit_gof,'confidence',false);
ylim([0,2.5])
xlabel('Normalized Space ({\itu})')
ylabel('Normalized ERK')
% title(['fit without time-average, hpa', num2str(t_plot)], ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b-fit_obj.a)])
%title(['fit rolling window timeave'], ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b-fit_obj.a)])

set(gca,'XDir','reverse')
config_plot(f,c);

%% Supp. Fig. 3E 
% f(u) on timeave matrix, aggregate single nuclei first - subtract later - color by time

%load "analysis_mat_timeaverage_042624.mat" from "SuppFig3" folder

E0 = 0.8;


L_reg = [analysis_mat_timeaverage.L_reg_timeave];
L_amp = [analysis_mat_timeaverage.L_amp];
phi = L_reg./L_amp;
timeHere = [analysis_mat_timeaverage.hpa];      



% bin average for single ray-time
averageKTR_u_cell = arrayfun(@(s)...
    adjust_neg(bin_average(s.u,...
    s.ktr,10,0,1)-E0),... %changed from 10 to 20
    analysis_mat_timeaverage,'UniformOutput',false);

mean_ktr_bin10_offsetE0 = cellfun(@nanmean,averageKTR_u_cell);



% normalize, e = Erk(u), m = mean(Erk)
norm_averageKTR_u_cell = arrayfun(@(e,m)(e{:})./(m),averageKTR_u_cell,mean_ktr_bin10_offsetE0,'UniformOutput',false);

norm_averageKTR_u = cell2mat(norm_averageKTR_u_cell);

%
% fit not-timeaveraged data with a line
% bin average
data.x = repmat([0.05:0.1:0.95]',[1,size(norm_averageKTR_u,2)]);
data.y = norm_averageKTR_u;
[binned.y,binned.x,binned.sem] = bin_average(data.x(:),data.y(:),10,0,1);
fit_x = linspace(0,1,100);


% fit
xdata = data.x(~isnan(data.y));
ydata = data.y(~isnan(data.y));
ft = fittype({'-x','1'});
foKTR = fitoptions('Method','LinearLeastSquares',...
    'Lower',[-Inf,0],...
    'Weights', []);
[fit_obj,fit_gof] = fit(xdata(:),ydata(:),ft,foKTR);

% plot
tLevels = 101;
color_map=colormap(cbrewer2('Blues',tLevels));

tt=[];
color_res = 100;
color_code = phi;
color_code = L_amp;

c_digit = -2;
cmax = round(max(color_code),c_digit);
cmin = round(min(color_code),c_digit);
color_unit = (cmax-cmin)./color_res; % config color

f = figure('Visible','on');hold on;

cm = colormap(flipud(cool(round(color_res)+1))); % config color


hold on;
for i = numel(analysis_mat_timeaverage):-1:1
%     s = analysis_mat_timeaverage(i);
%     sz = min([numel(s.averageKTR),numel(s.fractionGEM)]);

    cl_idx = round((color_code(i)-cmin)./color_unit)+1;
    if cl_idx <= 0
        cl_idx = 1;
    elseif cl_idx > round(color_res)+1
        cl_idx = round(color_res)+1;
    end
    cl = cm(cl_idx,:); % config color

    x_plot = data.x(:,i);
    y_plot = data.y(:,i);

    x_plot = x_plot(~isnan(y_plot));
    y_plot = y_plot(~isnan(y_plot));
%     if numel(y_plot)<10
%         tt = [tt,i];
%     end
    colorCode = round(((timeHere(i)-60)./264).*100+1);
    plot(x_plot,y_plot,'-','color',[color_map(colorCode,:),0.8],'markersize',5,'LineWidth',1.5);

end
hold off;
colormap(color_map);
%c = colorbar('Ticks',0:0.5:1,'TickLabels', round((linspace(cmin,cmax,3)),c_digit));
c = colorbar;
set(c,'Ticks',0:0.0455:1,'TickLabels',{'','3','','','','5','','','','7','','','','9','','','','11','','','','13',''});
c.Label.String = 'Time (dpa)';

hold on;
%plot([0.05:0.1:0.95]',mean(data.y,2),'-','color','k','markersize',5,'LineWidth',1.5);

fit_gof.rsquare = round(fit_gof.rsquare,1);

 plot_fit(fit_obj,fit_x,fit_gof,'confidence',false);
ylim([0,2.5])
xlabel('Normalized Space ({\itu})')
ylabel('Normalized ERK')
% title(['fit without time-average, hpa', num2str(t_plot)], ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b-fit_obj.a)])
%title(['fit rolling window timeave'], ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b-fit_obj.a)])

set(gca,'XDir','reverse')
config_plot(f,c);

