%% Fig. 4B - left
% Plot binned average ERK v. position - Used 1 May 24 - color Blue

%load "analysis_mat_WTdataSet.mat" from "Fig4" Folder


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

%% Fig. 4B - right
% Plot binned average ERK v. position - Used 1 May 24 - color Blue

%load "analysis_mat_WTdataSet.mat" from "Fig4" Folder


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
%% Fig. 4C & 4D

%load "analysis_mat_timeaverage_042624.mat" from "Fig4" folder

%step 1 - prep data
% aggregate single nuclei first - subtract later - bin data into 5 groups

% remove rays earlier than 72 hpa
% field_name = 'fit_ktr_linear';
% not_empty = arrayfun(@(x)~isempty(x.(field_name)),analysis_mat_timeaverage);

% E0 = 0.8483;
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

%separateData by Lamp
bin_x_1amp =  [];
bin_y_1amp = [];
bin_x_2amp = [];
bin_y_2amp = [];
bin_x_3amp = [];
bin_y_3amp = [];
bin_x_4amp = [];
bin_y_4amp = [];
bin_x_5amp = [];
bin_y_5amp = [];
for i = 1:size(data.x,2)
    if L_amp(i) < 1301
        bin_x_1amp= horzcat(bin_x_1amp,data.x(:,i));
        bin_y_1amp= horzcat(bin_y_1amp,data.y(:,i));
    elseif L_amp(i) < 2001
        bin_x_2amp = horzcat(bin_x_2amp,data.x(:,i));
        bin_y_2amp = horzcat(bin_y_2amp,data.y(:,i));
    elseif L_amp(i) < 2701
        bin_x_3amp = horzcat(bin_x_3amp,data.x(:,i));
        bin_y_3amp = horzcat(bin_y_3amp,data.y(:,i));
        elseif L_amp(i) < 3401
        bin_x_4amp = horzcat(bin_x_4amp,data.x(:,i));
        bin_y_4amp = horzcat(bin_y_4amp,data.y(:,i));
    else
        bin_x_5amp = horzcat(bin_x_5amp,data.x(:,i));
        bin_y_5amp = horzcat(bin_y_5amp,data.y(:,i));
    end
end

%separateData by time
bin_x_1time =  [];
bin_y_1time = [];
bin_x_2time = [];
bin_y_2time = [];
bin_x_3time = [];
bin_y_3time = [];
bin_x_4time = [];
bin_y_4time = [];
bin_x_5time = [];
bin_y_5time = [];
for i = 1:size(data.x,2)
    if timeHere(i) < 111
        bin_x_1time= horzcat(bin_x_1time,data.x(:,i));
        bin_y_1time= horzcat(bin_y_1time,data.y(:,i));
    elseif timeHere(i) < 171
        bin_x_2time = horzcat(bin_x_2time,data.x(:,i));
        bin_y_2time = horzcat(bin_y_2time,data.y(:,i));
    elseif timeHere(i) < 231
        bin_x_3time = horzcat(bin_x_3time,data.x(:,i));
        bin_y_3time = horzcat(bin_y_3time,data.y(:,i));
    elseif timeHere(i) < 291
        bin_x_4time = horzcat(bin_x_4time,data.x(:,i));
        bin_y_4time = horzcat(bin_y_4time,data.y(:,i));
    else
        bin_x_5time = horzcat(bin_x_5time,data.x(:,i));
        bin_y_5time = horzcat(bin_y_5time,data.y(:,i));
    end
end


%step 2
% color map
tLevels = 10;
color_map_blue=colormap(cbrewer2('Blues',tLevels));

color_map_red=colormap(cbrewer2('Reds',tLevels));



%step 3 - plot Fig. 4C
% plot f(u) average of groups


g = figure;

b = errorbar(mean(bin_x_1amp,2),nanmean(bin_y_1amp,2),(nanstd(bin_y_1amp,0,2)./sqrt(size(bin_y_1amp,2))),'-','color',[0.9966 0.8899 0.8395],'markersize',10,'LineWidth',4); hold on;
c = errorbar(mean(bin_x_2amp,2),nanmean(bin_y_2amp,2),(nanstd(bin_y_2amp,0,2)./sqrt(size(bin_y_2amp,2))),'-','color',[0.9884 0.6262 0.5059],'markersize',10,'LineWidth',4);
d = errorbar(mean(bin_x_3amp,2),nanmean(bin_y_3amp,2),(nanstd(bin_y_3amp,0,2)./sqrt(size(bin_y_3amp,2))),'-','color',[0.9746 0.3327 0.2330],'markersize',10,'LineWidth',4);
e = errorbar(mean(bin_x_4amp,2),nanmean(bin_y_4amp,2),(nanstd(bin_y_4amp,0,2)./sqrt(size(bin_y_4amp,2))),'-','color',[0.7664 0.0711 0.1052],'markersize',10,'LineWidth',4);
f = errorbar(mean(bin_x_5amp,2),nanmean(bin_y_5amp,2),(nanstd(bin_y_5amp,0,2)./sqrt(size(bin_y_5amp,2))),'-','color',[0.4039 0 0.0510],'markersize',10,'LineWidth',4);
a = errorbar([0.05:0.1:0.95]',nanmean(data.y,2),nanstd(data.y,0,2)./sqrt(size(data.y,2)),'-','color','k','markersize',10,'LineWidth',4); 
hold off;

ylim([0,2.5])
xlabel('Normalized Position ({\itu})')
ylabel('ERK Activity / Average ERK Activity')

legend_name={strcat("All Data (180)")...
    strcat("601-1300 \mum (37)"),...
        strcat("1301-2000 \mum (50)"),...
        strcat("2001-2700 \mum (14)"),...
        strcat("2701-3400 \mum (42)"),...
        strcat("3401-4100 \mum (37)")};
legend([a,b,c,d,e,f], legend_name,...
        'Location','northwest','color','none','box','off');


set(gca,'XDir','reverse','fontsize',16)

%step 4 - plot 4D
% plot f(u) average of groups - by Time
g2 = figure;

bb = errorbar(mean(bin_x_1time,2),nanmean(bin_y_1time,2),(nanstd(bin_y_1time,0,2)./sqrt(size(bin_y_1time,2))),'-','color',[0.8812    0.9286    0.9721],'markersize',5,'LineWidth',4); hold on;
cc = errorbar(mean(bin_x_2time,2),nanmean(bin_y_2time,2),(nanstd(bin_y_2time,0,2)./sqrt(size(bin_y_2time,2))),'-','color',[0.6744    0.8178    0.9018],'markersize',5,'LineWidth',4);
dd = errorbar(mean(bin_x_3time,2),nanmean(bin_y_3time,2),(nanstd(bin_y_3time,0,2)./sqrt(size(bin_y_3time,2))),'-','color',[0.3454    0.6344    0.8125],'markersize',5,'LineWidth',4);
ee = errorbar(mean(bin_x_4time,2),nanmean(bin_y_4time,2),(nanstd(bin_y_4time,0,2)./sqrt(size(bin_y_4time,2))),'-','color',[0.1116    0.4148    0.6906],'markersize',5,'LineWidth',4);
ff = errorbar(mean(bin_x_5time,2),nanmean(bin_y_5time,2),(nanstd(bin_y_5time,0,2)./sqrt(size(bin_y_5time,2))),'-','color',[0.0314    0.1882    0.4196],'markersize',5,'LineWidth',4);
aa = errorbar([0.05:0.1:0.95]',nanmean(data.y,2),(nanstd(data.y,0,2)./sqrt(size(data.y,2))),'-','color','k','markersize',5,'LineWidth',4); 
hold off

ylim([0,2.5])
xlabel('Normalized Position ({\itu})')
ylabel('ERK Activity / Average ERK Activity')
% title(['fit without time-average, hpa', num2str(t_plot)], ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b-fit_obj.a)])
%title(['fit rolling window timeave'], ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b-fit_obj.a)])
set(gca,'XDir','reverse','fontsize',16)

% legend_name={strcat("All Data (180)")...
%     strcat("0-3.5 dpa (42)"),...
%         strcat("4-5.5 dpa (28)"),...
%         strcat("6-7 dpa (36)"),...
%         strcat("7.5-10 dpa (42)"),...
%         strcat("10.5-14 dpa (32)")};
legend_name={strcat("All Data (180)")...
    strcat("2.5-4.5 dpa (58)"),...
        strcat("5-7 dpa (48)"),...
        strcat("7.5-9.5 dpa (36)"),...
        strcat("10-12 dpa (26)"),...
        strcat("12.5+ dpa (12)")};
legend([aa,bb,cc,dd,ee,ff], legend_name,...
        'Location','northwest','color','none','box','off');

%%  Fig. 4E

%load "analysis_mat_timeaverage_042624.mat" from "Fig4" folder

%step 1
% aggregate single nuclei first - subtract later - bin data into 5 groups

% remove rays earlier than 72 hpa
% field_name = 'fit_ktr_linear';
% not_empty = arrayfun(@(x)~isempty(x.(field_name)),analysis_mat_timeaverage);

% E0 = 0.8483;
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

averageKTR_u = cell2mat(averageKTR_u_cell);

%
% fit not-timeaveraged data with a line
% bin average
data.x = repmat([0.05:0.1:0.95]',[1,size(averageKTR_u,2)]);
data.y = averageKTR_u;
[binned.y,binned.x,binned.sem] = bin_average(data.x(:),data.y(:),10,0,1);
fit_x = linspace(0,1,100);

%separateData by time

bin_y_1phi = [];
bin_y_2phi = [];
bin_y_3phi = [];
bin_y_4phi = [];
bin_y_5phi = [];
bin_y_6phi = [];
bin_y_7phi = [];
bin_y_8phi = [];
bin_y_9phi = [];
bin_y_10phi = [];
for i = 1:size(data.x,2)
    if phi(i) < 0.1
        bin_y_1phi= horzcat(bin_y_1phi,data.y(:,i));
    elseif phi(i) < 0.2
        bin_y_2phi = horzcat(bin_y_2phi,data.y(:,i));
    elseif phi(i) < 0.3
        bin_y_3phi = horzcat(bin_y_3phi,data.y(:,i));
    elseif phi(i) < 0.4
        bin_y_4phi = horzcat(bin_y_4phi,data.y(:,i));
    elseif phi(i) < 0.5
        bin_y_5phi = horzcat(bin_y_5phi,data.y(:,i));
    elseif phi(i) < 0.6
        bin_y_6phi = horzcat(bin_y_6phi,data.y(:,i));
    elseif phi(i) < 0.7
        bin_y_7phi = horzcat(bin_y_7phi,data.y(:,i));
    elseif phi(i) < 0.8
        bin_y_8phi = horzcat(bin_y_8phi,data.y(:,i));
    elseif phi(i) < 0.9
        bin_y_9phi = horzcat(bin_y_9phi,data.y(:,i));
    else
        bin_y_10phi = horzcat(bin_y_10phi,data.y(:,i));
    end
end

%step2

%separateData by phi
phiCol1 =  [];
phiCol2 =  [];
phiCol3 =  [];
phiCol4 =  [];
phiCol5 =  [];
phiCol6 =  [];
phiCol7 =  [];
phiCol8 =  [];
phiCol9 =  [];
phiCol10 =  [];

for i = 1:size(data.x,2)
    if phi(i) < 0.1
        phiCol1= horzcat(phiCol1,phi(i));
    elseif phi(i) < 0.2
        phiCol2= horzcat(phiCol2,phi(i));
    elseif phi(i) < 0.3
        phiCol3= horzcat(phiCol3,phi(i));
    elseif phi(i) < 0.4
        phiCol4= horzcat(phiCol4,phi(i));
    elseif phi(i) < 0.5
        phiCol5= horzcat(phiCol5,phi(i));
    elseif phi(i) < 0.6
        phiCol6= horzcat(phiCol6,phi(i));
    elseif phi(i) < 0.7
        phiCol7= horzcat(phiCol7,phi(i));
    elseif phi(i) < 0.8
        phiCol8= horzcat(phiCol8,phi(i));
    elseif phi(i) < 0.9
        phiCol9= horzcat(phiCol9,phi(i));
    else
        phiCol10= horzcat(phiCol10,phi(i));
    end
end

meanPhiCol1 = mean(phiCol1);
meanPhiCol2 = mean(phiCol2);
meanPhiCol3 = mean(phiCol3);
meanPhiCol4 = mean(phiCol4);
meanPhiCol5 = mean(phiCol5);
meanPhiCol6 = mean(phiCol6);
meanPhiCol7 = mean(phiCol7);
meanPhiCol8 = mean(phiCol8);
meanPhiCol9 = mean(phiCol9);
meanPhiCol10 = mean(phiCol10);

medPhiCol1 = median(phiCol1);
medPhiCol2 = median(phiCol2);
medPhiCol3 = median(phiCol3);
medPhiCol4 = median(phiCol4);
medPhiCol5 = median(phiCol5);
medPhiCol6 = median(phiCol6);
medPhiCol7 = median(phiCol7);
medPhiCol8 = median(phiCol8);
medPhiCol9 = median(phiCol9);
medPhiCol10 = median(phiCol10);

%step3
% generate ERK v. x/Lamp - purple
g3 = figure;

tLevels = 20;
color_map_red=colormap(cbrewer2('Purples',tLevels));

bbb = errorbar(medPhiCol1-linspace(0,medPhiCol1,10),nanmean(bin_y_1phi,2),(nanstd(bin_y_1phi,0,2)./sqrt(size(bin_y_1phi,2))),'-','color',color_map_red(2,:),'markersize',5,'LineWidth',4); hold on;
ccc = errorbar(medPhiCol2-linspace(0,medPhiCol2,10),nanmean(bin_y_2phi,2),(nanstd(bin_y_2phi,0,2)./sqrt(size(bin_y_2phi,2))),'-','color',color_map_red(4,:),'markersize',5,'LineWidth',4);
ddd = errorbar(medPhiCol3-linspace(0,medPhiCol3,10),nanmean(bin_y_3phi,2),(nanstd(bin_y_3phi,0,2)./sqrt(size(bin_y_3phi,2))),'-','color',color_map_red(6,:),'markersize',5,'LineWidth',4);
eee = errorbar(medPhiCol4-linspace(0,medPhiCol4,10),nanmean(bin_y_4phi,2),(nanstd(bin_y_4phi,0,2)./sqrt(size(bin_y_4phi,2))),'-','color',color_map_red(8,:),'markersize',5,'LineWidth',4);
fff = errorbar(medPhiCol5-linspace(0,medPhiCol5,10),nanmean(bin_y_5phi,2),(nanstd(bin_y_5phi,0,2)./sqrt(size(bin_y_5phi,2))),'-','color',color_map_red(10,:),'markersize',5,'LineWidth',4);
ggg = errorbar(medPhiCol6-linspace(0,medPhiCol6,10),nanmean(bin_y_6phi,2),(nanstd(bin_y_6phi,0,2)./sqrt(size(bin_y_6phi,2))),'-','color',color_map_red(12,:),'markersize',5,'LineWidth',4);
hhh = errorbar(medPhiCol7-linspace(0,medPhiCol7,10),nanmean(bin_y_7phi,2),(nanstd(bin_y_7phi,0,2)./sqrt(size(bin_y_7phi,2))),'-','color',color_map_red(14,:),'markersize',5,'LineWidth',4);
iii = errorbar(medPhiCol8-linspace(0,medPhiCol8,10),nanmean(bin_y_8phi,2),(nanstd(bin_y_8phi,0,2)./sqrt(size(bin_y_8phi,2))),'-','color',color_map_red(16,:),'markersize',5,'LineWidth',4);
jjj = errorbar(medPhiCol9-linspace(0,medPhiCol9,10),nanmean(bin_y_9phi,2),(nanstd(bin_y_9phi,0,2)./sqrt(size(bin_y_9phi,2))),'-','color',color_map_red(18,:),'markersize',5,'LineWidth',4);
kkk = errorbar(medPhiCol10-linspace(0,medPhiCol10,10),nanmean(bin_y_10phi,2),(nanstd(bin_y_10phi,0,2)./sqrt(size(bin_y_10phi,2))),'-','color',color_map_red(20,:),'markersize',5,'LineWidth',4);

hold off
 
ylim([0,0.5])
xlabel('Average x/Length Amputated')
ylabel('ERK Activity')
yticks([0 0.1 0.2 0.3 0.4 0.5])
% title(['fit without time-average, hpa', num2str(t_plot)], ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b-fit_obj.a)])
%title(['fit rolling window timeave'], ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b-fit_obj.a)])
%set(gca,'XDir','reverse','fontsize',16)
set(gca,'fontsize',16)
xlim([0,1]);


% legend_name={strcat("All Data (180)")...
%     strcat("0-3.5 dpa (42)"),...
%         strcat("4-5.5 dpa (28)"),...
%         strcat("6-7 dpa (36)"),...
%         strcat("7.5-10 dpa (42)"),...
%         strcat("10.5-14 dpa (32)")};
legend_name={
    strcat("phi < 0.1 (5)"),...
        strcat("phi < 0.2 (25)"),...
        strcat("phi < 0.3 (15)"),...
        strcat("phi < 0.4 (4)"),...
        strcat("phi < 0.5 (22)"),...
        strcat("phi < 0.6 (38)"),...
        strcat("phi < 0.7 (22)"),...
        strcat("phi < 0.8 (24)"),...
        strcat("phi < 0.9 (15)"),...
        strcat("phi > 0.9 (10)")};
legend([bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk], legend_name,...
        'Location','northeast','color','none','box','off');



%% Fig. 4F