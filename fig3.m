%% Fig. 3b

%load "analysis_mat_WTdataSet.mat" from "Fig3" folder

% Plot lreg v. time - red - all fish
 
figure;

ERKall = [];
%for fish = 2
for fish = 1:55
    fishNow = fish;
    
    ERKcollect = [];
    xValuesCollect = [];
    xValuesMean = [];
    ERKcollectMean = [];

    for ray = 2:10
    %for ray = 3:7
        rayNow = ray;
        
        %for hpp = 72:336
            %hppHere = hpp;
            hppMin = 72;
            %figure;

            plotXcollect = [];
            plotYcollect = [];
            
            for i = 1:size(analysis_mat,2)
                if analysis_mat(i).fish == fishNow && analysis_mat(i).ray == rayNow &&  analysis_mat(i).hpa >= hppMin

                    fishHere = analysis_mat(i).fish;
                    rayHere = analysis_mat(i).ray;
                    hppHere = analysis_mat(i).hpa;

                    LregHere = analysis_mat(i).L_reg;
                    LampHere = analysis_mat(i).L_amp;
                    averageKTRhere = analysis_mat(i).averageKTR;
                    
                    BaseName='Fish';
                    BaseName2='Ray';
                    BaseName3='hpa';
                    FileName1=[BaseName,num2str(fishHere)];
                    FileName2=[BaseName2,num2str(rayHere)];
                    FileName3=[num2str(hppHere),BaseName3];

                avgKTRray2 = nanmean(averageKTRhere)-0.8;
            
                %Set up colors
                tLevels = 288;
                colormap(jet(tLevels));
                colors = jet(tLevels);
                tmin = 48;
                tmax = 336;
            
                %idxColor = round((tLevels-1)*(hppHere-tmin)/(tmax-tmin)+1);
                if rayHere < 4 
                    cl = [0.3010 0.7450 0.9330];
                else
                    cl = [1 .6 1];
                end

                plotXcollect = vertcat(plotXcollect,hppHere);
                plotYcollect = vertcat(plotYcollect,LregHere);

                end
            end

    if ~isempty(plotXcollect)
    %set up color map
    tLevels = 101;
    color_map=colormap(cbrewer2('Reds',tLevels));
    lampMinHere = 706;
    lampMaxHere = 4053;            
    colorCode = ((LampHere-706)./3353).*100+1; %%% 699 will need to be edited 
    
    plot(plotXcollect./24,plotYcollect,'.-','color',color_map((round(colorCode)),:),'MarkerSize',30,'linewidth',2); hold on;
    else
        continue
    end
    
    %ylim([0 1500]);
    %xlim([4 8]);
    xlabel('Time (dpa)');
    ylabel('Length Regenerated (\mum)');
    title(FileName1);
    set(gca, 'fontsize', 20);
    
    %c= colorbar('Ticks', 48:24:336,'TickLabels',{'48','72','96','120','144','168','192','216','240','264','288','312','336'});
    %caxis([tmin tmax])
    %c.FontSize = 10;
    %c.Label.String = 'Time (hpa)';
    set(gca, 'fontsize', 20);
    end
end

%% Fig. 3c

%load "analysis_mat_WTdataSet.mat" from "Fig3" folder

% A(phi) color by Lamp - 3 groups

E0 = 0.8;

mean_ktr = arrayfun(@(x)mean(x.ktr(x.trim_logical)),analysis_mat);
mean_ktr_bin = arrayfun(@(x)nanmean(x.averageKTR),analysis_mat);
averageKTR_u_cell = arrayfun(@(s)...
    bin_average(s.ccrot(s.trim_logical,1),...
    adjust_neg(s.ktr(s.trim_logical)),10,0,s.L_reg),...
    analysis_mat,'UniformOutput',false);

mean_ktr_bin10_offsetE0 = cellfun(@nanmean,averageKTR_u_cell);

L_reg = [analysis_mat.L_reg];
L_amp = [analysis_mat.L_amp];
phi = L_reg./[L_amp];
ray = [analysis_mat.ray];

ft = fittype('a*(1-x)+b');
ft = fittype({'1-x','1'});
foKTR = fitoptions('Method','LinearLeastSquares',...
    'Lower',[0,0],...
    'Weights', []);
xdata = phi;
ydata = mean_ktr-E0;
fit_x = linspace(0,1,100);


%Generate Lamp
LampCollect = [];
for qq = 1:size(analysis_mat(1:end),2)
    LampLocal = analysis_mat(qq).L_amp;
    LampCollect = vertcat(LampCollect,LampLocal);
end

LampCollect = LampCollect';


[fit_obj,fit_gof] = fit(xdata(:),ydata(:),ft,foKTR);

tLevels = 101;
color_map=colormap(cbrewer2('Reds',tLevels));
lampMinHere = 706;
lampMaxHere = 4053;

%cm = lines(10);
f = figure('Visible','on');
hold on;
xdata_long = xdata(ray<4);
ydata_long = ydata(ray<4);
xdata_short = xdata(ray>=4);
ydata_short = ydata(ray>=4);

bin_x_small =  [];
bin_y_small = [];
bin_x_med = [];
bin_y_med = [];
bin_x_large = [];
bin_y_large = [];

for zz = 1:size(xdata,2)
    colorCode = ((LampCollect(zz)-700)./3353).*100+1;
    p_all = plot(xdata(zz),ydata(zz),'o','color',color_map((round(colorCode)),:),'linewidth',2); hold on;
    if LampCollect(zz) < 1575
        bin_x_small = vertcat(bin_x_small,xdata(zz));
        bin_y_small = vertcat(bin_y_small,ydata(zz));
    elseif LampCollect(zz) < 3089
        bin_x_med = vertcat(bin_x_med,xdata(zz));
        bin_y_med = vertcat(bin_y_med,ydata(zz));
    else
        bin_x_large = vertcat(bin_x_large,xdata(zz));
        bin_y_large = vertcat(bin_y_large,ydata(zz));
    end
end
% p_long = plot(xdata_long,ydata_long,'o','Color',cm(6,:),'LineWidth',2);
% p_short = plot(xdata_short,ydata_short,'o','Color',[1,.6,1],'LineWidth',2);

[fit_obj_small,fit_gof_small] = fit(bin_x_small,bin_y_small,ft,foKTR);
[fit_obj_med,fit_gof_med] = fit(bin_x_med,bin_y_med,ft,foKTR);
[fit_obj_large,fit_gof_large] = fit(bin_x_large,bin_y_large,ft,foKTR);

p_fit_small = plot_fit(fit_obj_small,fit_x,fit_gof_small,confidence=false);
set(p_fit_small,'linewidth',2,'color',[1 0.5 0.5]);
p_fit_med = plot_fit(fit_obj_med,fit_x,fit_gof_med,confidence=false);
set(p_fit_med,'linewidth',2,'color',[1 0 0]);
p_fit_large = plot_fit(fit_obj_large,fit_x,fit_gof_large,confidence=false);
set(p_fit_large,'linewidth',2,'color',[0.5 0 0]);
p_fit = plot_fit(fit_obj,fit_x,fit_gof,confidence=false);


xlabel('Fraction Regenerated')
ylabel('Average ERK Activity (rescaled space,u)')
% title(['fit all data, E0 =', num2str(E0)],...
%     ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b)])
% title(['fit all data'],...
%     ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b)])

config_plot(f);

%if ~isempty(ft)
    %legend([fit_plot], {['R^2 = ' num2str(gof2.rsquare)]},'Location','best','color','none','box','off');
% % %     legend_name={strcat("R^2 = ",num2str(gof.rsquare)),...
% % %         strcat("R^2 = ", num2str(gof_mag.rsquare))};

    legend_name={strcat("R^2 = ",num2str(round(fit_gof.rsquare,1))),...
        strcat("R^2 = ", num2str(round(fit_gof_small.rsquare,1))),...
        strcat("R^2 = ", num2str(round(fit_gof_med.rsquare,1))),...
        strcat("R^2 = ", num2str(round(fit_gof_large.rsquare,1)))};
% % %     legend([fit_plot,fit_plot_mag], legend_name,...
% % %         'Location','best','color','none','box','off');

    legend([p_fit,p_fit_small,p_fit_med,p_fit_large], legend_name,...
        'Location','best','color','none','box','off');
    %legend([fit_plot,fit_plot_mag],["1","2"])
%end

color_map=colormap(cbrewer2('Reds',tLevels));
g = colorbar;
%set(g,'TickLabels',{'0','811','1621','2432','3242','4053'});%,'2740','3080','3420','3760','4100'})
set(g,'TickLabels',{'700','1380','2060','2740','3420','4100'});
g.Label.String = 'Length Amputated (um)';

%% Fig. 3D & 3E

%load "analysis_mat_timeaverage_042624.mat" from "Fig3" folder

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



%step 3 - plot Fig. 3D
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

%step 4 - plot 3E
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

%%  Fig. 3F

%load "analysis_mat_timeaverage_042624.mat" from "Fig3" folder

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
