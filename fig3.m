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