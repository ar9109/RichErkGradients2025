%% Supp. Fig. 3B & C

%load "analysis_mat_ERKthresh10.mat" from "SuppFig2" Folder

analysis_mat = analysis_mat_ERKthresh10;

%set up colors
tLevels = 12;
brew = cbrewer2('Reds',tLevels);
%color_map=colormap(cbrewer2('Reds',tLevels));
color_map=colormap(brew(2:12,:));
LampMin = 706;
LampMax = 4053;
%LampColorThresh = linspace(LampMin,LampMax,5);
%LampColorThresh = round(LampColorThresh);
close;

%fishUnique = [1 3 4 5];
fishUnique = [7];
%rayUnique = [2 3 6 7];
rayUnique = [3,7];

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
time = [analysis_mat.hpaTrue];

%ft = fittype('a*(1-x)+b');
% ft = fittype({'1-x','1'});
% foKTR = fitoptions('Method','LinearLeastSquares',...
%     'Lower',[0,-inf],...
%     'Weights', []);
ft = fittype('a*(1-x)+b');
ft = fittype({'1-x','1'});
foKTR = fitoptions('Method','LinearLeastSquares',...
    'Lower',[0,0],...
    'Weights', []);

% f = figure('Visible','on');
% hold on;

% plot individual fish
    for fishNum = 1:size(fishUnique,2)
        fishHere1 = fishUnique(fishNum);
        %fishHere1 = 7;
        for rayNum = 1:size(rayUnique,2)
            f = figure('Visible','on');
            hold on;
            rayHere1 = rayUnique(rayNum);
            %rayHere1 = 3;
            xdataCollect = [];
            ydataCollect = [];
            timeCollect = [];
%             for hpaNum = 1:size(hpaUnique,2)
%                 hpaHere = hpaUnique(hpaNum);
                for aa = 1:size(analysis_mat,2)
                    if analysis_mat(aa).fish == fishHere1 && analysis_mat(aa).ray == rayHere1
                    fishHere = analysis_mat(aa).fish;
                    rayHere = analysis_mat(aa).ray;
                    timeHere = analysis_mat(aa).hpa;
                    LampHere = analysis_mat(aa).L_amp;
                    xdataHere = phi(aa);
                    ydataHere = mean_ktr(aa)-E0;
                    xdataCollect = vertcat(xdataCollect,xdataHere);
                    ydataCollect = vertcat(ydataCollect,ydataHere);
                    timeCollect = vertcat(timeCollect,timeHere);
                    else
                        continue
                    end
                end
            %end

            %if ~isempty(xdataCollect)
            if size(xdataCollect,1)>1

            %xdata = xdataCollect;
            ydata = ydataCollect;
            xdata = timeCollect;
            xdata2 = xdataCollect;

%             fishHere
%             fishHere1
%             rayHere
%             rayHere1
%             timeHere

            [fit_obj,fit_gof] = fit(xdata(:)./24,ydata(:),ft,foKTR);
        
            %phiMax = max(xdata);
            timeMax = max(xdata);
            
            fit_x = linspace(0,timeMax,100);
            
%             f = figure('Visible','on');
%             hold on;
 
%             if rayHere < 4
%                 cl = [0.3010 0.7450 0.9330 0.5];
%             else
%                 cl = [1,.6,1 0.5];
%             end
            
            %%%phiFinal = xdata(size(xdata,1));



             colorCode = round((((LampHere-LampMin)./LampMax)).*10);
% 
             cl = color_map(colorCode+1,:);
% 
%             cl = horzcat(cl,0.5);

            p = figure;
            %left_color = [0 0 0];
            left_color = [0 0 0];
            right_color = [0.7 0.7 0.7];
            set(p,'defaultAxesColorOrder',[left_color; right_color]);
            

            yyaxis left

            if rayHere == 3
                clRight = [0.9256 0.2197 0.1649 0.5];
            elseif rayHere == 7
                clRight = [0.9928 0.8185 0.7419 0.5];
            end
            p = plot(xdata./24,ydata,'.','Color',clRight,'MarkerSize',30,'linewidth',2); hold on;
            pfit = plot(fit_obj);
            set(pfit,'linewidt',2,'color',clRight);
            
            ylim([0,0.6]);
            xlim([0,400]);
            ylabel('Average ERK Activity (rescaled space,u)')
            yticks([0 0.1 0.2 0.3 0.4 0.5]);

            yyaxis right
            plot(xdata./24,xdata2,'.--','color',[0.7 0.7 0.7],'MarkerSize',30,'linewidth',2);
            
% %             p_fit = plot_fit(fit_obj,fit_x,fit_gof,confidence=false); hold on;
% %             set(p_fit,'linestyle','-','color','k','lineWidth',3)
        
            else 
                continue
            end
            xdataCollect = [];
            ydataCollect = [];
% % %             close all;


            xlabel('Time (days post amp.)')
            ylabel('Fraction Regenerated')
            title(['fish' num2str(fishHere) ' ' 'ray' num2str(rayHere)]);
            
            config_plot(f);
            
            ylim([0,1]);
            xlim([0,18]);
            yticks([0 0.2 0.4 0.6 0.8 1])
            legend('off');
            set(gca,'fontsize',16);
            
            %close figure 1;
            %close figure 3;



        end
    end 

%% Supp. Fig. 3D

%load "analysis_mat_WTdataSet.mat" from "SuppFig2" folder

% Plot phi v. time - red - all fish

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
                    phiHere = LregHere./LampHere;
                    averageKTRhere = analysis_mat(i).averageKTR;
                    
                    BaseName='Fish';
                    BaseName2='Ray';
                    BaseName3='hpa';
                    FileName1=[BaseName,num2str(fishHere)];
                    FileName2=[BaseName2,num2str(rayHere)];
                    FileName3=[num2str(hppHere),BaseName3];

                avgKTRray2 = nanmean(averageKTRhere)-0.8;

                plotXcollect = vertcat(plotXcollect,hppHere);
                plotYcollect = vertcat(plotYcollect,phiHere);

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
    ylabel('Fraction Regenerated');
    %title(FileName1);
    set(gca, 'fontsize', 20);
    
    %c= colorbar('Ticks', 48:24:336,'TickLabels',{'48','72','96','120','144','168','192','216','240','264','288','312','336'});
    %caxis([tmin tmax])
    %c.FontSize = 10;
    %c.Label.String = 'Time (hpa)';
    set(gca, 'fontsize', 16);
    end
end

%% Supp. Fig. 3E
% Average ERK by time color by Lamp

%load "analysis_mat_WTdataSet.mat" from "SuppFig2" folder

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
time = [analysis_mat.hpaTrue];
ray = [analysis_mat.ray];
% End added 20 apr 23

ft = fittype('a*(1-x)+b');
ft = fittype({'1-x','1'});
foKTR = fitoptions('Method','LinearLeastSquares',...
    'Lower',[0,0],...
    'Weights', []);
xdata = time;
ydata = mean_ktr-E0;

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

fit_x = linspace(round(min(time))-1,round(max(time))+1,100);

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
xlabel('Time (hours post amp.)')
ylabel('Average ERK Activity (rescaled space,u)')
% title(['fit all data, E0 =', num2str(E0)],...
%     ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b)])
% title(['fit all data'],...
%     ['y = ',num2str(fit_obj.a),'(1-x) + ',num2str(fit_obj.b)])

config_plot(f);

legend_name={strcat("R^2 = ",num2str(round(fit_gof.rsquare,1))),...
    strcat("R^2 = ", num2str(round(fit_gof_small.rsquare,1))),...
    strcat("R^2 = ", num2str(round(fit_gof_med.rsquare,1))),...
    strcat("R^2 = ", num2str(round(fit_gof_large.rsquare,1)))};

legend([p_fit,p_fit_small,p_fit_med,p_fit_large], legend_name,...
    'Location','best','color','none','box','off');

color_map=colormap(cbrewer2('Reds',tLevels));
g = colorbar;
%set(g,'TickLabels',{'0','811','1621','2432','3242','4053'});%,'2740','3080','3420','3760','4100'})
set(g,'TickLabels',{'700','1380','2060','2740','3420','4100'});
g.Label.String = 'Length Amputated (um)';


ylim([0,0.5]);
%%%ylim([0.6,1.4])
%xlim([0,1.4])