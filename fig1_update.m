%% Fig. 1e
% Expansion
close all;

c = get(gca,'Children');
xdata = get(c, 'XData');
ydata = get(c, 'YData');

figure;

xdata = [0.9663 0.9574 0.9670 0.9581 0.9531 0.9105 0.8966 0.9182 0.8997 0.8737 0.6944 NaN ...
    0.7648 0.7293 NaN 0.4803 NaN 0.6086 0.5558 NaN, 0.2639 NaN 0.3296 0.3438 NaN 0.0469 ...
    0.1174 0.0438 0.1319 0.1027 0.0132 0.0481 0.0093 0.0561 0.0344];

ydata = [2.6017 1.4573 2.7757 1.2011 0.9905 1.8052 2.3000 2.8831 2.6522 0.9161 0.1046 NaN ...
    0.0245 0.1283 NaN 0.1857 NaN 0.1403 -0.0244 NaN -0.0087 NaN -0.0060 -0.0125 NaN 0 0.2195 ...
    0.0103 -0.1098 0.1132 -0.0435 -0.1351 0 0.0258 0.0280];

plot(xdata,ydata,'.','markersize',30,'color',[0.5 0.5 0.5]);

box on;

ylabel('Expansion Rate (1/day)')
xlabel('Normalized Position')
%xlim([-1,2]);
%ylim([-60,120]);
set(gca,'FontSize',18);

%% Fig. 1f
% Change in N

%a = get(gca,'Children');
%xdata = get(a, 'XData');
%ydata = get(a, 'YData');

distData = [0.6522 0.8378 0.6000 0.7200 0.5000]; 
distX = [3 3 3 3 3];
midData = [0 NaN 0 0.0294 NaN]; 
midX = [2 2 2 2 2];
proxData = [0 -0.0417 -0.0213 0 -0.0213];
proxX = [1 1 1 1 1];

yVals = horzcat(distData,midData,proxData);
xVals = horzcat(distX,midX,proxX);

figure;

b = boxchart(xVals,yVals);
b.BoxFaceColor = [0 0 0];
b.WhiskerLineColor = [0 0 0];
b.BoxFaceAlpha = 0;
%b.MarkerColor = 'k';
b.MarkerColor = [0 0 0];
b.LineWidth = 2;

hold on;
plot([-0.2,-0.1,0,0.1,0.2]+1,proxData,'.','MarkerSize',30,'color',[0.6350 0.0780 0.1840]);
plot([-0.2,-0.1,0,0.1,0.2]+2,midData,'.','MarkerSize',30,'color',[0.4660 0.6740 0.1880]);
plot([-0.2,-0.1,0,0.1,0.2]+3,distData,'.','MarkerSize',30,'color',[0 0.4470 0.7410]);

box on;

ylabel('\DeltaN/N_{i} (Converted Cells)')
%xlim([-1,2]);
%ylim([-60,120]);
set(gca,'FontSize',18);

ax = gca;
ax.XTick = [1, 2, 3];
ax.XTickLabels = {'Proximal','Middle','Distal'};

%% Fig. 1G
% delta N v. delta A

h = get(gca,'Children');
xdata = get(h, 'XData');
ydata = get(h, 'YData');

hEr = findobj(gca,'Type','ErrorBar');
eyNeg=get(hEr,'YNegativeDelta');
eyPos=get(hEr,'YPositiveDelta');

%xLine = cell2mat(xdata(1));
xLine = [0	0.0100000000000000	0.0200000000000000	0.0300000000000000	0.0400000000000000	0.0500000000000000	0.0600000000000000	0.0700000000000000	0.0800000000000000	0.0900000000000000	0.100000000000000	0.110000000000000	0.120000000000000	0.130000000000000	0.140000000000000	0.150000000000000	0.160000000000000	0.170000000000000	0.180000000000000	0.190000000000000	0.200000000000000	0.210000000000000	0.220000000000000	0.230000000000000	0.240000000000000	0.250000000000000	0.260000000000000	0.270000000000000	0.280000000000000	0.290000000000000	0.300000000000000	0.310000000000000	0.320000000000000	0.330000000000000	0.340000000000000	0.350000000000000	0.360000000000000	0.370000000000000	0.380000000000000	0.390000000000000	0.400000000000000	0.410000000000000	0.420000000000000	0.430000000000000	0.440000000000000	0.450000000000000	0.460000000000000	0.470000000000000	0.480000000000000	0.490000000000000	0.500000000000000	0.510000000000000	0.520000000000000	0.530000000000000	0.540000000000000	0.550000000000000	0.560000000000000	0.570000000000000	0.580000000000000	0.590000000000000	0.600000000000000	0.610000000000000	0.620000000000000	0.630000000000000	0.640000000000000	0.650000000000000	0.660000000000000	0.670000000000000	0.680000000000000	0.690000000000000	0.700000000000000	0.710000000000000	0.720000000000000	0.730000000000000	0.740000000000000	0.750000000000000	0.760000000000000	0.770000000000000	0.780000000000000	0.790000000000000	0.800000000000000	0.810000000000000	0.820000000000000	0.830000000000000	0.840000000000000	0.850000000000000	0.860000000000000	0.870000000000000	0.880000000000000	0.890000000000000	0.900000000000000];
%yLine = cell2mat(ydata(1));
yLine = [0	0.0145600000000000	0.0291200000000000	0.0436800000000000	0.0582400000000000	0.0728000000000000	0.0873600000000000	0.101920000000000	0.116480000000000	0.131040000000000	0.145600000000000	0.160160000000000	0.174720000000000	0.189280000000000	0.203840000000000	0.218400000000000	0.232960000000000	0.247520000000000	0.262080000000000	0.276640000000000	0.291200000000000	0.305760000000000	0.320320000000000	0.334880000000000	0.349440000000000	0.364000000000000	0.378560000000000	0.393120000000000	0.407680000000000	0.422240000000000	0.436800000000000	0.451360000000000	0.465920000000000	0.480480000000000	0.495040000000000	0.509600000000000	0.524160000000000	0.538720000000000	0.553280000000000	0.567840000000000	0.582400000000000	0.596960000000000	0.611520000000000	0.626080000000000	0.640640000000000	0.655200000000000	0.669760000000000	0.684320000000000	0.698880000000000	0.713440000000000	0.728000000000000	0.742560000000000	0.757120000000000	0.771680000000000	0.786240000000000	0.800800000000000	0.815360000000000	0.829920000000000	0.844480000000000	0.859040000000000	0.873600000000000	0.888160000000000	0.902720000000000	0.917280000000000	0.931840000000000	0.946400000000000	0.960960000000000	0.975520000000000	0.990080000000000	1.00464000000000	1.01920000000000	1.03376000000000	1.04832000000000	1.06288000000000	1.07744000000000	1.09200000000000	1.10656000000000	1.12112000000000	1.13568000000000	1.15024000000000	1.16480000000000	1.17936000000000	1.19392000000000	1.20848000000000	1.22304000000000	1.23760000000000	1.25216000000000	1.26672000000000	1.28128000000000	1.29584000000000	1.31040000000000];

dFluor = [-0.0248	% not used
0.0227	
0.1375	
0.1231	
0.0513	
0.1569	
0.0751	
0.1087	
0.0914	
0.1251	
0.7019	
1.3286	
0.6543	
1.5063	
0.5955]	;



erFluor = [0.1492 % not used
0.0751
0.1087
0.0914
0.1251
0.0468
NaN
0.216
0.2349
NaN
0.3154
0.3169
0.2421
0.5291
0.247];

dN =[0 -0.0417 -0.0213 0 -0.0213 ...
    0 NaN 0 0.0294 NaN ...
    0.6522 0.8378 0.6000 0.7200 0.5000]';

dArea = [0.1102
0.1204
0.0071
-0.0229
0.0335
0.2159
NaN
0.1861
0.2572
NaN
0.6371
1.2701
0.9075
1.4039
0.7718];

errArea = [0.0256
0.0194
0.0161
0.0187
0.0159
0.02345
NaN
0.024
0.0286
NaN
0.0396
0.0413
0.0382
0.0548
0.0285];

xIdent = [0 0.5 1];
yIdent = [0 0.5 1];

figure;
errorbar(dN(1:5),dArea(1:5),errArea(1:5),'.','markersize',30,'color',[0.6350 0.0780 0.1840],'linewidth',2); hold on;
errorbar(dN(6:10),dArea(6:10),errArea(6:10),'.','markersize',30,'color',[0.4660 0.6740 0.1880],'linewidth',2);
errorbar(dN(11:15),dArea(11:15),errArea(11:15),'.','markersize',30,'color',[0 0.4470 0.7410],'linewidth',2);
plot(xLine,yLine,'-','color','k','LineWidth',3);
plot(xIdent,yIdent,'--','color',[0.5 0.5 0.5],'LineWidth',3);

%[0.6350 0.0780 0.1840]
%[0.4660 0.6740 0.1880]
%[0 0.4470 0.7410]
box on
ylabel('\DeltaArea/Area_{i} (Converted Region)')
xlabel('\DeltaN/N_{i} (Converted Cells)')
set(gca,'FontSize',18);

ylim([-0.1 1.5]);
xlim([-0.1 0.9]);

%% fig. 1i 
% (dLdt/Lreg vs total fraction GEM - reds)

%load "analysis_mat_WTdataSet.mat" from "Fig1" folder

analysis_mat = analysis_mat([analysis_mat.hpa]>72);

hpa = [analysis_mat.hpa];
L_reg = [analysis_mat.L_reg];
L_amp = [analysis_mat.L_amp];
phi = L_reg./L_amp;
dLdt = [analysis_mat.dLdt].*24;
total_fractionGEM = arrayfun(@(x)sum(x.nucleiGEM)./sum(x.nucleiALL),analysis_mat); 

% plot
xdata = total_fractionGEM;
ydata = adjust_neg(dLdt)./L_reg;
[~,idx] = remove_nan([xdata(:),ydata(:)]);
xdata = xdata(idx);
ydata = ydata(idx);
cdata = L_amp;
% cdata = hpa./24;

[~,idx] = remove_nan([xdata(:),ydata(:)]);
xdata = xdata(idx);
ydata = ydata(idx);
cdata = cdata(idx);

% fit 
[fit_obj,fit_obj_gof] = fit(xdata(:),ydata(:),fittype({'x'}));


f = figure();
%cm = colormap(flipud(cool()));
tLevels = 101;
color_map=colormap(cbrewer2('Reds',tLevels));
lampMinHere = 706;
lampMaxHere = 4053;
colorCode = ((cdata-700)./3353).*100+1;

scatter(xdata,ydata,50,round(colorCode),'filled',LineWidth=2)
%c = colorbar();
% c = colorbar('Ticks',0:0.5:1,'TickLabels', (somerangeyoulike));
color_map=colormap(cbrewer2('Reds',tLevels));
c = colorbar();
%set(g,'TickLabels',{'0','811','1621','2432','3242','4053'});%,'2740','3080','3420','3760','4100'})
%c = colorbar('Ticks',0:0.5:1,'TickLabels', (somerangeyoulike));
%set(c,'Ticks',0:0.5:1,'TickLabels',{'700','1380','2060','2740','3420','4100'});
set(c,'Ticks',10:17:100,'TickLabels',{'700','1380','2060','2740','3420','4100'});
c.Label.String = 'Length Amputated (um)';
%c.Label.String = 'Lamp(\mum)';
% c.Label.String = 'dpa';

plot_fit(fit_obj,linspace(0,0.8,1000),fit_obj_gof,confidence=0)
hold on;
p_id = plot(0:0.1:1,0:0.1:1,'--','color',[0.5 0.5 0.5], LineWidth=2);
hold off;

legend_f = findobj(f.Children,'Type','Legend');
legend_f.String{2} = 'y=x';
% legend(p_id,'y=x')

config_plot(gcf,c)
%ylabel('(dL/dt)/L (1/day)')
ylabel('\DeltaLength/Length_{Reg} (1/day)')
xlabel('Fraction GEM+')
xlim([0,1])
ylim([0,1])

%% Fig. 1j 
% (Plot sumGEM+ vs hpp - from Alvin - modified by Ashley (to mess with) 18 Apr 23 Red - 3 groups 13 Feb 24 - changeFit)

%load "para_mat_gem_merged" from "Fig1" folder

para_mat_gem_merged_copy = para_mat_gem_merged;
%para_mat_gem_merged_copy(22:26) = [];

fy = @(s)((s.sumnucleiGEM));
fx = @(s)s.hppTrue;
f = figure('visible','on');
%plot_parameters_merged3_long_copy2(fx,fy,para_mat_gem_merged(1:end-5),'Fraction Regen. (Lreg/Lamp)','Fraction Proliferating (GEM+)')
[~,~,fit_l,fit_m] = plot_parameters_merged3_long_copy5_red_3groupsChangeFit(fx,fy,para_mat_gem_merged_copy(1:end),'Time (hours post amputation)','# Osteoblasts Cycling');
% xlim([0,Inf])
xlim([36,336])
ylim([0,600]);
%saveas(f,[paths.plotFolder,filesep,'sumGEM_vs_hpp_SS000001_noScale_13feb24_red_changeFit_timeOffset48.png']);

%% Fig. 1k
% Plot FractionGEM+ vs Lreg/Lamp - Color Red - change fit - 3 groups

%load "para_mat_gem_merged" from "Fig1" folder

para_mat_gem_merged_copy = para_mat_gem_merged;

fy = @(s)((s.sumnucleiGEM)./(s.sumnucleiALL));
fx = @(s)s.Lr./s.L_amp;
%%%fLamp = @(s)s.L_amp;
f = figure('visible','on');
%plot_parameters_merged3_long_copy2(fx,fy,para_mat_gem_merged(1:end-5),'Fraction Regen. (Lreg/Lamp)','Fraction Proliferating (GEM+)')
%%%[~,~,fit_b,fit_m] = plot_parameters_merged3_long_copy2(fx,fy,para_mat_gem_merged_copy(1:end),'Fraction Regen. (Lreg/Lamp)','Fraction Proliferating (GEM+)');
[~,~,fit_b,fit_m] = plot_parameters_merged3_long_copy2_9Mar23_2newFit_3groups(fx,fy,para_mat_gem_merged_copy(1:end),'Fraction Regen. (Lreg/Lamp)','Fraction Osteoblasts Cycling');
% xlim([0,Inf])
xlim([0,1.2])
ylim([0,0.8]);
saveas(f,[paths.plotFolder,filesep,'FractionGEM_vs_LregByLamp_definedFit_red2_updateFit2_3colors2.png']);


%% Fig. 1l
% plot lambda for geminin verses Lreg/Lamp - Clean up and Fit Lambda

%load "analysis_mat_WTdataSet.mat" from "Fig1" folder

%step1
% Collect the data you will need - written 15Apr22
counter1 = 0;

BuildLambdaCollect = [];
BuildLambdaHere = [];

for i = 1:size(analysis_mat,2) %from first row to last row (which = size of para_mat_gem_merged in 2nd dimension)
%for i = 101
    fishHere = analysis_mat(i).fish;
    rayHere = analysis_mat(i).ray;
    hppHere = analysis_mat(i).hpa;%
%    hppTrueHere = analysis_mat(i).hppTrue;%
    %lambdaHere = analysis_mat(i).fit_gem_heavi.c;
    if ~isempty(analysis_mat(i).fit_gem_heavi)
        lambdaParse1 = analysis_mat(i).fit_gem_heavi;
        lambdaParse2 = coeffvalues(lambdaParse1);
        lambdaHere = lambdaParse2(3);
    else
        lambdaHere = NaN;
    end

    
    fractionGEM = analysis_mat(i).fractionGEM;
    if ~isempty(analysis_mat(i).fit_gem_heavi)
        slopeParse1 = analysis_mat(i).fit_gem_heavi;
        slopeParse2 = coeffvalues(slopeParse1);
        slopeHere = slopeParse2(1);
    else
        slopeHere = NaN;
    end
    LampHere = analysis_mat(i).L_amp;
    LregHere = analysis_mat(i).L_reg;
    
    lowGEM = min(fractionGEM);   
    
        if lowGEM <= 0.1
            if slopeHere <= 0.000101
                lambdaUSE = NaN;
            elseif lambdaHere > LregHere
                lambdaUSE = LregHere;
            else
            lambdaUSE = lambdaHere;
            end
        else 
            lambdaUSE = LregHere;
        end
        

    BuildLambdaHere = horzcat(fishHere,rayHere,hppHere,lambdaUSE,LregHere,LampHere,slopeHere);

    BuildLambdaCollect = vertcat(BuildLambdaCollect, BuildLambdaHere);
end

% Plot lambda/LengthRegen. v. Frac. Regen.

g = figure;
BuildLambdaCollectLat = [];
BuildLambdaCollectMed = [];
for m = 1:size(BuildLambdaCollect,1)
    if BuildLambdaCollect(m,2) < 4
        plot((BuildLambdaCollect(m,5)./BuildLambdaCollect(m,6)),(BuildLambdaCollect(m,4)./BuildLambdaCollect(m,5)),'.','MarkerSize',10,'Color',[0.3010 0.7450 0.9330]); hold on;
        BuildLambdaCollectLat = vertcat(BuildLambdaCollectLat,BuildLambdaCollect(m,:));
    elseif BuildLambdaCollect(m,2) > 4
        plot((BuildLambdaCollect(m,5)./BuildLambdaCollect(m,6)),(BuildLambdaCollect(m,4)./BuildLambdaCollect(m,5)),'.','MarkerSize',10,'Color',[1 .6 1]); hold on;
        BuildLambdaCollectMed = vertcat(BuildLambdaCollectMed,BuildLambdaCollect(m,:));
    end
end


%step 2
%plot lambda for geminin verses Lreg/Lamp - Clean up and Fit Lambda - color code data points red by Lamp - edit 2.16.24 - 3 groups
% fit  w/ a*x^b*exp(-c*x)
%fit w/ F(x) = e^alpha*(beta/alpha)^alpha*x^alpha*e^(-beta*x)

%set up color map
tLevels = 101;
color_map=colormap(cbrewer2('Reds',tLevels));
lampMinHere = 706;
lampMaxHere = 4053;

%use BuildLambdaCollectLat and BuildLambdaCollectMed from previous code
%block

BuildLambdaCollectHere = vertcat(BuildLambdaCollectMed, BuildLambdaCollectLat);

% Plot lambda/LengthRegen. v. Frac. Regen.

figure;

for m = 1:size(BuildLambdaCollectHere,1)
    colorCode = ((BuildLambdaCollectHere(m,6)-706)./3353).*100+1; %%% 699 will need to be edited 
    plot((BuildLambdaCollectHere(m,5)./BuildLambdaCollectHere(m,6)),(BuildLambdaCollectHere(m,4)./BuildLambdaCollectHere(m,5)),'.','MarkerSize',20,'Color',color_map((round(colorCode)),:)); hold on;        
end

%fit all data
idxData = ~isnan(BuildLambdaCollectHere(:,4));
LregData = BuildLambdaCollectHere(idxData,5);
LampData = BuildLambdaCollectHere(idxData,6);
lambdaData = BuildLambdaCollectHere(idxData,4);

%subsetData
dataSmall = [];
dataMed = [];
dataLar = [];
    for qqq = 1:size(LampData)
        if LampData(qqq) < 1575
            dataHere = horzcat(LregData(qqq),LampData(qqq),lambdaData(qqq));
            dataSmall = vertcat(dataSmall,dataHere);
        elseif LampData(qqq)<3089
            dataHere = horzcat(LregData(qqq),LampData(qqq),lambdaData(qqq));
            dataMed = vertcat(dataMed,dataHere);
        else 
            dataHere = horzcat(LregData(qqq),LampData(qqq),lambdaData(qqq));
            dataLar = vertcat(dataLar,dataHere);
        end
    end

% ft = fittype( 'exp1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';

    ft = fittype( 'exp(a)*(b/a)^(a)*x^(a)*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.3 5];

%[ft, gof_lat] = fit((LregLat./LampLat),(lambdaLat./LregLat),'poly1');
[ft, gof_all] = fit((LregData./LampData),(lambdaData./LregData),ft,opts);

%[ft_all, gof_all] = fit((LregData./LampData),(lambdaData./LregData),ft,opts);
[ft_lar,gof_lar] = fit(dataLar(:,1)./dataLar(:,2),(dataLar(:,3)./dataLar(:,1)),ft,opts);
[ft_med,gof_med] = fit(dataMed(:,1)./dataMed(:,2),(dataMed(:,3)./dataMed(:,1)),ft,opts);
[ft_small,gof_small] = fit(dataSmall(:,1)./dataSmall(:,2),(dataSmall(:,3)./dataSmall(:,1)),ft,opts);


hold on;
fit_plot = plot(ft);
set(fit_plot,'linewidth',2,'color','k'); hold on;

fit_plot_lar = plot(ft_lar);
set(fit_plot_lar,'linewidth',2,'color',[0.5 0 0]);

fit_plot_med = plot(ft_med);
set(fit_plot_med,'linewidth',2,'color',[1 0 0]);

fit_plot_small = plot(ft_small);
set(fit_plot_small,'linewidth',2,'color',[1 0.5 0.5]);

%scale plot
xlabel('Fraction Regenerated');
ylabel('Fraction of Regenerate Occupied by Proliferative Zone');
xlim([0 1.2]);
ylim([0 1.2]);
legend('off');
set(gca, 'fontsize', 15,'linewidth',3'); 

legend_name={strcat("R^2 = ",num2str(gof_all.rsquare)),...
    strcat("R^2 = ",num2str(gof_lar.rsquare))...
    strcat("R^2 = ",num2str(gof_med.rsquare))...
    strcat("R^2 = ",num2str(gof_small.rsquare))};
legend([fit_plot,fit_plot_lar,fit_plot_med,fit_plot_small], legend_name,...
        'Location','best','color','none','box','off');
