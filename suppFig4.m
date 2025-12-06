%% Supp. Fig. 4A

%load "analysis_mat_movquant_8h.mat" from "SuppFig4"

% plot quantile correlation - 25th percentile v. average
close all

xAll = cat(2,analysis_mat_movquant.averageKTR);
%%%yAll = cat(2,analysis_mat_movquant.stdKTR);
yAll = cat(2,analysis_mat_movquant.quant25KTR);
% yAll = cat(2,analysis_mat_movquant.quant10KTR);
cAll = cat(2,analysis_mat_movquant.binvalue);
corr(xAll',yAll')
xAll = xAll -0.8;
yAll = yAll -0.8;
[fit_obj,gof] = fit(xAll',yAll',{'x','1'});
f = figure;
%%%cm = colormap(flipud(cool(numel(s.ktr_bindata))));
cm = colormap(flipud(cool(numel(cAll))));
c = colorbar;
c.Label.String = 'Distance from Amp. Site (\mum)';
hold on;
%scatter(xAll,yAll,30,cAll,LineWidth=2);
scatter(xAll,yAll,30,cAll,LineWidth=2);
%pf = plot_fit(fit_obj,xAll,gof,confidence=0);
pf = plot_fit(fit_obj,xAll,gof,confidence=0);
%xlim([0,0.4]);
%ylim([0,0.5]);
xlabel('Average Erk Activity (A.U.)');
%ylabel('25th Percentile Value (Erk Activity)')
ylabel('25th Percentile of Erk Activity (A.U.)')
config_plot(f,c)

%% Supp. Fig. 4B

%load "analysis_mat_movquant_8h.mat" from "SuppFig4"

% plot quantile correlation - std v. average
close all

xAll = cat(2,analysis_mat_movquant.averageKTR);
yAll = cat(2,analysis_mat_movquant.stdKTR);
%%%yAll = cat(2,analysis_mat_movquant.quant25KTR);
% yAll = cat(2,analysis_mat_movquant.quant10KTR);
cAll = cat(2,analysis_mat_movquant.binvalue);
corr(xAll',yAll')
xAll = xAll -0.8;
[fit_obj,gof] = fit(xAll',yAll',{'x','1'});
f = figure;
%%%cm = colormap(flipud(cool(numel(s.ktr_bindata))));
cm = colormap(flipud(cool(numel(cAll))));
c = colorbar;
c.Label.String = 'Distance from Amp. Site (\mum)';
hold on;
%scatter(xAll,yAll,30,cAll,LineWidth=2);
scatter(xAll,yAll,30,cAll,LineWidth=2);
%pf = plot_fit(fit_obj,xAll,gof,confidence=0);
pf = plot_fit(fit_obj,xAll,gof,confidence=0);
% xlim([0.7,1.21]);
xlim([0,0.4]);
% ylim([0.6,1.1]);
ylim([0,0.5]);
xlabel('Average Erk Activity (A.U.)');
%ylabel('25th Percentile Value (Erk Activity)')
ylabel('Standard Deviation of Erk Activity (A.U.)')
config_plot(f,c)

%% Supp. Fig. 4C
% plot osteoblast density in space

%load 1 of the following 4 matrices from "Supp Fig 4)
    %"analysis_mat_0to4_updateAug25.mat"
    %"analysis_mat_3to5_updateAug25_2.mat"
    %"analysis_mat_5to7_updateAug25.mat"
    %"analysis_mat_7to14_updateAug25.mat"

% Plot internuclear distance - drop bins with less than 35 nuclei

close all
figure;

analysis_mat = analysis_mat_0to4_updateAug25; %You will need to change 'analysis_mat_0to4_updateAug25' to the name of the matrix you have loaded

binCollect = [];
distCollect = [];

%find max number of bins
max = 0;
for i = 1:size(analysis_mat,2)
    sizeHere = size(analysis_mat(i).binvalue,2);
    if sizeHere > max
        max = sizeHere;
    else 
        continue
    end
end


for i = 1:size(analysis_mat,2)
    distScreen = [];
    idxHere = analysis_mat(i).nucleiNum > 34;
    for qq = 1:size(idxHere,2)
        if idxHere(qq) == 1
            distLocal = analysis_mat(i).interspace_mean5min(qq);
        elseif idxHere(qq) == 0 
            distLocal = NaN;
        end
    distScreen = horzcat(distScreen,distLocal);
    end

    plot(analysis_mat(i).binvalue,1./distScreen,...
        '-','Color', [0    0.4470    0.7410 0.3],'MarkerSize',10,'LineWidth',1);
    sizeBinMat = size(analysis_mat(i).binvalue,2);
    zeroMat = zeros(1,max-sizeBinMat);
    zeroMat(zeroMat == 0) = NaN;
    
    binColHere = horzcat(analysis_mat(i).binvalue,zeroMat);
    binCollect = vertcat(binCollect,binColHere);

    distColHere = horzcat(1./distScreen,zeroMat);  
    distCollect = vertcat(distCollect,distColHere);
    hold on;
end

%calc mean
semDist = [];
meanBin = nanmean(binCollect);
meanDist = nanmean(distCollect);
sdDist = nanstd(distCollect);
distIdx = ~isnan(distCollect);
numDist = sum(distIdx);
for q = 1:size(sdDist,2)
    semLoc = sdDist(q)./sqrt(numDist(q));
    semDist = horzcat(semDist,semLoc);
end

errorbar(meanBin,meanDist,semDist,'.-','markersize',10,'LineWidth',2,'Color',[0.8500    0.3250    0.0980])
ylim([0,0.1])

xlabel("Position (\mum)")
  

ylabel("1/Average Internuclear Distance (1/\mum)")
%%%title({'7-14 dpa Rays (n=110)';'Drop Bins w/ few Nuc.'})
set(gca,'fontsize',16)

%% Supp. Fig. 4D

%load data matrix:
%"/Users/ashleyrich/Desktop/analysis_mat_timeaverage_WTdataSet.mat" from
%"SuppFig4"

%step 1 - collect predicted proliferation profiles *

L_amp = [analysis_mat_timeaverage.L_amp];
L_reg = [analysis_mat_timeaverage.L_reg];
phi = L_reg./L_amp;
phiIdx = phi < 1;
phi_trim = phi(phiIdx);
gem_u = [analysis_mat_timeaverage.fractionGEM_u_timeave];
gem_u_trim = gem_u(:,phiIdx);
gem_u_flip = flip(gem_u,1);
gem_u_flip_trim = gem_u_flip(:,phiIdx);

fish = [analysis_mat_timeaverage.fish];
fish_trim = fish(phiIdx);
ray = [analysis_mat_timeaverage.ray];
ray_trim = ray(phiIdx);
hpa = [analysis_mat_timeaverage.hpa];
hpa_trim = hpa(phiIdx);

A0 = 0.35998; % A(phi)
af = 1.2; % f(u)
bf = 0.4; % f(u)
u = 0.05:0.1:0.95;
%ERK = A0.*(1-phi).*(af.*(u)+bf);

%alpha = 2.4097; % G(u) vs E(u)
alpha = 2.3981; % G(u) vs E(u)
%a2 = 7.06473694; % G(u) vs E(u)
a2 = 7.1977; % G(u) vs E(u)
%Prolif = a2*ERK.^alpha;

ProlifAll = [];
for i = 1:size(phi_trim,2)
%for i = 145
    ERK = [];
    for j = 1:size(u,2)
        ERKj = A0.*(1-phi_trim(i)).*(af.*(u(j))+bf);
        ERK = vertcat(ERK,ERKj);
    end
    
    Prolif = [];
    for k = 1:size(ERK)
        Prolif_i = a2*ERK(k).^alpha;
        Prolif = vertcat(Prolif,Prolif_i);
    end
    ProlifAll = horzcat(ProlifAll,Prolif);
    %test = 1;
end

% step 2 - reshape gem actual & gem predicted matrices *
gemActResh = [];
for i = 1:size(gem_u_flip_trim,2)
    valHere = gem_u_flip_trim(:,i);
    gemActResh = vertcat(gemActResh,valHere);
end

gemPredResh = [];
for j = 1:size(ProlifAll,2)
    valHere2 = ProlifAll(:,j);
    gemPredResh = vertcat(gemPredResh,valHere2);
end

hpaResh = [];
for k = 1:10
    hpaResh=vertcat(hpaResh,hpa_trim);
end

hpaResh2 = [];
for l = 1:size(hpaResh,2)
    valHere3 = hpaResh(:,l);
    hpaResh2 = vertcat(hpaResh2,valHere3);
end

% step 3 plot gem actual v. gem predicted - drop data from > 96 hpa - flip axes *

hpaIdx = hpaResh > 96;
gemActResh_trim = gemActResh(hpaIdx);
gemPredResh_trim = gemPredResh(hpaIdx);
hpaResh2_trim = hpaResh2(hpaIdx);

% Set up fittype and options.
%ft = fittype( 'poly1' );
yData = gemActResh_trim;
xData = gemPredResh_trim;
yDataIDX = ~isnan(yData);
xData = xData(yDataIDX);
yData = yData(yDataIDX);

%[yMean,xMean,ySEM] = bin_average(xData,yData,12,0,0.6064);
[yMean,xMean,ySEM] = bin_average(xData,yData,12,1.3655e-05,1);
%[yMean,xMean,ySEM] = bin_average(xData,yData,12,0,1);

f = figure;

for i = 1:size(gemActResh_trim,1)
    if hpaResh2_trim(i) < 108
        cl = [0.8500 0.3250 0.0980];
    else
        cl = [0 0.4470 0.7410];
    end
    %i
    %hpaResh2(i)
    plot(gemPredResh_trim(i),gemActResh_trim(i),'.','markersize',15,'color',cl); hold on;
end

% % % % Set up fittype and options.
ft = fittype( 'poly1' );
% % % xData = gemActResh_trim;
% % % yData = gemPredResh_trim;
% % % xDataIDX = ~isnan(xData);
% % % xData = xData(xDataIDX);
% % % yData = yData(xDataIDX);

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

xIdent = linspace(0,1,10);
yIdent = linspace(0,1,10);

mm = plot(xIdent,yIdent,'--','color',[0.5 0.5 0.5]);

% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
h = plot(fitresult);
set(h,'color','k','linewidth',2)

g = errorbar(xMean,yMean,ySEM,'.-','markersize',20,'linewidth',3,'color',[0.8500    0.3250    0.0980]);

xplot = [1 2 3 4 5 6 7 8 9 10];
yplot = [1 2 3 4 5 6 7 8 9 10];

legend([h,g,mm],{['R^2=' num2str(gof.rsquare)],'Binned Average','y=x'})

xlabel('Frac. Gem + (Predicted)')
ylabel('Frac. Gem + (Measured)')
set(gca,'fontsize',16)

ylim([0 1])
xlim([0 1])

% paths.plotFolder = '/Users/ashleyrich/Documents/DiTaliaLab/Analysis/5_2025_revisionWork/predGem_v_actGem';
% mkdir(paths.plotFolder);
% savename = ['gemPred_v_normSpace_over96hpa_flipAxis_2Sept25'];
% exportgraphics(f,[paths.plotFolder,filesep,savename,'.png']);



