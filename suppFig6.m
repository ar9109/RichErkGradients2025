%% Supplemental Figure 6

%% Supp. Fig. 6A
%load "analysis_mat_dbin150.mat" from "SuppFig4" folder

%step 1 - create analysis matrix that only contains 120 hpa (pre) and 144 hpa (24 hr post) *

analysis_mat = analysis_mat_dbin150;


analysis_mat_trim = [];
for i = 1:size(analysis_mat,2)
    if analysis_mat(i).hpa == 120 || analysis_mat(i).hpa == 132
        analysis_mat_trim = vertcat(analysis_mat_trim,analysis_mat(i));
    else
        continue
    end
    
end

%step 2 - create treatment mat *
treatmentMat = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1];

%step 3 - Plot difference at each position *

analysis_mat_combined = analysis_mat_trim;

%fishSet = [1,2,3,4,5,6,7,8,11,12,13,14,15,16,18];
preVals = linspace(1,15,8);
postVals = linspace(2,16,8);

%h4 = figure;
diffTreat = [];
diffCont = [];
xTreat = [];
xCont = [];

for i = 1:size(preVals,2)

    q=preVals(i);
    j=postVals(i);
    treatment = treatmentMat(q);

    KTRi = analysis_mat_combined(q).ktr;
    KTRf = analysis_mat_combined(j).ktr;
    centersi = analysis_mat_combined(q).ccrot;
    centersf = analysis_mat_combined(j).ccrot;
    centersXi = analysis_mat_combined(q).ccrot(:,1);
    centersXf = analysis_mat_combined(j).ccrot(:,1);
    
    maxXi = max(centersXi);
    minXi = min(centersXi);
    flipCentersXi = (centersXi-maxXi)*-1;
    flipMaxXi = max(flipCentersXi);
    flipMinXi = min(flipCentersXi);
    normCentersXi = (flipCentersXi-flipMinXi)/(flipMaxXi-flipMinXi);

    maxXf = max(centersXf);
    minXf = min(centersXf);
    flipCentersXf = (centersXf-maxXf)*-1;
    flipMaxXf = max(flipCentersXf);
    flipMinXf = min(flipCentersXf);
    normCentersXf = (flipCentersXf-flipMinXf)/(flipMaxXf-flipMinXf);
    
    dbins=5;
    
    [yMeani,xMeani,ySEMi,yNi] = bin_average(normCentersXi,KTRi,dbins,min(normCentersXi),max(normCentersXi));

    %set up plot label
    fish = analysis_mat_combined(q).fish;
    fish2 = analysis_mat_combined(j).fish;
    ray = analysis_mat_combined(q).ray;
    ray2 = analysis_mat_combined(j).ray;
    hpa = analysis_mat_combined(q).hpa;
    hpa2 = analysis_mat_combined(j).hpa;

    [yMeanf,xMeanf,ySEMf,yNf] = bin_average(normCentersXf,KTRf,dbins,min(normCentersXf),max(normCentersXf));

    diffHere = yMeanf-yMeani;

    if treatment == 1
        diffTreat = vertcat(diffTreat,diffHere');
        xTreat = vertcat(xTreat,xMeani');
    elseif treatment == 0
        diffCont = vertcat(diffCont,diffHere');
        xCont = vertcat(xCont,xMeani');
    end
    
end

%Prep data for box plot
plotData = horzcat(diffCont(:,1),diffTreat(:,1),diffCont(:,2),diffTreat(:,2),diffCont(:,3),diffTreat(:,3),diffCont(:,4),diffTreat(:,4),diffCont(:,5),diffTreat(:,5));


xMat = zeros(size(plotData,1),1);



h4 = figure;
for i = 1:size(diffCont,2)*2
    if i == 1 || i == 3 || i    == 5 || i == 7 || i == 9
        cl = [0 0 0 0.5];
    elseif i == 2 || i == 4 || i == 6 || i == 8 || i == 10
        cl = [1 0 0 0.5];
    end
    boxchart(xMat+i,plotData(:,i),'BoxFaceColor',cl,'markercolor',cl); hold on;
end

xticks([1 2 3 4 5 6 7 8 9 10])
xticklabels({'Amp. Site','Amp. Site','1/4 Regen.','1/4 Regen.','Middle','Middle','3/4 Regen.','3/4 Regen.','Tip','Tip'})

ylabel('Difference in Avg. Erk Activity (Pre-Post)')
xlabel('Position')
set(gca,'fontsize',16)

%% Supp. Fig. 6B
%% Plot change in integrated density - use this one

close all
figure;

xVals = [ 2 2 2 2 2 2 1 1 1 1 1 1];
yVals = [33300.276 187217.631 412525.253 128668.503 122534.436 200587.695 677157.025 472295.684 488427.916 361966.942 338343.435 675032.14];

[h,p,ci,stats] = ttest2(yVals(7:12),yVals(1:6));

% yValsCyclo = [33300.276 187217.631 412525.253 128668.503 122534.436 200587.695];
% yValsWater = [677157.025 472295.684 488427.916 361966.942 338343.435 675032.14];
% 
% [h2,p2,ci2,stats2] = ttest2(yVals(7:12),yVals(1:6));

%bar(xVals,yVals)
boxchart(xVals,yVals); hold on
plot(xVals,yVals,'.','MarkerSize',20,'color','k')

xticks([1 2])
xticklabels({'+H2O','+Cyclo'})
%xtickangle(45)

xstat = [1,2];
ystat = [690000,690000];
plot(xstat,ystat,'-','color','k','linewidth',1)
txt = ['p = ' num2str(round(p,1,"significant"))];

text(1.4,670000,txt,'FontSize',16)

ylabel('Change in Int. Density')

set(gca,'fontsize',16)

%yline(200000,'--')
%xline(6.5,'--')