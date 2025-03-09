%% Fig. 2B (run all code blocks until Fig. 2D)

%step1
%import data in excel file,
%"1_9_23_MekInhibition_Measurement_Collection.xlsx" in "Fig2" folder
%Adjust code below for particular file location

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 8);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:H32";

% Specify column names and types
opts.VariableNames = ["Experiment", "Fish", "Ray", "Pre", "Post", "Treatment", "umChange", "Change"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
MekInhibitionMeasurementCollection = readtable("/Users/ashleyrich/Documents/DiTaliaLab/Experiments/1_9_23_MekInhibition_Measurements_Collection/1_9_23_MekInhibition_Measurement_Collection.xlsx", opts, "UseExcel", false);
mekData = MekInhibitionMeasurementCollection;

%step2
% Convert to output type
mekData = table2array(mekData);
clear opts

%step3
% Organize Data

%subsetData
conData = [];
treatData = [];

for i = 1:size(mekData,1)
    if mekData(i,6) == 0
        conData = vertcat(conData,mekData(i,:));
    elseif mekData(i,6) == 1
        treatData = vertcat(treatData,mekData(i,:));
    end
end

%generate averages

conAvg = nanmean(conData(:,7));
conStd = nanstd(conData(:,7));
conSem = conStd./sqrt(numel(conData));

treatAvg = nanmean(treatData(:,7));
treatStd = nanstd(treatData(:,7));
treatSem = treatStd./sqrt(numel(treatData));

%step4
% Plot data
close all;


cHere = [0 0.4470 0.7410];
scatter0 = [-.5 -.4 -.3 -.2 -.1 0 .1 .2 .3 .4 .5];
scatter1 = [-.6 -.5 -.4 -.3 -.2 -.1 0 .1 .2 .3 .4 .5 .6 .7];
scatter0 = scatter0./4;
scatter1 = scatter1./4;

boxData = vertcat(conData(:,6:7),treatData(:,6:7));
xgroupData = vertcat(conData(:,6),treatData(:,6));
yData = vertcat(conData(:,7),treatData(:,7));

%stats
[h,p,ci,stats] = ttest2(conData(:,7),treatData(:,7));

figure; 

b = boxchart(xgroupData,yData);
b.BoxFaceColor = [0 0 0];
b.BoxFaceAlpha = 0;
%b.MarkerColor = 'k';
b.MarkerColor = [1 1 1];
hold on;
plot(scatter0 + 0,conData(:,7),'.','color',cHere,'Markersize',20); hold on;
plot(scatter1 + 1,treatData(:,7),'.','color',cHere,'Markersize',20); hold on;
box on;

plot([0,1],[220,220],'-','linewidth',1,'color','k');
t = text(0.2,230,['\it \fontsize{16} p=' num2str((round(p,11)))]);

ylabel('Growth (\mum)')
xlim([-1,2]);
%ylim([-60,120]);
set(gca,'FontSize',18);

ax = gca;
ax.XTick = [0, 1];
ax.XTickLabels = {'+DMSO','+PD03'};

%% Fig. 2D 

%step1
%import data in excel file,
%"1_9_23_MekInhibition_Measurement_Collection.xlsx" in "Fig2" folder
%Adjust code below for particular file location


% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 13);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A3:M32";

% Specify column names and types
opts.VariableNames = ["Experiment", "Fish", "Ray", "Pre", "Post", "Treatment", "umChange", "Change", "VarName9", "GemPos_pre", "TotalCell_pre", "GemPos_post", "TotalCell_post"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
MekInhibitionMeasurementCollection = readtable("/Users/ashleyrich/Documents/DiTaliaLab/Experiments/1_9_23_MekInhibition_Measurements_Collection/1_9_23_MekInhibition_Measurement_Collection.xlsx", opts, "UseExcel", false);

clear opts
mekData = MekInhibitionMeasurementCollection;
mekData = table2array(mekData);

%step2
% 3.28.24 - collect data for below plots
dataAll = [];
for i = 1:size(mekData)
    fracPre = mekData(i,10)./mekData(i,11);
    fracPost = mekData(i,12)./mekData(i,13);
    numPre = mekData(i,10);
    numPost = mekData(i,12);
    fish = mekData(i,2);
    ray = mekData(i,3);
    treat = mekData(i,6);
    dataHere = horzcat(fish,ray,numPre,numPost,fracPre,fracPost,treat);
    dataAll = vertcat(dataAll,dataHere);
end

%step3
% Plot #gem+

figure;

conDataHere = [];
treatDataHere = [];
for i = 1:size(dataAll)
    if dataAll(i,7) < 1
        plot([0.1,0.2],[dataAll(i,3),dataAll(i,4)],'.-','MarkerSize',20,'linewidth',2,'color',[0 0.4470 0.7410]); hold on;
        dataHere = horzcat(dataAll(i,3),dataAll(i,4));
        conDataHere = vertcat(conDataHere,dataHere);
    elseif dataAll(i,7) > 0
        plot([0.3,0.4],[dataAll(i,3),dataAll(i,4)],'.-','MarkerSize',20,'linewidth',2,'color',[0 0.4470 0.7410]); hold on;
        dataHere = horzcat(dataAll(i,3),dataAll(i,4));
        treatDataHere = vertcat(treatDataHere,dataHere);
    end
end

hold on;
plot([0.1,0.2],[mean(conDataHere(:,1)),mean(conDataHere(:,2))],'.-','markersize',40,'linewidth',4,'color','k')
plot([0.3,0.4],[mean(treatDataHere(:,1)),mean(treatDataHere(:,2))],'.-','markersize',40,'linewidth',4,'color','k')

[hCon,pCon] = ttest(conDataHere(:,1),conDataHere(:,2));
[hTreat,pTreat] = ttest(treatDataHere(:,1),treatDataHere(:,2));

plot([0.1,0.2],[310,310],'-','linewidth',1,'color','k');
t = text(0.1,320,['\it \fontsize{16} p=' num2str(round(pCon,1))]);

plot([0.3,0.4],[310,310],'-','linewidth',1,'color','k');
t = text(0.3,320,['\it \fontsize{16} p=' num2str(round(pTreat,9))]);

ylabel('# Cycling Cells');
ylim([0,340]);
xlim([0,0.5]);


ax0 = gca;
ax0.XTick = [0.1, 0.2,0.3,0.4];
ax0.XTickLabels = {'0 hpt','24 hpt','0 hpt','24 hpt'};

set(gca,'fontsize',16);


%% Fig. 2G

%load "analysis_mat_timeaverage_WTdataSet.mat" from "Fig2" folder

% G(x) vs E(x) - this is actually u, not x (RED) - 3 groups

analysis_mat = analysis_mat_timeaverage;

close all;
f = figure();
[c,fit_obj,fit_obj_gof] = plot_Gx_by_Ex_noODR_red_3groups(analysis_mat,100,1); % use this line for figure in paper
config_plot(f,c);



