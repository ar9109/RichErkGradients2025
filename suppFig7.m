%% Supp. Fig. 7A
% Plot GFP domain length / Lamp v. Lamp - here
chToResize = [1]; %reference channel to resize 

%paths.objFolder is located in "SuppFig7" folder under the name
%"objects_old_22Apr24_new_copy_no0"
%Please adjust the paths.masterFolder and paths.objFolder lines (12 & 13)
%to where you store the data


paths=[];

paths.masterFolder='/Users/ashleyrich/Documents/submissionCode/SuppFig7/'; %folder where data is stored
paths.objFolder= [paths.masterFolder 'objects_old_22Apr24_new_copy_no0/']; %objects folder

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='72';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

opts=[];
%opts.resizeFactor = 0.25; %resizing factor

optsSave=[];
optsSave.compress='lzw';
optsSave.overwrite=true;

h = figure;

lampCollect = [];
lGFPcollect = [];

for i = 1:size(st_dir)
    display([paths.objFolder filesep st_dir(i).name]);
    load([paths.objFolder filesep st_dir(i).name]);
    fish = myScale.fish;
    ray = myScale.ray;
    hpa = myScale.hpp;
    lGFP = myScale.GFPlength;
    lamp = myScale.lamp;

    if hpa > 70
        if fish < 5
            cl = 'k';
            plot(lamp,lGFP,'.','MarkerSize',15,'color',cl); hold on;
            lampCollect = vertcat(lampCollect,lamp);
            lGFPcollect = vertcat(lGFPcollect,lGFP);
%         elseif fish <20
%             cl = 'r';
        elseif fish > 20
            cl = 'k';
            plot(lamp,lGFP,'.','MarkerSize',15,'color',cl); hold on;
            lampCollect = vertcat(lampCollect,lamp);
            lGFPcollect = vertcat(lGFPcollect,lGFP);
        end

        %plot(lamp,lGFP./lamp,'.','MarkerSize',15,'color',cl); hold on;
    else
        continue
    end
end

% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit( lampCollect, lGFPcollect, ft );
% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
%h = 
f = plot(fitresult);
set(f,'linewidth',2,'color','k')
%legend(f, 'R^{2}', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );

legend_name={strcat("R^2 = ",num2str(round(gof.rsquare,1)))};

legend(f, legend_name,...
        'Location','best','color','none','box','off');

xlabel('Length Amputated (\mum)');
ylabel({'Length Fgf20a Expression Domain (\mum)'})
%ylabel('Length Fgf20a Expression Domain (\mum) / Length Amputated (\mum)');
set(gca,'fontsize',16);

%% Supp. Fig. 7B
% Box Plot - Actual Fraction Regenerated

% load "fgf20a_exp1_2_allWTdata_Trim.mat" and "fgf20a_exp2_allLOFdata.mat"
% from "SuppFig7" folder

% Subset Data into WTlat,med and LOF lat,med - Ashley use

%subsetLOF data

LOFlat = [];
LOFmed = [];
for n = 1:size(exp2DataLOFNew)
    if exp2DataLOFNew(n,3) < 4
        LOFlat = vertcat(LOFlat,exp2DataLOFNew(n,:));
    elseif exp2DataLOFNew(n,3) > 4
        LOFmed = vertcat(LOFmed,exp2DataLOFNew(n,:));
    end
end

%subsetWT data

WTlat = [];
WTmed = [];
for n = 1:size(allWTnew)
    if allWTnew(n,3) < 4 || allWTnew(n,3) > 10
        WTlat = vertcat(WTlat,allWTnew(n,:));
    elseif allWTnew(n,3) > 4
        WTmed = vertcat(WTmed,allWTnew(n,:));
    end
end

close all;

Lfin_Lamp_WTlat = (WTlat(:,7))./(WTlat(:,6));
Lfin_Lamp_WTmed = (WTmed(:,7))./(WTmed(:,6));
Lfin_Lamp_LOFlat = (LOFlat(:,7))./(LOFlat(:,6));
Lfin_Lamp_LOFmed = (LOFmed(:,7))./(LOFmed(:,6));

WTlatLabel = zeros(size(Lfin_Lamp_WTlat));
WTmedLabel = zeros(size(Lfin_Lamp_WTmed))+1;
LOFlatLabel = zeros(size(Lfin_Lamp_LOFlat))+2;
LOFmedLabel = zeros(size(Lfin_Lamp_LOFmed))+3;

%box chart script
xgroupData2 = vertcat(WTlatLabel,WTmedLabel,LOFlatLabel,LOFmedLabel);
yData2 = vertcat(Lfin_Lamp_WTlat,Lfin_Lamp_WTmed,Lfin_Lamp_LOFlat,Lfin_Lamp_LOFmed);

figure;
b = boxchart(xgroupData2,yData2);
b.BoxFaceColor = [0 0 0];
b.BoxFaceAlpha = 0;
%b.MarkerColor = 'k';
b.MarkerColor = [1 1 1];
hold on;

%cHere = [0 0.4470 0.7410];
scatter0 = [-.5 -.4 -.3 -.2 -.1 0 .1 .2 .3 .4 .5];
scatter1 = [-.6 -.5 -.4 -.3 -.2 -.1 0 .1 .2 .3 .4 .5 .6 .7];
scatter0 = scatter0./4;
scatter1 = scatter1./4;

scatter0 = linspace(-0.25,0.25,size(WTlatLabel,1));
scatter1 = linspace(0.75,1.25,size(WTmedLabel,1));
scatter2 = linspace(1.75,2.25,size(LOFlatLabel,1));
scatter3 = linspace(2.75,3.25,size(LOFmedLabel,1));

plot(scatter0,Lfin_Lamp_WTlat,'.','color',[0 0 0],'Markersize',20); hold on;
plot(scatter1,Lfin_Lamp_WTmed,'.','color',[0 0 0],'Markersize',20); hold on;
plot(scatter2,Lfin_Lamp_LOFlat,'.','color',[0.6350 0.0780 0.1840],'Markersize',20); hold on;
plot(scatter3,Lfin_Lamp_LOFmed,'.','color',[0.6350 0.0780 0.1840],'Markersize',20); hold on;
box on;





%Linf stats work
[p2_WT_lat_med,h2_WT_lat_med,stats2_WT_lat_med] = ranksum(Lfin_Lamp_WTlat,Lfin_Lamp_WTmed); %p = 0.1255
[p2_LOF_lat_med,h2_LOF_lat_med,stats2_LOF_lat_med] = ranksum(Lfin_Lamp_LOFlat,Lfin_Lamp_LOFmed); %p = 0.0095
[p2_WTlat_LOFlat,h2_WTlat_LOFlat,stats2_WTlat_LOFlat] = ranksum(Lfin_Lamp_WTlat,Lfin_Lamp_LOFlat); %p = 0.0022
[p2_WTmed_LOFmed,h2_WTmed_LOFmed,stats2_WTmed_LOFmed] = ranksum(Lfin_Lamp_WTmed,Lfin_Lamp_LOFmed); %p = 0.3277

plot([0,1],[2.15,2.15],'-','color','k','linewidth',1);
%txt1 = ['{\it p=}' num2str(p2_WT_lat_med,'%.4f')];
txt1 = ['{\it p=}' num2str(round(p2_WT_lat_med,5))];
text(0.2,2.25,txt1,'fontsize',16,'color',[0 0 0])

plot([2,3],[2.15,2.15],'-','color','k','linewidth',1);
%txt2 = ['{\it p=}' num2str(p2_LOF_lat_med,'%.4f')];
txt2 = ['{\it p=}' num2str(round(p2_LOF_lat_med,3))];
text(2.2,2.25,txt2,'fontsize',16,'color',[0 0 0])

plot([0,2],[2.45,2.45],'-','color','k','linewidth',1);
%txt3 = ['{\it p=}' num2str(p2_WTlat_LOFlat,'%.4f')];
txt3 = ['{\it p=}' num2str(round(p2_WTlat_LOFlat,3))];
text(0.5,2.55,txt3,'fontsize',16,'color','k')

plot([1,3],[2.75,2.75],'-','color','k','linewidth',1);
%txt3 = ['{\it p=}' num2str(p2_WTmed_LOFmed,'%.4f')];
txt3 = ['{\it p=}' num2str(round(p2_WTmed_LOFmed,4))];
text(1.5,2.85,txt3,'fontsize',16,'color','k')



ylabel({'Length_{final}/Length_{amp}'})
%xlim([-1,2]);
ylim([0,3]);
set(gca,'FontSize',18);

ax = gca;
ax.XTick = [0, 1, 2, 3];
ax.XTickLabels = {'WT Lat.','WT Med.','LOF Lat.','LOF Med.'};
xtickangle(60)

%% Supp. Fig. 7C

% load "fgf20a_exp1_2_allWTdata_Trim.mat" and "fgf20a_exp2_allLOFdata.mat"
% from "SuppFig7" folder

% Subset Data into WTlat,med and LOF lat,med - Ashley use

%subsetLOF data

LOFlat = [];
LOFmed = [];
for n = 1:size(exp2DataLOFNew)
    if exp2DataLOFNew(n,3) < 4
        LOFlat = vertcat(LOFlat,exp2DataLOFNew(n,:));
    elseif exp2DataLOFNew(n,3) > 4
        LOFmed = vertcat(LOFmed,exp2DataLOFNew(n,:));
    end
end

%subsetWT data

WTlat = [];
WTmed = [];
for n = 1:size(allWTnew)
    if allWTnew(n,3) < 4 || allWTnew(n,3) > 10
        WTlat = vertcat(WTlat,allWTnew(n,:));
    elseif allWTnew(n,3) > 4
        WTmed = vertcat(WTmed,allWTnew(n,:));
    end
end

% Box plot Linf/20a data
close all;
Linf_NumEx_WTlat = (WTlat(:,7))./WTlat(:,5);
Linf_NumEx_WTmed = (WTmed(:,7))./WTmed(:,5);
Linf_NumEx_LOFlat = (LOFlat(:,7))./LOFlat(:,5);
Linf_NumEx_LOFmed = (LOFmed(:,7))./LOFmed(:,5);

WTlatLabel = zeros(size(Linf_NumEx_WTlat));
WTmedLabel = zeros(size(Linf_NumEx_WTmed))+1;
LOFlatLabel = zeros(size(Linf_NumEx_LOFlat))+2;
LOFmedLabel = zeros(size(Linf_NumEx_LOFmed))+3;

%box chart script
xgroupData = vertcat(WTlatLabel,WTmedLabel,LOFlatLabel,LOFmedLabel);
yData = vertcat(Linf_NumEx_WTlat,Linf_NumEx_WTmed,Linf_NumEx_LOFlat,Linf_NumEx_LOFmed);

figure;
b = boxchart(xgroupData,yData);
b.BoxFaceColor = [0 0 0];
b.BoxFaceAlpha = 0;
%b.MarkerColor = 'k';
b.MarkerColor = [1 1 1];
hold on;

%cHere = [0 0.4470 0.7410];
scatter0 = [-.5 -.4 -.3 -.2 -.1 0 .1 .2 .3 .4 .5];
scatter1 = [-.6 -.5 -.4 -.3 -.2 -.1 0 .1 .2 .3 .4 .5 .6 .7];
scatter0 = scatter0./4;
scatter1 = scatter1./4;

scatter0 = linspace(-0.25,0.25,size(WTlatLabel,1));
scatter1 = linspace(0.75,1.25,size(WTmedLabel,1));
scatter2 = linspace(1.75,2.25,size(LOFlatLabel,1));
scatter3 = linspace(2.75,3.25,size(LOFmedLabel,1));

plot(scatter0,Linf_NumEx_WTlat,'.','color',[0 0 0],'Markersize',20); hold on;
plot(scatter1,Linf_NumEx_WTmed,'.','color',[0 0 0],'Markersize',20); hold on;
plot(scatter2,Linf_NumEx_LOFlat,'.','color',[0.6350 0.0780 0.1840],'Markersize',20); hold on;
plot(scatter3,Linf_NumEx_LOFmed,'.','color',[0.6350 0.0780 0.1840],'Markersize',20); hold on;
box on;

%Linf stats work
[p_WT_lat_med,h_WT_lat_med,stats_WT_lat_med] = ranksum(Linf_NumEx_WTlat,Linf_NumEx_WTmed); %p = 0.0695
[p_LOF_lat_med,h_LOF_lat_med,stats_LOF_lat_med] = ranksum(Linf_NumEx_LOFlat,Linf_NumEx_LOFmed); %p = 0.9143
[p_WTlat_LOFlat,h_WTlat_LOFlat,stats_WTlat_LOFlat] = ranksum(Linf_NumEx_WTlat,Linf_NumEx_LOFlat); %p = 0.0118
[p_WTmed_LOFmed,h_WTmed_LOFmed,stats_WTmed_LOFmed] = ranksum(Linf_NumEx_WTmed,Linf_NumEx_LOFmed); %p = 0.0176

plot([0,1],[3.15,3.15],'-','color','k','linewidth',1);
%txt1 = ['{\it p=}' num2str(p_WT_lat_med,'%.4f')];
txt1 = ['{\it p=}' num2str(round(p_WT_lat_med,2))];
text(0.2,3.25,txt1,'fontsize',16,'color',[0 0 0])

plot([2,3],[3.15,3.15],'-','color','k','linewidth',1);
%txt2 = ['{\it p=}' num2str(p_LOF_lat_med,'%.4f')];
txt2 = ['{\it p=}' num2str(round(p_LOF_lat_med,1))];
text(2.2,3.25,txt2,'fontsize',16,'color',[0 0 0])

plot([0,2],[3.45,3.45],'-','color','k','linewidth',1);
%txt3 = ['{\it p=}' num2str(p_WTlat_LOFlat,'%.4f')];
txt3 = ['{\it p=}' num2str(round(p_WTlat_LOFlat,2))];
text(0.5,3.55,txt3,'fontsize',16,'color','k')

plot([1,3],[3.75,3.75],'-','color','k','linewidth',1);
%txt3 = ['{\it p=}' num2str(p_WTmed_LOFmed,'%.4f')];
txt3 = ['{\it p=}' num2str(round(p_WTmed_LOFmed,2))];
text(1.5,3.85,txt3,'fontsize',16,'color','k')


ylabel({'Length_{reg} (\mum, Predicted) /';'# fgf20a Expressors'})
%xlim([-1,2]);
%ylim([-60,120]);
set(gca,'FontSize',18);

ax = gca;
ax.XTick = [0, 1, 2, 3];
ax.XTickLabels = {'WT Lat.','WT Med.','LOF Lat.','LOF Med.'};
xtickangle(60)
ylim([0,4]);

%% Supp. Fig. 7D

%load "fgf20a_timecourse.mat" from "SuppFig7"

close all;

data = exp2DataFull;

fig = figure;

fishSet = [1 2 3 4 11 12 13 14];
raySet = [2 3 6 7];

allData = {};


for p= 1:size(fishSet,2)
    fishHere = fishSet(p);
    for t = 1:size(raySet,2)
        rayHere = raySet(t);
        timeCollect = [];
        countCollect = [];  
        for q = 1:size(exp2DataFull,1)
            

            if data(q,2) == fishHere && data(q,3) == rayHere
                timeHere = data(q,4);
                countHere = data(q,5);
                timeCollect = horzcat(timeCollect,timeHere);
                countCollect = horzcat(countCollect,countHere);
            else 
                continue
            end


        end
        if ~isempty(countCollect)
            oneRayData = {fishHere,rayHere,timeCollect,countCollect};
            allData = vertcat(allData,oneRayData);
            oneRayData = {};
        else
            continue
        end
    end
end

% Plot GFP expression in WT and LOF over time
close all;
f = figure;
allCounts = [];
WTcounts = [];
LOFcounts = [];
for i = 1:size(allData)
    if cell2mat(allData(i,1)) < 10
        cl = [0 0 0];
    else
        cl = [1 0 0];
    end
    cl = horzcat(cl,0.3);
    plot(cell2mat(allData(i,3))./24,cell2mat(allData(i,4))./max(cell2mat(allData(i,4))),'.-','color',cl,'markersize',10,'linewidth',1); hold on;

    %collect data for averaging
    tempTime = cell2mat(allData(i,3))./24;
    tempCount = cell2mat(allData(i,4));
    if size(tempTime,2) == 7
        adjCount = tempCount;
    elseif size(tempTime,2) == 6
        adjCount = horzcat(tempCount,NaN);
    elseif size(tempTime,2) == 5
        adjCount = horzcat(tempCount,NaN,NaN);
    end
        allCounts = vertcat(allCounts,adjCount);
end

    %separate collected values into WT and LOF
    for qq = 1:size(allData)
        if cell2mat(allData(qq,1)) < 10
            WTcounts = vertcat(WTcounts,allCounts(qq,:)./max(allCounts(qq,:)));
        else
            LOFcounts = vertcat(LOFcounts,allCounts(qq,:)./max(allCounts(qq,:)));
        end
    end

    timeVals = [1 2 3 5 7 10 14];
    errorbar(timeVals,nanmean(WTcounts,1),nanstd(WTcounts,0,1)./sqrt(size(WTcounts,1)),'.-','markersize',15,'linewidth',3,'color',[0 0 0])
    errorbar(timeVals,nanmean(LOFcounts,1),nanstd(LOFcounts,0,1)./sqrt(size(LOFcounts,1)),'.-','markersize',15,'linewidth',3,'color',[1 0 0])
ylabel('Normalized #fgf20a Expressors')
xlabel('Time (days post amp.)')
yticks([0:0.2:1])
% ylim([0,1])
% xlim([0,15])
set(gca,'fontsize',16)

%% Supp. Fig. 7E

% Compare decay rate in WT and LOF fish
close all;
f = figure;
plot(timeWTnoAvg./24,meanWTnoAvg,'.','markersize',15,'color',[0 0 0 0.5]); hold on;
plot(timeLOF./24,meanLOF,'.','markersize',15,'color',[1 0 0 0.5]);

ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.StartPoint = [0.6 0.1 0.75];
opts.Upper = [1 1 1];

% Fit model to data.
[fitresultWT, gofWT] = fit(timeWTnoAvg./24, meanWTnoAvg, ft, opts );
[fitresultLOF, gofLOF] = fit(timeLOF./24, meanLOF, ft, opts );

xWT = [48:10:348];

h = plot(fitresultWT);
set(h,'color','k','linewidth',3);
h2 = plot(fitresultLOF);
set(h2,'color','r','linewidth',3);

ylabel('Average ERK Activity @ Amp. Site (AU)')
xlabel('Time (dpa)')
ylim([0.6,1.4])
%yticks([0:0.1:0.5])
set(gca,'fontsize',16)

legend_name={
    strcat("WT 1/b=",num2str(1/fitresultWT.b)),...
        strcat("LOF 1/b=",num2str(1/fitresultLOF.b))};
legend([h,h2], legend_name,...
        'Location','northeast','color','none','box','off');


