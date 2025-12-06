%% Fig. 6C

%load "fgf20a_exp2Full_data.mat", "fgf20aConcatData.mat", and
%"fgf20aMeanConcatData.mat"

% Collect WT GFP expression data

WTgfp = [];
for i = 1:size(exp2DataFull)
    if exp2DataFull(i,2) < 10
        dataHere = exp2DataFull(i,:);
        WTgfp = vertcat(WTgfp,dataHere);
    else
        continue
    end
end

% Collect GFP counts from WT timecourse data
fishUnique = unique(WTgfp(:,2));
rayUnique = unique(WTgfp(:,3));

GFPplotData = {};

    for fishHere = 1:size(fishUnique)
        fish = fishUnique(fishHere);
        singleRay = {};
        for rayHere = 1:size(rayUnique)
            ray = rayUnique(rayHere);
            
            timeCollect = [];
            countCollect = [];
            for q = 1:size(WTgfp)

                if fish == WTgfp(q,2) && ray == WTgfp(q,3)
                    timeTemp = WTgfp(q,4);
                    countTemp = WTgfp(q,5);
                    timeCollect = horzcat(timeCollect,timeTemp);
                    countCollect = horzcat(countCollect,countTemp);
                %else
                %    a=1
                end
            end
            
            if ~isempty(timeCollect)
                singleRay = [{fish},{ray},{timeCollect},{countCollect}];
                GFPplotData = vertcat(GFPplotData,singleRay);
            else
                continue
            end
        end
    end

fgfSeq = concatData;
fgfSeqMean = meanConcatData;

% Make figure - include average of expression data - flip green and black

fgfSeqStd = std(fgfSeq);
fgfSeqSem = fgfSeqStd./sqrt(3);

fig = figure;
right_color = [0 0 0];
left_color = [0.4660 0.6740 0.1880];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

colGFP = [];
for cc = 1:size(GFPplotData)
    GFPhere = cell2mat(GFPplotData(cc,4))./max(cell2mat(GFPplotData(cc,4)));
    if size(GFPhere,2) == 3
        GFPhere = horzcat(GFPhere,NaN,NaN,NaN,NaN);
    elseif size(GFPhere,2) == 6
        GFPhere = horzcat(GFPhere,NaN);
%     else
%         continue
    end
    colGFP = vertcat(colGFP,GFPhere);
end

yyaxis left
seqTime = [0 3 7 10 14];
for g = 1:size(GFPplotData)
    plot(cell2mat(GFPplotData(g,3))./24,cell2mat(GFPplotData(g,4))./max(cell2mat(GFPplotData(g,4))),'.-','color',[0.4660 0.6740 0.1880 0.3],'linewidth',1,'markersize',10); hold on;
end
timeVal = [24	48	72	120	168	240	336];
semGFP = nanstd(colGFP)./sqrt(numel(colGFP));
errorbar(timeVal./24,nanmean(colGFP),semGFP,'.-','color',[0.4660 0.6740 0.1880 1],'linewidth',3,'markersize',30); hold on;


xlabel('Time (days post amp.)');
ylabel('Normalized #fgf20a Expressors');
ylim([0,1.05])


yyaxis right

maxFgfSeq = max(fgfSeq,[],2);
fgfSeqNorm = fgfSeq./maxFgfSeq;

fgfSeqStdNorm = std(fgfSeqNorm);
fgfSeqSemNorm = fgfSeqStdNorm./sqrt(3);

plot(seqTime,fgfSeqNorm,'.','MarkerSize',35,'color',[0 0 0])
%plot(seqTime,fgfSeqNorm,'.','MarkerSize',25,'color',[0 0 0])
errorbar(seqTime,mean(fgfSeqNorm),fgfSeqSemNorm,'-','MarkerSize',30,'color',[0 0 0],'lineWidth',4)
%errorbar(seqTime,mean(fgfSeqNorm),fgfSeqSemNorm,'-','MarkerSize',30,'color',[0 0 0],'lineWidth',2)



ylabel('Normalized fgf20a Expression (TPM)');

xlim([-1,15])
ylim([0,1.05])

set(gca,'fontsize',16)

%% Fig. 6D
% This code generates # expressors v. Lamp for WT data set

%load "fgf20a_exp1_2_allWTdata_Trim.mat" from "Fig6" folder

figure; plot(allWTnew(:,6),allWTnew(:,5),'.','MarkerSize',20,'color','k'); hold on;

xDataHere = allWTnew(:,6);
yDataHere = allWTnew(:,5);

%ft = fittype( 'a*x+0', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
[fitresult, gof] = fit( xDataHere, yDataHere, ft, opts );

%q = plot([0,6000],[0,6000],'--','color',[0.7 0.7 0.7],'lineWidth',2);
h = plot((0:100:5000),fitresult(0:100:5000));
%fit_plot_l = plot((0:2:350)+48,fitresult_l(0:2:350));
set(h,'linestyle','-','linewidth',2,'color','k');



ylabel('# Fgf20a Expressors');
xlabel('Length_{amp} (\mum)');
xlim([0 6000]);
ylim([0 6000]);

legend_name={strcat("R^2 = ",num2str(round(gof.rsquare,1)))};%,...
    %strcat("y = x")};%,...
    %strcat("R^2 = ", num2str(gof_s.rsquare))};
legend([h], legend_name,...
    'Location','best','color','none','box','off');
% legend([h,q], legend_name,...
%     'Location','best','color','none','box','off');
    

set(gca,'fontsize',16);

%% Fig. 6E

%load "cellCollect.mat" from "Fig6" folder

% Plot length GFP+ / Lamp v. time (grays)

tLevels = 10;
color_map=colormap(cbrewer2('Greys',tLevels));
lampMinHere = 706;
lampMaxHere = 4053;

figure;
for i = 1:size(cellCollect)
    timeHere = cell2mat(cellCollect(i,3));
    lGFPhere = cell2mat(cellCollect(i,4));
    lAmpHere = cell2mat(cellCollect(i,5));
    
    [sortTime,sortIdx] = sort(timeHere);
    sortLength = lGFPhere(sortIdx);
    sortAmp = lAmpHere(sortIdx);
    
    
    %colorCode = ((sortAmp(1)-700)./3353).*100+1;
    %plot(sortTime(3:6),sortLength(3:6)./sortAmp(3:6),'.-','markersize',15,'linewidth',2,'color',color_map(i+2,:)); hold on;
    plot(sortTime./24,sortLength./sortAmp,'.-','markersize',15,'linewidth',2,'color',color_map(i+2,:)); hold on;
end

xlabel('Time (days post amp.)');
ylabel({['Length Fgf20a Exp. Domain / Length Amp.']})
%ylabel('Length Fgf20a Expression Domain (\mum) / Length Amputated (\mum)');
set(gca,'fontsize',16);
ylim([0,0.3])
xlim([0,12])

%% Fig. 6F

% load "fgf20a_exp1_2_allWTdata_Trim.mat" and "fgf20a_exp2_allLOFdata.mat"
% from "Fig6" folder

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

% Box plot 20a/Lamp data

Lamp_NumEx_WTlat = (WTlat(:,5))./WTlat(:,6);
Lamp_NumEx_WTmed = (WTmed(:,5))./WTmed(:,6);
Lamp_NumEx_LOFlat = (LOFlat(:,5))./LOFlat(:,6);
Lamp_NumEx_LOFmed = (LOFmed(:,5))./LOFmed(:,6);

WTlatLabel = zeros(size(Lamp_NumEx_WTlat));
WTmedLabel = zeros(size(Lamp_NumEx_WTmed))+1;
LOFlatLabel = zeros(size(Lamp_NumEx_LOFlat))+2;
LOFmedLabel = zeros(size(Lamp_NumEx_LOFmed))+3;

%box chart script
xgroupData2 = vertcat(WTlatLabel,WTmedLabel,LOFlatLabel,LOFmedLabel);
yData2 = vertcat(Lamp_NumEx_WTlat,Lamp_NumEx_WTmed,Lamp_NumEx_LOFlat,Lamp_NumEx_LOFmed);

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

plot(scatter0,Lamp_NumEx_WTlat,'.','color',[0 0 0],'Markersize',20); hold on;
plot(scatter1,Lamp_NumEx_WTmed,'.','color',[0 0 0],'Markersize',20); hold on;
plot(scatter2,Lamp_NumEx_LOFlat,'.','color',[0.6350 0.0780 0.1840],'Markersize',20); hold on;
plot(scatter3,Lamp_NumEx_LOFmed,'.','color',[0.6350 0.0780 0.1840],'Markersize',20); hold on;
box on;





%Linf stats work
[p2_WT_lat_med,h2_WT_lat_med,stats2_WT_lat_med] = ranksum(Lamp_NumEx_WTlat,Lamp_NumEx_WTmed); %p = 0.1255
[p2_LOF_lat_med,h2_LOF_lat_med,stats2_LOF_lat_med] = ranksum(Lamp_NumEx_LOFlat,Lamp_NumEx_LOFmed); %p = 0.0095
[p2_WTlat_LOFlat,h2_WTlat_LOFlat,stats2_WTlat_LOFlat] = ranksum(Lamp_NumEx_WTlat,Lamp_NumEx_LOFlat); %p = 0.0022
[p2_WTmed_LOFmed,h2_WTmed_LOFmed,stats2_WTmed_LOFmed] = ranksum(Lamp_NumEx_WTmed,Lamp_NumEx_LOFmed); %p = 0.3277

plot([0,1],[2.3,2.3],'-','color','k','linewidth',1);
%txt1 = ['{\it p=}' num2str(p2_WT_lat_med,'%.4f')];
txt1 = ['{\it p=}' num2str(round(p2_WT_lat_med,1))];
text(0.2,2.4,txt1,'fontsize',16,'color',[0 0 0])

plot([2,3],[2.3,2.3],'-','color','k','linewidth',1);
%txt2 = ['{\it p=}' num2str(p2_LOF_lat_med,'%.4f')];
txt2 = ['{\it p=}' num2str(round(p2_LOF_lat_med,3))];
text(2.2,2.4,txt2,'fontsize',16,'color',[0 0 0])

plot([0,2],[0.3,0.3],'-','color','k','linewidth',1);
%txt3 = ['{\it p=}' num2str(p2_WTlat_LOFlat,'%.4f')];
txt3 = ['{\it p=}' num2str(round(p2_WTlat_LOFlat,3))];
text(0.5,0.24,txt3,'fontsize',16,'color','k')

plot([1,3],[0.15,0.15],'-','color','k','linewidth',1);
%txt3 = ['{\it p=}' num2str(p2_WTmed_LOFmed,'%.4f')];
txt3 = ['{\it p=}' num2str(round(p2_WTmed_LOFmed,1))];
text(1.5,0.09,txt3,'fontsize',16,'color','k')



ylabel({'# fgf20a Expressors';'/ Length_{amp} (\mum)'})
%xlim([-1,2]);
ylim([0,2.5]);
set(gca,'FontSize',18);

ax = gca;
ax.XTick = [0, 1, 2, 3];
ax.XTickLabels = {'WT Lat.','WT Med.','LOF Lat.','LOF Med.'};
xtickangle(60)

%% Fig. 6G

% load "fgf20a_exp1_2_allWTdata_Trim.mat" and "fgf20a_exp2_allLOFdata.mat"
% from "Fig6" folder

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

% Box plot 20a/Lreg & 20a/Linf data for LOF
close all;
Linf_NumEx_WTlat = (WTlat(:,5))./WTlat(:,7);
Linf_NumEx_WTmed = (WTmed(:,5))./WTmed(:,7);
Linf_NumEx_LOFlat = (LOFlat(:,5))./LOFlat(:,7);
Linf_NumEx_LOFmed = (LOFmed(:,5))./LOFmed(:,7);

Lamp_NumEx_WTlat = (WTlat(:,5))./WTlat(:,6);
Lamp_NumEx_WTmed = (WTmed(:,5))./WTmed(:,6);
Lamp_NumEx_LOFlat = (LOFlat(:,5))./LOFlat(:,6);
Lamp_NumEx_LOFmed = (LOFmed(:,5))./LOFmed(:,6);

%WTlatLabel = zeros(size(Linf_NumEx_WTlat));
%WTmedLabel = zeros(size(Linf_NumEx_WTmed))+1;
LOFlatLabel1 = zeros(size(Linf_NumEx_LOFlat));
LOFmedLabel1 = zeros(size(Linf_NumEx_LOFmed))+1;
LOFlatLabel2 = zeros(size(Linf_NumEx_LOFlat))+2;
LOFmedLabel2 = zeros(size(Linf_NumEx_LOFmed))+3;

%box chart script
xgroupData = vertcat(LOFlatLabel1,LOFmedLabel1,LOFlatLabel2,LOFmedLabel2);
%xgroupData = vertcat(WTlatLabel,WTmedLabel,LOFlatLabel,LOFmedLabel);
%yData = vertcat(Linf_NumEx_WTlat,Linf_NumEx_WTmed,Linf_NumEx_LOFlat,Linf_NumEx_LOFmed);
yData = vertcat(Lamp_NumEx_LOFlat,Lamp_NumEx_LOFmed,Linf_NumEx_LOFlat,Linf_NumEx_LOFmed);

q = figure;

left_color = [0.6350 0.0780 0.1840];
right_color = [0.4940 0.1840 0.5560]; %purple
right_color = [0.3010 0.7450 0.9330]; %blue
right_color = [0 0.4470 0.7410]; %darker blue
set(q,'defaultAxesColorOrder',[left_color; right_color]);


yyaxis left

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

scatter0 = linspace(-0.25,0.25,size(LOFlatLabel1,1));
scatter1 = linspace(0.75,1.25,size(LOFmedLabel1,1));
scatter2 = linspace(1.75,2.25,size(LOFlatLabel2,1));
scatter3 = linspace(2.75,3.25,size(LOFmedLabel2,1));

plot(scatter0,Lamp_NumEx_LOFlat,'.','color',[0.6350 0.0780 0.1840],'Markersize',20); hold on;
plot(scatter1,Lamp_NumEx_LOFmed,'.','color',[0.6350 0.0780 0.1840],'Markersize',20); hold on;
plot(scatter2,Linf_NumEx_LOFlat,'.','color',right_color,'Markersize',20); hold on;
plot(scatter3,Linf_NumEx_LOFmed,'.','color',right_color,'Markersize',20); hold on;
box on;

%Linf stats work
[p_LOF_Lamp,h_LOF_Lamp,stats_LOF_Lamp] = ranksum(Lamp_NumEx_LOFlat,Lamp_NumEx_LOFmed); %p = 0.0695
[p_LOF_Linf,h_LOF_Linf,stats_LOF_Linf] = ranksum(Linf_NumEx_LOFlat,Linf_NumEx_LOFmed); %p = 0.9143
%[p_WTlat_LOFlat,h_WTlat_LOFlat,stats_WTlat_LOFlat] = ranksum(Linf_NumEx_WTlat,Linf_NumEx_LOFlat); %p = 0.0118
%[p_WTmed_LOFmed,h_WTmed_LOFmed,stats_WTmed_LOFmed] = ranksum(Linf_NumEx_WTmed,Linf_NumEx_LOFmed); %p = 0.0176

plot([0,1],[2.3,2.3],'-','color','k','linewidth',1);
%txt1 = ['{\it p=}' num2str(p_LOF_Lamp,'%.4f')];
txt1 = ['{\it p=}' num2str(round(p_LOF_Lamp,3))];
text(0.2,2.4,txt1,'fontsize',16,'color',[0 0 0])

plot([2,3],[2.3,2.3],'-','color','k','linewidth',1);
txt2 = ['{\it p=}' num2str(round(p_LOF_Linf,1))];
text(2.2,2.4,txt2,'fontsize',16,'color',[0 0 0])

% plot([0,2],[3.45,3.45],'-','color','k','linewidth',1);
% txt3 = ['{\it p=}' num2str(p_WTlat_LOFlat,'%.4f')];
% text(0.5,3.55,txt3,'fontsize',16,'color','k')
% 
% plot([1,3],[3.75,3.75],'-','color','k','linewidth',1);
% txt3 = ['{\it p=}' num2str(p_WTmed_LOFmed,'%.4f')];
% text(1.5,3.85,txt3,'fontsize',16,'color','k')

set(q,'defaultAxesColorOrder',[left_color; right_color]);
ylabel({'# fgf20a Expressors /' ; 'Length_{amp}'})
ylim([0,2.5])

yyaxis right
ylabel({'# fgf20a Expressors' ; 'Length_{reg} (\mum, Predicted)'})
%xlim([-1,2]);
%ylim([-60,120]);
set(gca,'FontSize',18);
set(q,'defaultAxesColorOrder',[left_color; right_color]);

ax = gca;
ax.XTick = [0, 1, 2, 3];
ax.XTickLabels = {'LOF Lat.','LOF Med.','LOF Lat.','LOF Med.'};
xtickangle(60)
ylim([0,2.5]);