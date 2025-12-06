%% Fig. 5B

%load "analysis_mat_Aug_mod", "analysis_mat_juneERK", and "treatmentMat"
%from "Fig5"

% Merge matrices
analysis_mat_allERK = [analysis_mat_juneERK,analysis_mat_Aug_mod];

% Generate ERK matrix

fishSet = [1,2,3,4,5,6,7,8,11,12,13,14,15,16,18];
raySet = [2,3];

ERKmat = [];
for q = 1:size(fishSet,2)
    fishHere = fishSet(q);
    for k = 1:size(raySet,2)
        rayHere = raySet(k);
        
        avgKTRcol= [];
        semKTRcol = [];
        timeCollect = [];
        treatCollect = [];
        for i = 1:size(analysis_mat_allERK,2)
            if analysis_mat_allERK(i).fish == fishHere && analysis_mat_allERK(i).ray == rayHere
                ktrHere = analysis_mat_allERK(i).ktr-0.8;
                avgKTRhere = nanmean(ktrHere);
                stdKTRhere = nanstd(ktrHere);
                semKTRhere = stdKTRhere./sqrt(numel(ktrHere));
                timeHere = analysis_mat_allERK(i).hpa;
                treatmentHere = treatmentMat(i,2);
                avgKTRcol = horzcat(avgKTRcol,avgKTRhere);
                semKTRcol = horzcat(semKTRcol,semKTRhere);
                timeCollect = horzcat(timeCollect,timeHere);
                treatCollect = horzcat(treatCollect,treatmentHere);
            end
        end
        if ~isempty(avgKTRcol)
            ERKmatHere = horzcat(fishHere,rayHere,mean(treatCollect),avgKTRcol,semKTRcol,timeCollect(2)-timeCollect(1));
            ERKmat = vertcat(ERKmat,ERKmatHere);
        else
            continue
        end
    end
end

% subsetData

conData = [];
treatData = [];
for i = 1:size(ERKmat)
    if ERKmat(i,3) == 0
        conDataHere = ERKmat(i,:);
        conData = vertcat(conData,conDataHere);
    elseif ERKmat(i,3) == 1
        treatDataHere = ERKmat(i,:);
        treatData = vertcat(treatData,treatDataHere);
    end
end

% Plot new
figure;
for i = 1:size(ERKmat,1)
    if ERKmat(i,3)==0
        errorbar([0.1,0.2],[ERKmat(i,4),ERKmat(i,5)],[ERKmat(i,6),ERKmat(i,7)],'.-','linewidth',2,'markersize',10,'color',[0 0.4470 0.7410]); hold on;

    elseif ERKmat(i,3)==1
        errorbar([0.3,0.4],[ERKmat(i,4),ERKmat(i,5)],[ERKmat(i,6),ERKmat(i,7)],'.-','linewidth',2,'markersize',10,'color',[0 0.4470 0.7410]); hold on;

    end
end

hold on;

[hCon,pCon] = ttest2(conData(:,4),conData(:,5));
[hTreat,pTreat] = ttest2(treatData(:,4),treatDataHere(:,5));

plot([0.1,0.2],[mean(conData(:,4)),mean(conData(:,5))],'.-','linewidth',4,'markersize',20,'color','k'); hold on;
plot([0.3,0.4],[mean(treatData(:,4)),mean(treatData(:,5))],'.-','linewidth',4,'markersize',20,'color','k')

plot([0.1,0.2],[0.5,0.5],'-','linewidth',1,'color','k');
%t = text(0.1,0.525,['\it \fontsize{16} p=' num2str(pCon)]);
t = text(0.12,0.525,['\it \fontsize{16} p=' round(num2str(pCon,1))]);

plot([0.3,0.4],[0.5,0.5],'-','linewidth',1,'color','k');
t = text(0.31,0.525,['\it \fontsize{16} p=' round(num2str(pTreat,1))]);


ylabel('Average ERK Activity (AU)')
xlabel('Treatment');
ylim([0,0.55]);
xlim([0,0.5]);
set(gca,'fontsize',16);

ax = gca;
ax.XTick = [0.15, 0.35];
ax.XTickLabels = {'+DMSO','+BGJ'};


%% Fig. 5C

%load "analysis_mat_juneGEM.mat", "analysis_mat_Aug_mod.mat", and
%"treatMat.mat" from "Fig5" folder

% combine analysis mat.

analysis_mat_allGEM = [analysis_mat_juneGEM,analysis_mat_Aug_mod];

% Generate GEM matrix

fishSet = [1,2,3,4,5,6,7,8,11,12,13,14,15,16,18];
raySet = [2,3];

GEMmat = [];
for q = 1:size(fishSet,2)
    fishHere = fishSet(q);
    for k = 1:size(raySet,2)
        rayHere = raySet(k);
        
        noGEMcol= [];
        noCELLScol = [];
        timeCollect = [];
        treatCollect = [];
        for i = 1:size(analysis_mat_allGEM,2)
            if analysis_mat_allGEM(i).fish == fishHere && analysis_mat_allGEM(i).ray == rayHere
                gemHere = analysis_mat_allGEM(i).nucleiGEM;
                gemHere = sum(gemHere);
                cellsHere = analysis_mat_allGEM(i).nucleiALL;
                cellsHere = sum(cellsHere);
                timeHere = analysis_mat_allGEM(i).hpa;
                treatmentHere = treatmentMat(i,2);
                noGEMcol = horzcat(noGEMcol,gemHere);
                noCELLScol = horzcat(noCELLScol,cellsHere);
                timeCollect = horzcat(timeCollect,timeHere);
                treatCollect = horzcat(treatCollect,treatmentHere);
            end
        end
        if ~isempty(treatCollect)
            GEMmatHere = horzcat(fishHere,rayHere,mean(treatCollect),noGEMcol,noCELLScol,timeCollect(2)-timeCollect(1));
            GEMmat = vertcat(GEMmat,GEMmatHere);
        else
            continue
        end
    end
end

% subsetData

conData = [];
treatData = [];
for i = 1:size(GEMmat)
    if GEMmat(i,3) == 0
        conDataHere = GEMmat(i,:);
        conData = vertcat(conData,conDataHere);
    elseif GEMmat(i,3) == 1
        treatDataHere = GEMmat(i,:);
        treatData = vertcat(treatData,treatDataHere);
    end
end

% Plot new
figure;
for i = 1:size(GEMmat,1)
    if GEMmat(i,3)==0
        plot([0.1,0.2],[GEMmat(i,4),GEMmat(i,5)],'.-','linewidth',2,'markersize',10,'color',[0 0.4470 0.7410]); hold on;

    elseif GEMmat(i,3)==1
        plot([0.3,0.4],[GEMmat(i,4),GEMmat(i,5)],'.-','linewidth',2,'markersize',10,'color',[0 0.4470 0.7410]); hold on;

    end
end

[hCon,pCon] = ttest(conData(:,4),conData(:,5));
[hTreat,pTreat] = ttest(treatData(:,4),treatDataHere(:,5));

hold on;
plot([0.1,0.2],[mean(conData(:,4)),mean(conData(:,5))],'.-','linewidth',4,'markersize',30,'color','k');
plot([0.3,0.4],[mean(treatData(:,4)),mean(treatData(:,5))],'.-','linewidth',4,'markersize',30,'color','k');



plot([0.1,0.2],[280,280],'-','linewidth',1,'color','k');
%t = text(0.1,290,['\it \fontsize{16} p=' num2str(pCon)]);
t = text(0.12,290,['\it \fontsize{16} p=' round(num2str(pCon,1))]);

plot([0.3,0.4],[280,280],'-','linewidth',1,'color','k');
%t = text(0.3,290,['\it \fontsize{16} p=' num2str(pTreat)]);
t = text(0.31,290,['\it \fontsize{16} p=' round(num2str(pTreat,1))]);


ylabel('# Gem.+ Cells');
xlabel('Treatment');
%ylim([0,0.55]);
xlim([0,0.5]);
set(gca,'fontsize',16);

ax = gca;
ax.XTick = [0.15, 0.35];
ax.XTickLabels = {'+DMSO','+BGJ'};


%% Fig. 5D
%load "growthMat.mat" from "Fig5" folder

% Plot data - FIJI measurements
close all;


cHere = [0 0.4470 0.7410];
scatter0 = [-.5 -.4 -.3 -.2 -.1 0 .1 .2 .3 .4 .5 .6];
scatter1 = [-.4 -.3 -.2 -.1 0 .1 .2 .3 .4];
scatter0 = scatter0./4;
scatter1 = scatter1./4;

boxData = growthMat(:,[3,6]);
xgroupData = growthMat(:,3);
yData = growthMat(:,6);

%subsetData
conData = [];
treatData = [];
for i = 1:size(growthMat)
    if growthMat(i,3) == 0
        conDataHere = growthMat(i,6);
        conData = vertcat(conData,conDataHere);
    elseif growthMat(i,3) == 1
        treatDataHere = growthMat(i,6);
        treatData = vertcat(treatData,treatDataHere);
    end
end

%stats
[h,p,ci,stats] = ttest2(conData,treatData);

figure; 

b = boxchart(xgroupData,yData);
b.BoxFaceColor = [0 0 0];
b.BoxFaceAlpha = 0;
%b.MarkerColor = 'k';
b.MarkerColor = [1 1 1];
hold on;
plot(scatter0 + 0,conData,'.','color',cHere,'Markersize',20); hold on;
plot(scatter1 + 1,treatData,'.','color',cHere,'Markersize',20); hold on;
box on;

plot([0,1],[220,220],'-','linewidth',1,'color','k');
%t = text(0.2,230,['\it \fontsize{16} p=' num2str(p)]);
t = text(0.2,230,['\it \fontsize{16} p=' round(num2str(p,1))]);

ylabel('Growth (\mum)')
xlim([-1,2]);
%ylim([-60,120]);
set(gca,'FontSize',18);

ax = gca;
ax.XTick = [0, 1];
ax.XTickLabels = {'+DMSO','+BGJ'};

%% Fig. 5E

%load "tipData" and "ampData" from "Fig5" folder

tipTime = [0,3,7,10,14];
ampTime = [0,4,7,14,21];

% fgf10a

fgf10a_amp = cell2mat(ampData(1,2));
fgf10a_tip = cell2mat(tipData(1,2));

% fgf10b

fgf10b_amp = cell2mat(ampData(2,2));
fgf10b_tip = cell2mat(tipData(2,2));

% fgf3

fgf3_amp = cell2mat(ampData(22,2));
fgf3_tip = cell2mat(tipData(23,2));

% fgf20a

fgf20a_amp = cell2mat(ampData(17,2));
fgf20a_tip = cell2mat(tipData(18,2));

% fgf20b

fgf20b_amp = cell2mat(ampData(18,2));
fgf20b_tip = cell2mat(tipData(19,2));

% vertical bars - make bars wider, lines thicker
xyz = figure;
%fgf3
a = bar(0,fgf3_tip(1),0.3,'Facecolor',[0 0.4470 0.7410],'linewidth',2);hold on;
a.FaceAlpha = 1;
b = bar(0.3,fgf3_amp(1),0.3,'Facecolor',[0 0.4470 0.7410],'linewidth',2);hold on;
b.FaceAlpha = 0.5;
c = bar(0.7,fgf3_tip(2),0.3,'Facecolor',[0 0.4470 0.7410],'linewidth',2);hold on;
c.FaceAlpha = 1;
d = bar(1.0,fgf3_amp(2),0.3,'Facecolor',[0 0.4470 0.7410],'linewidth',2);hold on;
d.FaceAlpha = 0.5;
e = bar(1.4,fgf3_tip(3),0.3,'Facecolor',[0 0.4470 0.7410],'linewidth',2);hold on;
e.FaceAlpha = 1;
f = bar(1.7,fgf3_amp(3),0.3,'Facecolor',[0 0.4470 0.7410],'linewidth',2);hold on;
f.FaceAlpha = 0.5;
g = bar(2.1,fgf3_tip(5),0.3,'Facecolor',[0 0.4470 0.7410],'linewidth',2);hold on;
g.FaceAlpha = 1;
h = bar(2.4,fgf3_amp(4),0.3,'Facecolor',[0 0.4470 0.7410],'linewidth',2);hold on;
h.FaceAlpha = 0.5;

%fgf10a
aa = bar(3.0,fgf10a_tip(1),0.3,'Facecolor',[0.8500 0.3250 0.0980],'linewidth',2);hold on;
aa.FaceAlpha = 1;
bb = bar(3.3,fgf10a_amp(1),0.3,'Facecolor',[0.8500 0.3250 0.0980],'linewidth',2);hold on;
bb.FaceAlpha = 0.5;
cc = bar(3.7,fgf10a_tip(2),0.3,'Facecolor',[0.8500 0.3250 0.0980],'linewidth',2);hold on;
cc.FaceAlpha = 1;
dd = bar(4.0,fgf10a_amp(2),0.3,'Facecolor',[0.8500 0.3250 0.0980],'linewidth',2);hold on;
dd.FaceAlpha = 0.5;
ee = bar(4.4,fgf10a_tip(3),0.3,'Facecolor',[0.8500 0.3250 0.0980],'linewidth',2);hold on;
ee.FaceAlpha = 1;
ff = bar(4.7,fgf10a_amp(3),0.3,'Facecolor',[0.8500 0.3250 0.0980],'linewidth',2);hold on;
ff.FaceAlpha = 0.5;
gg = bar(5.1,fgf10a_tip(5),0.3,'Facecolor',[0.8500 0.3250 0.0980],'linewidth',2);hold on;
gg.FaceAlpha = 1;
hh = bar(5.4,fgf10a_amp(4),0.3,'Facecolor',[0.8500 0.3250 0.0980],'linewidth',2);hold on;
hh.FaceAlpha = 0.5;

%fgf10b
aaa = bar(6.0,fgf10b_tip(1),0.3,'Facecolor',[0.4940 0.1840 0.3560],'linewidth',2);hold on;
aaa.FaceAlpha = 1;
bbb = bar(6.3,fgf10b_amp(1),0.3,'Facecolor',[0.4940 0.1840 0.3560],'linewidth',2);hold on;
bbb.FaceAlpha = 0.5;
ccc = bar(6.7,fgf10b_tip(2),0.3,'Facecolor',[0.4940 0.1840 0.3560],'linewidth',2);hold on;
ccc.FaceAlpha = 1;
ddd = bar(7.0,fgf10b_amp(2),0.3,'Facecolor',[0.4940 0.1840 0.3560],'linewidth',2);hold on;
ddd.FaceAlpha = 0.5;
eee = bar(7.4,fgf10b_tip(3),0.3,'Facecolor',[0.4940 0.1840 0.3560],'linewidth',2);hold on;
eee.FaceAlpha = 1;
fff = bar(7.7,fgf10b_amp(3),0.3,'Facecolor',[0.4940 0.1840 0.3560],'linewidth',2);hold on;
fff.FaceAlpha = 0.5;
ggg = bar(8.1,fgf10b_tip(5),0.3,'Facecolor',[0.4940 0.1840 0.3560],'linewidth',2);hold on;
ggg.FaceAlpha = 1;
hhh = bar(8.4,fgf10b_amp(4),0.3,'Facecolor',[0.4940 0.1840 0.3560],'linewidth',2);hold on;
hhh.FaceAlpha = 0.5;

%fgf20a
aaaa = bar(9.0,fgf20a_tip(1),0.3,'Facecolor',[0.4660 0.6740 0.1880],'linewidth',2);hold on;
aaaa.FaceAlpha = 1;
bbbb = bar(9.3,fgf20a_amp(1),0.3,'Facecolor',[0.4660 0.6740 0.1880],'linewidth',2);hold on;
bbbb.FaceAlpha = 0.5;
cccc = bar(9.7,fgf20a_tip(2),0.3,'Facecolor',[0.4660 0.6740 0.1880],'linewidth',2);hold on;
cccc.FaceAlpha = 1;
dddd = bar(10.0,fgf20a_amp(2),0.3,'Facecolor',[0.4660 0.6740 0.1880],'linewidth',2);hold on;
dddd.FaceAlpha = 0.5;
eeee = bar(10.4,fgf20a_tip(3),0.3,'Facecolor',[0.4660 0.6740 0.1880],'linewidth',2);hold on;
eeee.FaceAlpha = 1;
ffff = bar(10.7,fgf20a_amp(3),0.3,'Facecolor',[0.4660 0.6740 0.1880],'linewidth',2);hold on;
ffff.FaceAlpha = 0.5;
gggg = bar(11.1,fgf20a_tip(5),0.3,'Facecolor',[0.4660 0.6740 0.1880],'linewidth',2);hold on;
gggg.FaceAlpha = 1;
hhhh = bar(11.4,fgf20a_amp(4),0.3,'Facecolor',[0.4660 0.6740 0.1880],'linewidth',2);hold on;
hhhh.FaceAlpha = 0.5;


plot([0.2,2.5],[25,25],'-','LineWidth',2,'color','k')
plot([3.2,5.5],[25,25],'-','LineWidth',2,'color','k')
plot([6.2,8.5],[25,25],'-','LineWidth',2,'color','k')
plot([9.2,11.5],[25,25],'-','LineWidth',2,'color','k')

txt3 = 'fgf3';
text(0.9,26,txt3,'Rotation',0,'FontSize',16)

txt10a = 'fgf10a';
text(3.7,26,txt10a,'Rotation',0,'FontSize',16)

txt10b = 'fgf10b';
text(6.7,26,txt10b,'Rotation',0,'FontSize',16)

txt20a = 'fgf20a';
text(9.7,26,txt20a,'Rotation',0,'FontSize',16)

xticks([0.15 0.85 1.55 2.25 3.15 3.85 4.55 5.25 6.15 6.85 7.55 8.25 9.15 9.85 10.55 11.25]);%	2.05000000000000	2.35000000000000	2.65000000000000	2.95000000000000	4.05000000000000	4.35000000000000	4.65000000000000	4.95000000000000	6.05000000000000	6.35000000000000	6.65000000000000	6.95000000000000]);
xticklabels({'0','3.5','7','14','0','3.5','7','14','0','3.5','7','14','0','3.5','7','14'});
%yticklabels({'14','7','3.5','0','14','7','3.5','0','14','7','3.5','0','14','7','3.5','0'})

% yticks([-8.95,-8.65,-8.35,-8.05,-6.95,-6.65,-6.35,-6.05,-4.95,-4.65,-4.35,-4.05,-2.95,-2.65,-2.35,-2.05,-0.95,-0.65,-0.35,-0.05])
% yticklabels({'14','7','3.5','0','14','7','3.5','0','14','7','3.5','0','14','7','3.5','0','14','7','3.5','0'})
ylabel('RNA Expression (normalized, TPM)')
xlabel({'Time (dpa)'})
xtickangle(90);
%ylabel({'Time (dpa)',''})
set(gca,'fontsize',16)

ylim([0,27])


%% Fig. 5G

%load "analysis_mat_timeaverage_042624.mat" from "Fig5" folder
%load "xdataSim.mat" and "yDataSim.mat" from "Fig5" folder

%step 1
% aggregate single nuclei first - subtract later - bin data into 5 groups
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
set(gca,'fontsize',16)
xlim([0,1]);



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


% Plot Alvin's simulations with my spatial profile

tLevels = 20;
color_map_red=colormap(cbrewer2('Purples',tLevels));
b = [18,16,14,12,10,8,6,4];

g3= figure;
for q = 1:size(XData,1)
    plot(XData(q,:),YData(q,:),'-','linewidth',4,'color',color_map_red(b(q),:)); hold on;
end


% Combine Alvin's simulations with my spatial profile - 2

tLevels = 20;
color_map_red=colormap(cbrewer2('Purples',tLevels));
b = [18,16,14,12,10,8,6,4];

g3= figure;
for q = 1:size(XData,1)
    plot(XData(q,:),YData(q,:),'-','linewidth',4,'color',color_map_red(b(q),:)); hold on;
end

tLevels = 20;
color_map_red=colormap(cbrewer2('Purples',tLevels));

ccc = plot(medPhiCol2-linspace(0,medPhiCol2,10),nanmean(bin_y_2phi,2),'o','color',color_map_red(4,:),'markersize',5,'linewidth',2);
ddd = plot(medPhiCol3-linspace(0,medPhiCol3,10),nanmean(bin_y_3phi,2),'o','color',color_map_red(6,:),'markersize',5,'linewidth',2);
eee = plot(medPhiCol4-linspace(0,medPhiCol4,10),nanmean(bin_y_4phi,2),'o','color',color_map_red(8,:),'markersize',5,'linewidth',2);
fff = plot(medPhiCol5-linspace(0,medPhiCol5,10),nanmean(bin_y_5phi,2),'o','color',color_map_red(10,:),'markersize',5,'linewidth',2);
ggg = plot(medPhiCol6-linspace(0,medPhiCol6,10),nanmean(bin_y_6phi,2),'o','color',color_map_red(12,:),'markersize',5,'linewidth',2);
hhh = plot(medPhiCol7-linspace(0,medPhiCol7,10),nanmean(bin_y_7phi,2),'o','color',color_map_red(14,:),'markersize',5,'linewidth',2);
iii = plot(medPhiCol8-linspace(0,medPhiCol8,10),nanmean(bin_y_8phi,2),'o','color',color_map_red(16,:),'markersize',5,'linewidth',2);
jjj = plot(medPhiCol9-linspace(0,medPhiCol9,10),nanmean(bin_y_9phi,2),'o','color',color_map_red(18,:),'markersize',5,'linewidth',2);

hold off
 
ylim([0,0.5])
xlabel('Average x/Length Amputated')
ylabel('ERK Activity')
yticks([0 0.1 0.2 0.3 0.4 0.5])
set(gca,'fontsize',16)
xlim([0,1]);

legend_name={
    
        strcat("phi < 0.2 (25)"),...
        strcat("phi < 0.3 (15)"),...
        strcat("phi < 0.4 (4)"),...
        strcat("phi < 0.5 (22)"),...
        strcat("phi < 0.6 (38)"),...
        strcat("phi < 0.7 (22)"),...
        strcat("phi < 0.8 (24)"),...
        strcat("phi < 0.9 (15)"),...
        };
legend([ccc,ddd,eee,fff,ggg,hhh,iii,jjj], legend_name,...
        'Location','northeast','color','none','box','off');

%% Fig. 5h

%% A(phi) color by Lamp - 3 groups

E0 = 0.8;

mean_ktr = arrayfun(@(x)mean(x.ktr(x.trim_logical)),analysis_mat);
mean_ktr_bin = arrayfun(@(x)nanmean(x.averageKTR),analysis_mat);
averageKTR_u_cell = arrayfun(@(s)...
    bin_average(s.ccrot(s.trim_logical,1),...
    adjust_neg(s.ktr(s.trim_logical)),10,0,s.L_reg),...
    analysis_mat,'UniformOutput',false);

mean_ktr_bin10_offsetE0 = cellfun(@nanmean,averageKTR_u_cell);

L_reg = [analysis_mat.L_reg];
%L_amp = [analysis_mat.L_amp];
%phi = L_reg./[L_amp];
ray = [analysis_mat.ray];
fish = [analysis_mat.fish];
time = [analysis_mat.hpa];

ydata = mean_ktr-E0;

%% generate treatment matrix

treat = [];
for i = 1:size(fish,2)
    if fish(i) < 3 || fish(i) > 8
        treatHere = 0;
    else
        treatHere = 1;
    end
    
    treat = vertcat(treat,treatHere);
end

%% collect data of interest

dataHere = horzcat(fish',time',ydata',treat);

%%
plotData = [];
for fish = 1:10

    for i = 1:size(dataHere,1)
        if dataHere(i,1) == fish
            treatVal = dataHere(i,4);
            if dataHere(i,2) == 120
                preERK = dataHere(i,3);
            elseif dataHere(i,2) == 132
                postERK = dataHere(i,3);
            end

        end

    end
    ERKhere = horzcat(preERK,postERK);
    plotDataHere = horzcat(fish,treatVal,ERKhere);
    plotData = vertcat(plotData,plotDataHere);

end

%%

figure;

conData = [];
treatData = [];
for i = 1:size(plotData)
    if plotData(i,2) == 0
        plot([0,12],[plotData(i,3),plotData(i,4)],'.-','color','k','linewidth',2); hold on;
        conData = vertcat(conData,[plotData(i,3),plotData(i,4)]);
    elseif plotData(i,2) == 1
        plot([24,36],[plotData(i,3),plotData(i,4)],'.-','color','k','linewidth',2); hold on;
        treatData = vertcat(treatData,[plotData(i,3),plotData(i,4)]);
    end
end

%% plot average ERK - compare Pre - compare Post
f=figure;

meanCon = nanmean(conData);
stdCon = nanstd(conData);
semCon = stdCon./sqrt(size(conData,1));

meanTreat = nanmean(treatData);
stdTreat = nanstd(treatData);
semTreat = stdTreat./sqrt(size(treatData,1));

[h,p,ci,stats] = ttest2(conData(:,1),treatData(:,1));
[h2,p2,ci2,stats2] = ttest2(conData(:,2),treatData(:,2));

xgroupData = [1 1 1 1 2 2 2 2 2 2 3 3 3 3 4 4 4 4 4 4];
yData = vertcat(conData(:,1),treatData(:,1),conData(:,2),treatData(:,2));

b = boxchart(xgroupData,yData); hold on;
b.BoxFaceColor = [0 0 0];
b.BoxFaceAlpha = 0;
%b.MarkerColor = 'k';
b.MarkerColor = [1 1 1];
hold on;
box on;

plot(linspace(0.8,1.2,4),[conData(:,1)],'.','color',[0 0.4470 0.7410],'markersize',20); hold on;
plot(linspace(1.8,2.2,6),[treatData(:,1)],'.','color',[0 0.4470 0.7410],'markersize',20); hold on;
plot(linspace(2.8,3.2,4),[conData(:,2)],'.','color',[0 0.4470 0.7410],'markersize',20); hold on;
plot(linspace(3.8,4.2,6),[treatData(:,2)],'.','color',[0 0.4470 0.7410],'markersize',20); hold on;

plot([1,2],[0.43,0.43],'-','linewidth',1,'color','k');
%t = text(1.2,0.45,['\it \fontsize{16} p=' num2str(p)]);
t = text(1.3,0.45,['\it \fontsize{16} p=' round(num2str(p,1))]);
plot([3,4],[0.43,0.43],'-','linewidth',1,'color','k');
%t = text(3.2,0.45,['\it \fontsize{16} p=' num2str(p2)]);
t = text(3.3,0.45,['\it \fontsize{16} p=' round(num2str(p2,1))]);

ylim([0,0.5]);
%xlim([-1,4]);

xticks([1 2 3 4])
xticklabels({'Pre(+Water)','Pre(+Cyclohex.)','Post(+Water)','Post(+Cyclohex.)'})
ylabel('Average ERK Activity (AU)');
set(gca,'fontsize',16)

%% Fig. 5H
% plot average ERK - compare Pre - compare Post

%load "analysis_mat_cyclo.mat" from "Fig5" folder

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
%L_amp = [analysis_mat.L_amp];
%phi = L_reg./[L_amp];
ray = [analysis_mat.ray];
fish = [analysis_mat.fish];
time = [analysis_mat.hpa];

ydata = mean_ktr-E0;

% generate treatment matrix

treat = [];
for i = 1:size(fish,2)
    if fish(i) < 3 || fish(i) > 8
        treatHere = 0;
    else
        treatHere = 1;
    end
    
    treat = vertcat(treat,treatHere);
end

% collect data of interest

dataHere = horzcat(fish',time',ydata',treat);

%
plotData = [];
for fish = 1:10

    for i = 1:size(dataHere,1)
        if dataHere(i,1) == fish
            treatVal = dataHere(i,4);
            if dataHere(i,2) == 120
                preERK = dataHere(i,3);
            elseif dataHere(i,2) == 132
                postERK = dataHere(i,3);
            end

        end

    end
    ERKhere = horzcat(preERK,postERK);
    plotDataHere = horzcat(fish,treatVal,ERKhere);
    plotData = vertcat(plotData,plotDataHere);

end


%subsetData
conData = [];
treatData = [];
for i = 1:size(plotData)
    if plotData(i,2) == 0
        plot([0,12],[plotData(i,3),plotData(i,4)],'.-','color','k','linewidth',2); hold on;
        conData = vertcat(conData,[plotData(i,3),plotData(i,4)]);
    elseif plotData(i,2) == 1
        plot([24,36],[plotData(i,3),plotData(i,4)],'.-','color','k','linewidth',2); hold on;
        treatData = vertcat(treatData,[plotData(i,3),plotData(i,4)]);
    end
end

%figure
f=figure;

meanCon = nanmean(conData);
stdCon = nanstd(conData);
semCon = stdCon./sqrt(size(conData,1));

meanTreat = nanmean(treatData);
stdTreat = nanstd(treatData);
semTreat = stdTreat./sqrt(size(treatData,1));

[h,p,ci,stats] = ttest2(conData(:,1),treatData(:,1));
[h2,p2,ci2,stats2] = ttest2(conData(:,2),treatData(:,2));

xgroupData = [1 1 1 1 2 2 2 2 2 2 3 3 3 3 4 4 4 4 4 4];
yData = vertcat(conData(:,1),treatData(:,1),conData(:,2),treatData(:,2));

b = boxchart(xgroupData,yData); hold on;
b.BoxFaceColor = [0 0 0];
b.BoxFaceAlpha = 0;
%b.MarkerColor = 'k';
b.MarkerColor = [1 1 1];
hold on;
box on;

plot(linspace(0.8,1.2,4),[conData(:,1)],'.','color',[0 0.4470 0.7410],'markersize',20); hold on;
plot(linspace(1.8,2.2,6),[treatData(:,1)],'.','color',[0 0.4470 0.7410],'markersize',20); hold on;
plot(linspace(2.8,3.2,4),[conData(:,2)],'.','color',[0 0.4470 0.7410],'markersize',20); hold on;
plot(linspace(3.8,4.2,6),[treatData(:,2)],'.','color',[0 0.4470 0.7410],'markersize',20); hold on;

plot([1,2],[0.43,0.43],'-','linewidth',1,'color','k');
%t = text(1.2,0.45,['\it \fontsize{16} p=' num2str(p)]);
t = text(1.3,0.45,['\it \fontsize{16} p=' round(num2str(p,1))]);
plot([3,4],[0.43,0.43],'-','linewidth',1,'color','k');
%t = text(3.2,0.45,['\it \fontsize{16} p=' num2str(p2)]);
t = text(3.3,0.45,['\it \fontsize{16} p=' round(num2str(p2,1))]);

ylim([0,0.5]);
%xlim([-1,4]);

xticks([1 2 3 4])
xticklabels({'Pre(+Water)','Pre(+Cyclohex.)','Post(+Water)','Post(+Cyclohex.)'})
ylabel('Average ERK Activity (AU)');
set(gca,'fontsize',16)

