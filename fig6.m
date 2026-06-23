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

%load "fgf20a_exp1_2_allWTdata_Trim.mat" and "dataWT.mat" from "Fig6" folder


% fit is ALL (paper + new) data

xTem = linspace(0,5000,100);
yTem = linspace(0,5000,100);

figure; plot(allWTnew(:,6),allWTnew(:,5),'.','MarkerSize',20,'color','k'); hold on;
%plot(expData(:,6),expData(:,5),'.','MarkerSize',20,'color','m'); hold on;
plot(dataWT(:,5),dataWT(:,4),'.','MarkerSize',20,'color','k'); hold on;

%plot(xTem,yTem,'--','Color',[0.5 0.5 0.5],'LineWidth',2);

xDataHerePaper = vertcat(allWTnew(:,6));
yDataHerePaper = vertcat(allWTnew(:,5));

xDataHereNew = vertcat(dataWT(:,5));
yDataHereNew = vertcat(dataWT(:,4));

xDataHereAll = vertcat(allWTnew(:,6),dataWT(:,5));
yDataHereAll = vertcat(allWTnew(:,5),dataWT(:,4));

%fit combined data
ftA = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
excludedPoints = yDataHereAll > 6000;
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Exclude = excludedPoints;
[fitresultA, gofA] = fit( xDataHereAll, yDataHereAll, ftA, opts );

%q = plot([0,6000],[0,6000],'--','color',[0.7 0.7 0.7],'lineWidth',2);
%hP = plot((0:100:5000),fitresultP(0:100:5000));
%hN = plot((0:100:5000),fitresultN(0:100:5000));
hA = plot((1000:100:5000),fitresultA(1000:100:5000));
%fit_plot_l = plot((0:2:350)+48,fitresult_l(0:2:350));
%set(hP,'linestyle','-','linewidth',2,'color','k');
%set(hN,'linestyle','-','linewidth',2,'color','c');
set(hA,'linestyle','-','linewidth',2,'color','k');



ylabel('# Fgf20a Expressors');
xlabel('Length_{amp} (\mum)');
xlim([0 6000]);
ylim([0 6000]);

legend_name={strcat("R^2 = ",num2str(round(gofA.rsquare,4)))};%,...
    %strcat("y = x")};%,...
    %strcat("R^2 = ", num2str(round(gofN.rsquare,4)))
    %strcat("R^2 = ", num2str(round(gofA.rsquare,4)))};
%legend([hP,hN,hA], legend_name,...
legend([hA], legend_name,...
    'Location','best','color','none','box','off');
% legend([h,q], legend_name,...
%     'Location','best','color','none','box','off');
    

set(gca,'fontsize',16);

%"GFPplus_lamp.jpg"

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
% Plot lamp v. #exp - old & new LOF data
% fit is ALL (paper + new) WT data

%load "fgf20a_exp1_2_allWTdata_Trim.mat" and "dataWT.mat" from "Fig6" folder
%load "fgf20a_exp2_allLOFdata.mat" and "dataLOF.mat" from "Fig6" folder

figure; plot(allWTnew(:,6),allWTnew(:,5),'.','MarkerSize',20,'color','k'); hold on;
plot(dataWT(:,5),dataWT(:,4),'.','MarkerSize',20,'color','k','linewidth',2); hold on;
plot(exp2DataLOFNew(:,6),exp2DataLOFNew(:,5),'.','MarkerSize',20,'color','r'); hold on;
plot(dataLOF(:,5),dataLOF(:,4),'.','MarkerSize',20,'color','r','linewidth',2); hold on;

xDataHerePaper = vertcat(allWTnew(:,6));
yDataHerePaper = vertcat(allWTnew(:,5));
xDataHereNew = vertcat(dataWT(:,5));
yDataHereNew = vertcat(dataWT(:,4));
xDataHereAll = vertcat(allWTnew(:,6),dataWT(:,5));
yDataHereAll = vertcat(allWTnew(:,5),dataWT(:,4));

xDataHereLOF = vertcat(dataLOF(:,5),exp2DataLOFNew(:,6));
yDataHereLOF = vertcat(dataLOF(:,4),exp2DataLOFNew(:,5));

ftA = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
excludedPoints = yDataHereAll > 6000;
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Exclude = excludedPoints;
[fitresultA, gofA] = fit( xDataHereAll, yDataHereAll, ftA, opts );

% ftLOF = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% [fitresultLOF, gofLOF] = fit( xDataHereLOF, yDataHereLOF, ftLOF, opts );

ftLOF2 = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.Lower = [0.5945 124.1];
opts2.StartPoint = [0.7895 714.3];
opts2.Upper = [0.9845 1305];
excludedPoints = yDataHereLOF > 6000;
opts2.Exclude = excludedPoints;
[fitresultLOF2, gofLOF2] = fit( xDataHereLOF, yDataHereLOF, ftLOF2, opts2 );


hA = plot((1000:100:15000),fitresultA(1000:100:15000));
set(hA,'linestyle','-','linewidth',2,'color','k');
% hLOF2 = plot((1000:100:15000),fitresultLOF2(1000:100:15000));
% set(hLOF2,'linestyle','-','linewidth',2,'color','r');



ylabel('# Fgf20a Expressors');
%ylabel('Length_{reg} (\mum)');
xlabel('Length_{amp} (\mum)');
%xlim([0 6000]);
ylim([0 15000]);

legend_name={strcat("R^2 = ",num2str(round(gofA.rsquare,4)))};%,...
    %strcat("y = x")};%,...
    %strcat("R^2 = ", num2str(gof2.rsquare))
    %strcat("R^2 = ", num2str(gofLOF2.rsquare))};
legend([hA], legend_name,...
    'Location','best','color','none','box','off');
% legend([h,q], legend_name,...
%     'Location','best','color','none','box','off');
    

set(gca,'fontsize',16);

%"GFPplus_lamp_LOF_and_WT.jpg"


%% Fig. 6G

% Plot lreg v. #exp - new WT data (black) + new LOF data (red)
% fit is ALL (paper + new) data

%load "dataLOF2.mat" and "dataWT2.mat" from "Fig6" folder

%figure; plot(allWTnew(:,6),allWTnew(:,5),'.','MarkerSize',20,'color','k'); hold on;
%plot(expData(:,6),expData(:,5),'.','MarkerSize',20,'color','m'); hold on;
figure; plot(dataWT(:,4),dataWT(:,7),'.','MarkerSize',20,'color','k'); hold on;
plot(dataLOF(:,4),dataLOF(:,7),'.','MarkerSize',20,'color','r'); hold on;

%plot(xTem,yTem,'--','Color',[0.5 0.5 0.5],'LineWidth',2);



xDataHere2 = vertcat(dataWT(:,4));
yDataHere2 = vertcat(dataWT(:,7));

% idx2 = ~isnan(yDataHere2);
% xDataHere2 = xDataHere2(idx2);
% yDataHere2 = yDataHere2(idx2);

xDataHereLOF = vertcat(dataLOF(:,4));
yDataHereLOF = vertcat(dataLOF(:,7));

% idxLOF = ~isnan(yDataHereLOF);
% xDataHereLOF = xDataHereLOF(idxLOF);
% yDataHereLOF = yDataHereLOF(idxLOF);


%ft = fittype( 'a*x+0', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
excludedPoints = xDataHere2 > 6000;
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Exclude = excludedPoints;
[fitresult2, gof2] = fit( xDataHere2, yDataHere2, ft, opts );

ftLOF = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
[fitresultLOF, gofLOF] = fit( xDataHereLOF, yDataHereLOF, ftLOF, opts );
%[fitresult3, gof3] = fit( xDataHere3, yDataHere3, ft, opts );

excludedPoints = xDataHereLOF > 6000;




%q = plot([0,6000],[0,6000],'--','color',[0.7 0.7 0.7],'lineWidth',2);
%h1 = plot((0:100:5000),fitresult1(0:100:5000));
h2 = plot((1000:100:5000),fitresult2(1000:100:5000));
hLOF = plot((1000:100:6000),fitresultLOF(1000:100:6000));
%h3 = plot((0:100:5000),fitresult3(0:100:5000));
%fit_plot_l = plot((0:2:350)+48,fitresult_l(0:2:350));
%set(h1,'linestyle','-','linewidth',2,'color','k');
set(h2,'linestyle','-','linewidth',2,'color','k');
set(hLOF,'linestyle','-','linewidth',2,'color','r');



xlabel('# Fgf20a Expressors');
ylabel('Length_{reg} (\mum)');
%xlabel('Length_{amp} (\mum)');
xlim([0 6000]);
ylim([0 8000]);

legend_name={strcat("R^2 = ",num2str(round(gof2.rsquare,1)))%};%,...
    %strcat("y = x")};%,...
    %strcat("R^2 = ", num2str(gof2.rsquare))
    strcat("R^2 = ", num2str(round(gofLOF.rsquare,1)))};
legend([h2,hLOF], legend_name,...
    'Location','best','color','none','box','off');
% legend([h,q], legend_name,...
%     'Location','best','color','none','box','off');
    

set(gca,'fontsize',16);
