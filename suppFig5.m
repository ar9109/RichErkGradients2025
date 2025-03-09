%% Supp. Fig. 5B & C

% read matrix
%Excel file is found in "SuppFig5" folder. 
%Please update the path information (line 7) for importing the matrix

LOFdata = readmatrix('/Volumes/AshleyData3/31Oct22_LongfinTimecourse/31Oct22_LongfinTimecourse.xlsx');

% Collect data

% for i = 1:size(LOFdata)
%     selectData = LOFdata(i,[1,2,5,9,12,15,18,21]);
% end

selectData = LOFdata(:,[1,2,5,9,12,15,18,21,24,27,30,33,36,39,42,45]);

selectDataAdjust = selectData(:,3:12).*(96000/158);

selectData = horzcat(selectData(:,1:2),selectDataAdjust);

timePoints = [1,2,4,7,15,21,28,35,42,49,56,64,70];

%Generate Supp. Fig. 5B
% Generate WT timecourse plot

% for i = 1:size(LOFdata)
%     selectData = LOFdata(i,[1,2,5,9,12,15,18,21]);
% end



selectData = LOFdata(:,[1,2,5,9,12,15,18,21,24,27,30,33,36,39, 42, 45]);

selectDataAdjust = selectData(:,4:16)./selectData(:,3);

LampMeasurements = selectData(:,3).*(96000/158);

selectData = horzcat(selectData(:,1:2),selectDataAdjust);

timePoints = [1,2,4,7,15,21,28,35,42,49,56, 64, 70];

LampMeasurementsWT = LampMeasurements(1:36,:);

% tLevels = 5;
% color_map=colormap(cbrewer2('Reds',tLevels));
% LampColorThresh = linspace(LampMin,LampMax,5);
% LampColorThresh = round(LampColorThresh);

tLevels = 12;
brew = cbrewer2('Reds',tLevels);
color_map=colormap(brew(2:12,:));
LampMin = min(LampMeasurementsWT);
LampMax = max(LampMeasurementsWT);
%close;

figure;
%currently NOT plotting dorsal lateral rays
for i = 1:36
    if selectData(i,2) > 1
%             if LampMeasurements(i) < LampColorThresh(2)
%                 cl = color_map(2,:);
%             elseif LampMeasurements(i) < LampColorThresh(3)
%                 cl = color_map(3,:);
%             elseif LampMeasurements(i) < LampColorThresh(4)
%                 cl = color_map(4,:);
%             elseif LampMeasurements(i) < LampColorThresh(5)
%                 cl = color_map(5,:);
%             else
%                 cl = 'g';
%             end
        LampHere = LampMeasurements(i);
        colorCode = round((((LampHere-LampMin)./(LampMax-LampMin))).*10);
        cl = color_map(colorCode+1,:);
        cl = horzcat(cl);
    plot(timePoints,selectData(i,3:15),'-','LineWidth',1,'color',cl); hold on;
    %colorCode
    %pause(1)
    end

    

end
ylabel('Fraction Regenerated');
xlabel('Time (days post amp.)');
set(gca,'fontsize',16);
yline(1,'--');
%ylim([0,3]);
ylim([0,1.2]); %changed 9Aug24

%linspace(LampMin,LampMax,11);


color_map=colormap(brew(2:12,:));
g = colorbar;
caxis([LampMin,LampMax]);
g.Label.String = 'Length Amputated (\mum)';

%Generate Supp. Fig. 5C
% Generate LOF timecourse plot

figure;

LampMeasurementsLOF = LampMeasurements(37:72,:);
% tLevels = 5;
% color_map=colormap(cbrewer2('Reds',tLevels));
% LampColorThresh = linspace(LampMin,LampMax,5);
% LampColorThresh = round(LampColorThresh);

%set up colors
tLevels = 12;
brew = cbrewer2('Reds',tLevels);
color_map=colormap(brew(2:12,:));
LampMin = min(LampMeasurementsLOF);
LampMax = max(LampMeasurementsLOF);
close;

figure;
for q = 37:72
    if ~isnan(LampMeasurements(q))
        if selectData(q,2) > 1
    %             if LampMeasurements(q) < LampColorThresh(2)
    %                 cl = color_map(2,:);
    %             elseif LampMeasurements(q) < LampColorThresh(3)
    %                 cl = color_map(3,:);
    %             elseif LampMeasurements(q) < LampColorThresh(4)
    %                 cl = color_map(4,:);
    %             elseif LampMeasurements(q) < LampColorThresh(5)
    %                 cl = color_map(5,:);
    %             else
    %                 cl = 'g';
    %                 selectData(q,1)
    %                 selectData(q,2)
    %             end
            LampHere = LampMeasurements(q);
            colorCode = round((((LampHere-LampMin)./(LampMax-LampMin))).*10);
            cl = color_map(colorCode+1,:);
            cl = horzcat(cl);
        plot(timePoints,selectData(q,3:15),'--','LineWidth',1,'color',cl); hold on;
        end 
    else
        continue
    end
end

ylabel('Fraction Regenerated');
xlabel('Time (days post amp.)');
set(gca,'fontsize',16);
yline(1,'--');
ylim([0,3]);

color_map=colormap(brew(2:12,:));
g = colorbar;
caxis([LampMin,LampMax]);
g.Label.String = 'Length Amputated (\mum)';


%% Supp. Fig. 5D

%load "analysis_mat_oldLOFerkOnlythresh15__newLOFthresh27.mat" from "SuppFig5"
%folder

% duplicate analysis matrix and rename as analysis_mat_timeaverage
analysis_mat_timeaverage = analysis_mat;

% add u values to analysis matrix
sAll = {};

for i = 1:size(analysis_mat_timeaverage,2)

    xVals = analysis_mat_timeaverage(i).ccrot;
    LregHere = analysis_mat_timeaverage(i).L_reg;
    u = xVals(:,1)./LregHere;
    sAll = [sAll,u];
end

newstruct = 'u';

[analysis_mat_timeaverage.(newstruct)] = sAll{:};


% aggregate single nuclei first - subtract later - bin data into 5 groups

% remove rays earlier than 72 hpa
% field_name = 'fit_ktr_linear';
% not_empty = arrayfun(@(x)~isempty(x.(field_name)),analysis_mat_timeaverage);

% E0 = 0.8483;
E0 = 0.8;

L_reg = [analysis_mat_timeaverage.L_reg];
L_amp = [analysis_mat_timeaverage.L_amp];
phi = L_reg./L_amp;
timeHere = [analysis_mat.hpa];

% bin average for single ray-time
averageKTR_u_cell = arrayfun(@(s)...
    adjust_neg(bin_average(s.u(s.trim_logical),...
    s.ktr(s.trim_logical),10,0,1)-E0),... %changed from 10 to 20
    analysis_mat_timeaverage,'UniformOutput',false);

mean_ktr_bin10_offsetE0 = cellfun(@nanmean,averageKTR_u_cell);

% normalize, e = Erk(u), m = mean(Erk)
norm_averageKTR_u_cell = arrayfun(@(e,m)(e{:})./(m),averageKTR_u_cell,mean_ktr_bin10_offsetE0,'UniformOutput',false);

norm_averageKTR_u = cell2mat(norm_averageKTR_u_cell);

%
% fit not-timeaveraged data with a line
% bin average
data.x = repmat([0.05:0.1:0.95]',[1,size(norm_averageKTR_u,2)]);
data.y = norm_averageKTR_u;
[binned.y,binned.x,binned.sem] = bin_average(data.x(:),data.y(:),10,0,1);
fit_x = linspace(0,1,100);

% Keep only data for first 14 dpas

%separateData by Lamp
bin_x_1amp =  [];
bin_y_1amp = [];
bin_x_2amp = [];
bin_y_2amp = [];
bin_x_3amp = [];
bin_y_3amp = [];
bin_x_4amp = [];
bin_y_4amp = [];
bin_x_5amp = [];
bin_y_5amp = [];
bin_x_allAmp = [];
bin_y_allAmp = [];
for i = 1:size(data.x,2)
    if timeHere(i) < 337
            bin_x_allAmp = horzcat(bin_x_allAmp,data.x(:,i));
            bin_y_allAmp = horzcat(bin_y_allAmp,data.y(:,i));
        if L_amp(i) < 3381
            bin_x_1amp= horzcat(bin_x_1amp,data.x(:,i));
            bin_y_1amp= horzcat(bin_y_1amp,data.y(:,i));
        elseif L_amp(i) < 5261
            bin_x_2amp = horzcat(bin_x_2amp,data.x(:,i));
            bin_y_2amp = horzcat(bin_y_2amp,data.y(:,i));
        elseif L_amp(i) < 7141
            bin_x_3amp = horzcat(bin_x_3amp,data.x(:,i));
            bin_y_3amp = horzcat(bin_y_3amp,data.y(:,i));
            elseif L_amp(i) < 9021
            bin_x_4amp = horzcat(bin_x_4amp,data.x(:,i));
            bin_y_4amp = horzcat(bin_y_4amp,data.y(:,i));
        else %10833
            bin_x_5amp = horzcat(bin_x_5amp,data.x(:,i));
            bin_y_5amp = horzcat(bin_y_5amp,data.y(:,i));
        end
    else
        continue
    end
end

%separateData by time
bin_x_1time =  [];
bin_y_1time = [];
bin_x_2time = [];
bin_y_2time = [];
bin_x_3time = [];
bin_y_3time = [];
bin_x_4time = [];
bin_y_4time = [];
bin_x_5time = [];
bin_y_5time = [];
for i = 1:size(data.x,2)
    if timeHere(i) < 111 %158
        bin_x_1time= horzcat(bin_x_1time,data.x(:,i));
        bin_y_1time= horzcat(bin_y_1time,data.y(:,i));
    elseif timeHere(i) < 171 %245
        bin_x_2time = horzcat(bin_x_2time,data.x(:,i));
        bin_y_2time = horzcat(bin_y_2time,data.y(:,i));
    elseif timeHere(i) < 231 %331
        bin_x_3time = horzcat(bin_x_3time,data.x(:,i));
        bin_y_3time = horzcat(bin_y_3time,data.y(:,i));
    elseif timeHere(i) < 291 %418
        bin_x_4time = horzcat(bin_x_4time,data.x(:,i));
        bin_y_4time = horzcat(bin_y_4time,data.y(:,i));
    elseif timeHere(i) < 340
        bin_x_5time = horzcat(bin_x_5time,data.x(:,i));
        bin_y_5time = horzcat(bin_y_5time,data.y(:,i));
    end
end

% plot LOF f(u) - color by Time - only use data for first 14 dpa
g2 = figure;

bb = errorbar(mean(bin_x_1time,2),nanmean(bin_y_1time,2),(nanstd(bin_y_1time,0,2)./sqrt(size(bin_y_1time,2))),'-','color',[0.8812    0.9286    0.9721],'markersize',5,'LineWidth',4); hold on;
cc = errorbar(mean(bin_x_2time,2),nanmean(bin_y_2time,2),(nanstd(bin_y_2time,0,2)./sqrt(size(bin_y_2time,2))),'-','color',[0.6744    0.8178    0.9018],'markersize',5,'LineWidth',4);
dd = errorbar(mean(bin_x_3time,2),nanmean(bin_y_3time,2),(nanstd(bin_y_3time,0,2)./sqrt(size(bin_y_3time,2))),'-','color',[0.3454    0.6344    0.8125],'markersize',5,'LineWidth',4);
ee = errorbar(mean(bin_x_4time,2),nanmean(bin_y_4time,2),(nanstd(bin_y_4time,0,2)./sqrt(size(bin_y_4time,2))),'-','color',[0.1116    0.4148    0.6906],'markersize',5,'LineWidth',4);
ff = errorbar(mean(bin_x_5time,2),nanmean(bin_y_5time,2),(nanstd(bin_y_5time,0,2)./sqrt(size(bin_y_5time,2))),'-','color',[0.0314    0.1882    0.4196],'markersize',5,'LineWidth',4);
aa = errorbar([0.05:0.1:0.95]',nanmean(bin_y_allAmp,2),(nanstd(bin_y_allAmp,0,2)./sqrt(size(bin_y_allAmp,2))),'-','color','k','markersize',5,'LineWidth',4); 
hold off

ylim([0,3])
xlabel('u (Normalized Position)')
ylabel('ERK Activity / Average ERK Activity')

set(gca,'XDir','reverse','fontsize',16)


legend_name={strcat("All Data (228)")...
    strcat("3-4.5 dpa (93)"),...
        strcat("5-7 dpa (55)"),...
        strcat("7.5-9.5 dpa (39)"),...
        strcat("10-12 dpa (25)"),...
        strcat("12.5+ dpa (16)")};
legend([aa,bb,cc,dd,ee,ff], legend_name,...
        'Location','northwest','color','none','box','off');