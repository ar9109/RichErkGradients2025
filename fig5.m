%% Fig. 5D

%load "binnedAvg.mat", "binned_wt15.mat", "binned_lofOld15.mat", and
%"binned_lofNew27.mat" from "Fig5" folder

figure;
wt15 = errorbar(binned_wt15.x(4:15)-0.8,binned_wt15.y(4:15),binned_wt15.sem(4:15),'-o','color','k','linewidth',2,'CapSize',10); hold on;
AVGnew27old15 = errorbar(binnedAvg.x(4:15)-0.8,binnedAvg.y(4:15),binnedAvg.sem(4:15),'-o','color','r','linewidth',2,'CapSize',10); hold on;
%new27 = errorbar(binned_lofNew27.x-0.8,binned_lofNew27.y,binned_lofNew27.sem,'-o','color',[0.3010 0.7450 0.9330],'linewidth',2,'CapSize',10); hold on;
%old15 = errorbar(binned_lofOld15.x-0.8,binned_lofOld15.y,binned_lofOld15.sem,'-o','color',[0.9290 0.6940 0.1250],'linewidth',2,'CapSize',10); hold on;

hold off;
%xlabel('averageERK(u)')
xlabel('Binned Average ERK Activity (A.U.)')
%ylabel('%GEM+(u)')
ylabel('Binned Fraction Cycling')
xlim([0,0.6]);
legend([wt15,AVGnew27old15], {['wildtype'],['longfin']},'Location','best','color','none','box','off');
% f.OuterPosition = [100,100,650,600];
set(gca, 'fontsize', 20,'linewidth',3')%,'Position' ,[0.1269    0.1561    0.6053    0.7689]);

%% Fig. 5E

%load "analysis_mat_oldLOFerkOnlythresh15__newLOFthresh27.mat" from "Fig5"
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


% plot LOF f(u) - color by Lamp - only use data for first 14 dpa


g = figure;

b = errorbar(mean(bin_x_1amp,2),nanmean(bin_y_1amp,2),(nanstd(bin_y_1amp,0,2)./sqrt(size(bin_y_1amp,2))),'-','color',[0.9966 0.8899 0.8395],'markersize',10,'LineWidth',4); hold on;
c = errorbar(mean(bin_x_2amp,2),nanmean(bin_y_2amp,2),(nanstd(bin_y_2amp,0,2)./sqrt(size(bin_y_2amp,2))),'-','color',[0.9884 0.6262 0.5059],'markersize',10,'LineWidth',4);
d = errorbar(mean(bin_x_3amp,2),nanmean(bin_y_3amp,2),(nanstd(bin_y_3amp,0,2)./sqrt(size(bin_y_3amp,2))),'-','color',[0.9746 0.3327 0.2330],'markersize',10,'LineWidth',4);
e = errorbar(mean(bin_x_4amp,2),nanmean(bin_y_4amp,2),(nanstd(bin_y_4amp,0,2)./sqrt(size(bin_y_4amp,2))),'-','color',[0.7664 0.0711 0.1052],'markersize',10,'LineWidth',4);
f = errorbar(mean(bin_x_5amp,2),nanmean(bin_y_5amp,2),(nanstd(bin_y_5amp,0,2)./sqrt(size(bin_y_5amp,2))),'-','color',[0.4039 0 0.0510],'markersize',10,'LineWidth',4);
a = errorbar([0.05:0.1:0.95]',nanmean(bin_y_allAmp,2),nanstd(bin_y_allAmp,0,2)./sqrt(size(bin_y_allAmp,2)),'-','color','k','markersize',10,'LineWidth',4); 
hold off;

ylim([0,3])
xlabel('u (Normalized Position)')
ylabel('ERK Activity / Average ERK Activity')

legend_name={strcat("All Data (228)")...
    strcat("1501-3380 \mum (100)"),...
        strcat("3381-5260 \mum (18)"),...
        strcat("5261-7240 \mum (29)"),...
        strcat("7141-9200 \mum (40)"),...
        strcat("9021-10900 \mum (41)")};
legend([a,b,c,d,e,f], legend_name,...
        'Location','northwest','color','none','box','off');


set(gca,'XDir','reverse','fontsize',16)


%% Fig. 5F

%load "analysis_mat_oldLOFerkOnlythresh15__newLOFthresh27.mat" from "Fig5"
%folder

% A(phi) color by L_AMP - fit individual fish - plot all on 1 plot

%set up colors
tLevels = 12;
brew = cbrewer2('Reds',tLevels);
%color_map=colormap(cbrewer2('Reds',tLevels));
color_map=colormap(brew(2:12,:));
LampMin = 1500;
LampMax = 10900;
close;

%newData
hpaUnique = [72,84,96,108,120,168,180,192,204,216,264,312,336,360,396,408,456,504];
fishUnique = [1,2,3,5,6,11,12,14,15];
rayUnique = [2,3,6,7,8];

% Added 20 Apr 23
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
% End added 20 apr 23

%ft = fittype('a*(1-x)+b');
ft = fittype({'1-x','1'});
foKTR = fitoptions('Method','LinearLeastSquares',...
    'Lower',[0,-inf],...
    'Weights', []);

f = figure('Visible','on');
hold on;

% plot individual fish
    for fishNum = 1:size(fishUnique,2)
        fishHere1 = fishUnique(fishNum);
        for rayNum = 1:size(rayUnique,2)
            rayHere1 = rayUnique(rayNum);
            xdataCollect = [];
            ydataCollect = [];
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
                    else
                        continue
                    end
                end
            %end

            %if ~isempty(xdataCollect)
            if size(xdataCollect,1)>1

            xdata = xdataCollect;
            ydata = ydataCollect;


            [fit_obj,fit_gof] = fit(xdata(:),ydata(:),ft,foKTR);
        
            phiMax = max(xdata);
            
            fit_x = linspace(0,phiMax,100);
            
            
            phiFinal = xdata(size(xdata,1));

            colorCode = round((((LampHere-LampMin)./LampMax)).*10);

            cl = color_map(colorCode+1,:);

            cl = horzcat(cl,0.5);

            p = plot(xdata,ydata,'.','Color',cl,'MarkerSize',15); hold on;
            
            p_fit = plot_fit(fit_obj,fit_x,fit_gof,confidence=false); hold on;
            set(p_fit,'linestyle','-','color',cl,'lineWidth',3)
        
            else 
                continue
            end
            xdataCollect = [];
            ydataCollect = [];
% % %             close all;
        end
    end

xlabel('Fraction Regenerated')
ylabel('Average ERK Activity (rescaled space,u)')
config_plot(f);
ylim([0,0.6]);
xlim([0,2]);
legend('off');


color_map=colormap(brew(2:12,:));



g = colorbar;
caxis([LampMin,LampMax]);

g.Label.String = 'Length Amputated (\mum)';

%% Fig. 5G

%load "analysis_mat_oldLOFerkOnlythresh15__newLOFthresh27.mat" from "Fig5"
%folder

%step 1
% Extract Lfinal measurements

L_reg = [analysis_mat.L_reg];
L_amp = [analysis_mat.L_amp];
phi = L_reg./[L_amp];
ray = [analysis_mat.ray];
fish = [analysis_mat.fish];
time = [analysis_mat.hpa];

mergedData = horzcat(fish', ray', time', L_amp', L_reg');

fishValues = [1 2 3 5 6 11 12 14 15]';
rayValues = [2 3 6 7 8]';

finalCollect = [];
for q = 1:size(fishValues)
    finalValues = [];
    

    for j = 1:size(rayValues)
        LregOneRay = [];
        LampOneRay = [];
        timeOneRay = [];
        for i = 1:size(mergedData)
        if mergedData(i,1) == fishValues(q) && mergedData(i,2) == rayValues(j)
            LregHere = mergedData(i,5);
            LampHere = mergedData(i,4);
            timeHere = mergedData(i,3);
            LregOneRay = vertcat(LregOneRay,LregHere);
            LampOneRay = vertcat(LampOneRay,LampHere);
            timeOneRay = vertcat(timeOneRay,timeHere);
        end

        end
        if ~isempty(LregOneRay) && max(timeOneRay) > 300
            Lfinal = max(LregOneRay);
            fishFinal = fishValues(q);
            rayFinal = rayValues (j);


        finalValues = horzcat(fishFinal,rayFinal,Lfinal);
        finalCollect = vertcat(finalCollect,finalValues);
        else
            continue
        end

    end
end

%step 2
% A(Lreg/Lfinal) color by L_AMP - fit individual fish - plot all on 1 plot  No Lines

paths.plotFolder = '/Volumes/AshleyData2/29Jan23_lof_combined_24May23/plots';
%mkdir(paths.plotFolder);

%set up colors
tLevels = 12;
brew = cbrewer2('Reds',tLevels);
%color_map=colormap(cbrewer2('Reds',tLevels));
color_map=colormap(brew(2:12,:));
LampMin = 1500;
LampMax = 10900;
close;

%newData
hpaUnique = [72,84,96,108,120,168,180,192,204,216,264,312,336,360,396,408,456,504];

fishUnique = [1,2,3,5,6,11,12,14,15];

rayUnique = [2,3,6,7,8];

% Added 20 Apr 23
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
% End added 20 apr 23

%ft = fittype('a*(1-x)+b');
ft = fittype({'1-x','1'});
foKTR = fitoptions('Method','LinearLeastSquares',...
    'Lower',[0,-inf],...
    'Weights', []);

f = figure('Visible','on');
hold on;

xdataFit = [];
ydataFit = [];

% plot individual fish
    for fishNum = 1:size(fishUnique,2)
        fishHere1 = fishUnique(fishNum);
        for rayNum = 1:size(rayUnique,2)
            rayHere1 = rayUnique(rayNum);
            xdataCollect = [];
            ydataCollect = [];
%             for hpaNum = 1:size(hpaUnique,2)
%                 hpaHere = hpaUnique(hpaNum);
                for aa = 1:size(analysis_mat,2)
                    if analysis_mat(aa).fish == fishHere1 && analysis_mat(aa).ray == rayHere1
                    fishHere = analysis_mat(aa).fish;
                    rayHere = analysis_mat(aa).ray;
                    timeHere = analysis_mat(aa).hpa;
                    LampHere = analysis_mat(aa).L_amp;
                    LregHere = analysis_mat(aa).L_reg;
                    %xdataHere = phi(aa);
                    
                    LfinalHere = [];
                    for qqq = 1:size(finalCollect)
                        
                        if fishHere == finalCollect(qqq,1) && rayHere == finalCollect(qqq,2)
                            LfinalHere = finalCollect(qqq,3);
                        else
                            continue
                        end
                    end


                    xdataHere = LregHere./LfinalHere;
                    
                    ydataHere = mean_ktr(aa)-E0;
                    xdataCollect = vertcat(xdataCollect,xdataHere);
                    ydataCollect = vertcat(ydataCollect,ydataHere);
                    else
                        continue
                    end
                end
            %end

            %if ~isempty(xdataCollect)
            if size(xdataCollect,1)>1

            xdata = xdataCollect;
            ydata = ydataCollect;
        
            phiMax = max(xdata);
            
            %fit_x = linspace(0,phiMax,100);
            fit_x = linspace(0,1,100);
            

            
            phiFinal = xdata(size(xdata,1));

            colorCode = round((((LampHere-LampMin)./LampMax)).*10);

            cl = color_map(colorCode+1,:);

            cl = horzcat(cl,0.5);

            p = plot(xdata,ydata,'.','Color',cl,'MarkerSize',15); hold on;
            
        
            else 
                continue
            end
            xdataFit = vertcat(xdataFit,xdataCollect);
            ydataFit = vertcat(ydataFit,ydataCollect);
            xdataCollect = [];
            ydataCollect = [];

% % %             close all;
        end
    end

[fit_obj,fit_gof] = fit(xdataFit(:),ydataFit(:),ft,foKTR);
fit_x = linspace(0,1,100);   
p_fit = plot_fit(fit_obj,fit_x,fit_gof,confidence=false); hold on;
set(p_fit,'linestyle','-','color','k','lineWidth',3)

xlabel('L_{reg} / L_{final}')
ylabel('Average ERK Activity (rescaled space,u)')
config_plot(f);
ylim([0,0.6]);

color_map=colormap(brew(2:12,:));



g = colorbar;
caxis([LampMin,LampMax]);

g.Label.String = 'Length Amputated (\mum)';

legend_name={strcat("R^2=",num2str(round(fit_gof.rsquare,1)))};
legend([p_fit], legend_name,...
        'Location','northwest','color','none','box','off');

