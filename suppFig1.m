%% Supp Fig 1a

%load "timeDataEarly", "timeDataLate", "lateralDataEarly",
%"lateralDataLate", "medialDataEarly", and "medialDataLate" from "SuppFig1"
%folder

% Plot average curves for lateral and medial rays from early and late experiment - Updated

paths=[];
% paths.masterFolder='/Volumes/AshleyData2/23Jan23_CombineCellCycleData/'; %folder where data is stored
% paths.outFolder=   [paths.masterFolder 'figures/']; % output folder quantifications

figure;

errorbar(timeDataEarly-96, nanmean(lateralDataEarly,2),nanstd(lateralDataEarly,0,2)./sqrt(numel(lateralDataEarly)),'.','linewidth',4,'markersize',30,'color',[0.4039 0 0.0510]); hold on;
errorbar(timeDataEarly-96, nanmean(medialDataEarly,2),nanstd(medialDataEarly,0,2)./sqrt(numel(medialDataEarly)),'.','linewidth',4,'markersize',30,'color',[0.9897 0.4865 0.3549]); hold on;

errorbar(timeDataLate-192, nanmean(lateralDataLate,2),nanstd(lateralDataLate,0,2)./sqrt(numel(lateralDataLate)),'x','linewidth',4,'markersize',15,'color',[0.4039 0 0.0510]); hold on;
errorbar(timeDataLate-192, nanmean(medialDataLate,2),nanstd(medialDataLate,0,2)./sqrt(numel(medialDataLate)),'x','linewidth',4,'markersize',15,'color',[0.9897 0.4865 0.3549]); hold on;



% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%opts.StartPoint = [0.196595250431208 0.251083857976031];
% Fit model to data.
[fitresult_latEarly, gof_latEarly] = fit((timeDataEarly-96)', nanmean(lateralDataEarly,2), ft, opts );
[fitresult_medEarly, gof_medEarly] = fit((timeDataEarly-96)', nanmean(medialDataEarly,2), ft, opts );
[fitresult_latLate, gof_latLate] = fit((timeDataLate-192)', nanmean(lateralDataLate,2), ft, opts );
[fitresult_medLate, gof_medLate] = fit((timeDataLate-192)', nanmean(medialDataLate,2), ft, opts );

a = plot(fitresult_latEarly);
set(a,'linestyle','-','linewidth',2,'color',[0.4039 0 0.0510]);
b = plot(fitresult_medEarly);
set(b,'linestyle','-','linewidth',2,'color',[0.9897 0.4865 0.3549]);
c= plot(fitresult_latLate);
set(c,'linestyle','--','linewidth',2,'color',[0.4039 0 0.0510]);
d = plot(fitresult_medLate);
set(d,'linestyle','--','linewidth',2,'color',[0.9897 0.4865 0.3549]);

set(gca, 'fontsize',16);
xlabel('Time Post Inhibition (hours)');
ylabel('Average Fraction Proliferating');
xlim([0,24]);

legend('4dpa Lateral','4dpa Medial','8dpa Lateral','8dpa Medial');

%% Supp Fig 1d
% load analysis_mat_threshold.mat from "Supp Fig 1" folder
t_bins = [];
for numfish = unique([analysis_mat.fish])
    for numray = unique([analysis_mat.ray])
        % create cropped myScale_merged struct
        fishNray_data = analysis_mat([analysis_mat.fish]==numfish & [analysis_mat.ray]==numray);
        if numel(fishNray_data)<2
            continue
        end

        disp('------------');
        display(['fish' num2str(numfish) '_ray' num2str(numray)]);
        
        single_ray = [];
        single_ray.fish = numfish;
        single_ray.ray = numray;

        if single_ray.fish == 2 && single_ray.ray ==2

        
        if ~isempty(fishNray_data)
            single_ray.hpa = [fishNray_data.hpa];
            single_ray.hpaTrue = [fishNray_data.hpaTrue];
            single_ray.L_reg = [fishNray_data.L_reg];
            single_ray.L_amp = [fishNray_data.L_amp];
            % time course matrix
            % binValue, average, sem
            t_mat_cell = arrayfun(@(x)flip(x.averageKTR),fishNray_data,...
                'UniformOutput', false);
            t_bins_mat = padcat(t_mat_cell{:});

            x_mat_cell = arrayfun(@(x)flip(x.binvalue),fishNray_data,...
                'UniformOutput', false);
            x_bins_mat = padcat(x_mat_cell{:});
            
            % plot            
            f = figure('visible','on');
            cm = colormap(lines(size(t_bins_mat,2)));
            tLevels = 40;
            color_map=colormap(cbrewer2('Greys',tLevels));
            for i = 1:size(t_bins_mat,2)

            %plot(single_ray.hpa,t_bins_mat(:,i)./mean(t_bins_mat(:,i)),'-','color',cm(i,:),'linewidth',2,'markersize',10);hold on;
            %plot(single_ray.hpa,(t_bins_mat(:,i)./mean(t_bins_mat(:,i)))-0.8,'.-','color','k','linewidth',2,'markersize',10);hold on;
            plot(single_ray.hpa,(t_bins_mat(:,i)./mean(t_bins_mat(:,i)))-0.8,'.-','color',color_map(i+10,:),'linewidth',2,'markersize',10);hold on;
            end
            %colorbar('Ticks',[0,1],'TickLabels', {'Amp','End'});
%             legend({'t1','t2','t3','t4'},...
%                    'FontSize',8,'Location','best')

%             xticks([120:12:180])
            ylabel('Average KTR');
            xlabel('Time (hpa)');
%             xlim([0 2000]);
%             ylim([0 1]);
            title(['fish' num2str(numfish) '_ray' num2str(numray)],'Interpreter', 'none');
%             set(gca,'XDir','reverse')
            set(gca,'fontsize',16)
            hold off;

            %saveas(f,[paths.plotFolder,filesep,['fish' num2str(numfish) '_ray' num2str(numray)],'_tbins_normalized.png']);
            %close(f);
            
            single_ray.t_bins_mat = t_bins_mat;
            single_ray.t_mat_cell = t_mat_cell;
            single_ray.x_bins_mat = x_bins_mat;
            single_ray.x_mat_cell = x_mat_cell;
            t_bins = [t_bins single_ray];

        else 
            continue
        end
        end
    end
end

%% Supp. Fig. 1E
% load "analysis_mat_threshold_8h.mat" from folder "SuppFig1"

%step 1
% oscillation of bin ktr against time

t_bins = [];
for numfish = unique([analysis_mat.fish])
    for numray = unique([analysis_mat.ray])
        % create cropped myScale_merged struct
        fishNray_data = analysis_mat([analysis_mat.fish]==numfish & [analysis_mat.ray]==numray);
        if numel(fishNray_data)<2
            continue
        end

        disp('------------');
        display(['fish' num2str(numfish) '_ray' num2str(numray)]);
        
        single_ray = [];
        single_ray.fish = numfish;
        single_ray.ray = numray;

        
        if ~isempty(fishNray_data)
            single_ray.hpa = [fishNray_data.hpa];
            single_ray.hpaTrue = [fishNray_data.hpaTrue];
            single_ray.L_reg = [fishNray_data.L_reg];
            single_ray.L_amp = [fishNray_data.L_amp];
            % time course matrix
            % binValue, average, sem
            t_mat_cell = arrayfun(@(x)flip(x.averageKTR),fishNray_data,...
                'UniformOutput', false);
            t_bins_mat = padcat(t_mat_cell{:});

            x_mat_cell = arrayfun(@(x)flip(x.binvalue),fishNray_data,...
                'UniformOutput', false);
            x_bins_mat = padcat(x_mat_cell{:});
            
            
            single_ray.t_bins_mat = t_bins_mat;
            single_ray.t_mat_cell = t_mat_cell;
            single_ray.x_bins_mat = x_bins_mat;
            single_ray.x_mat_cell = x_mat_cell;
            t_bins = [t_bins single_ray];
        end
    end
end

%step 2
% 2 create lag matrix REAL intervals
fields = {'lag','ERK_pairs','ERK_pairs_detrend_meanERK','ERK_pairs_detrend_linearERK'};
lag_mat = cell2struct(cell(length(fields),1),fields);
lag_mat = [];

% dt = 7;
% 
% ERK_lag_here = [];
% ERK_lag_here_detrend_meanERK = [];
% ERK_lag_here_detrend_linearERK = [];
% ERK_lag_here_detrend_AphiFu = [];
for i = 1:numel(t_bins)
    t_bins_mat_here = t_bins(i).t_bins_mat;
    x_bins_mat_here = t_bins(i).x_bins_mat;
    hpa_here = t_bins(i).hpa;
    hpaTrue_here = t_bins(i).hpaTrue;

%         t_bins_mat_here = t_bins_mat_here(:,[1:6]);
    
    t_bins_detrend_meanERK = t_bins_mat_here-nanmean(t_bins_mat_here,2);
    t_bins_detrend_linearERK = nan(size(t_bins_mat_here));
    for rowERK = 1:size(t_bins_mat_here,1)
        yfit = t_bins_mat_here(rowERK,:);
        xfit = x_bins_mat_here(rowERK,:);
        datafit = remove_nan([xfit(:),yfit(:)]);
        if size(datafit,1)>2
            fit_obj = fit(datafit(:,1),datafit(:,2),fittype({'x','1'}));
%                 fit_obj = fit(datafit(:,1),datafit(:,2),'smoothingspline',SmoothingParam = 0.0000001);
%                 %%
%                 figure;
%                 plot(fit_obj,datafit(:,1),datafit(:,2))
%                 %%
            t_bins_detrend_linearERK(rowERK,~isnan(t_bins_mat_here(rowERK,:))) = ...
                t_bins_mat_here(rowERK,~isnan(t_bins_mat_here(rowERK,:)))...
                -fit_obj(datafit(:,1))';
        end
    end


    for col = 1:size(t_bins_mat_here,2)
        colHere = t_bins_mat_here(:,col);
        colHere_detrend_meanERK = t_bins_detrend_meanERK(:,col);
        colHere_detrend_linearERK = t_bins_detrend_linearERK(:,col);
        for row = 1:numel(colHere)
            for rowT = row:numel(colHere)

                
            ERK_t0 = colHere(row);
            ERK_tdt = colHere(rowT);
            lagHerePrecise = hpaTrue_here(rowT)-hpaTrue_here(row);
            lagHere = round(lagHerePrecise);
            if isempty(ERK_t0)||isnan(ERK_t0) % skip empty
                continue
            end
            if isempty(ERK_tdt)||isnan(ERK_tdt) % skip empty
                continue
            end
            % detrend by subtracting A(phi)f(u)
%                 LregHere_t0 = t_bins(i).L_reg(row);
%                 LampHere_t0 = t_bins(i).L_amp(row);
%                 phiHere_t0 = LregHere_t0./LampHere_t0;
%                 xHere_t0 = t_bins(i).x_bins_mat(row,col);
%                 uHere_t0 = (LregHere_t0-xHere_t0)./LregHere_t0;
%                 ERK_t0_detrend_AphiFu = ERK_t0-...
%                     (0.80824+0.1829.*(1-phiHere_t0)+0.39449.*(1-phiHere_t0).*uHere_t0);
            
%                 LregHere_tdt = t_bins(i).L_reg(row+lag);
%                 LampHere_tdt = t_bins(i).L_amp(row+lag);
%                 phiHere_tdt = LregHere_tdt./LampHere_tdt;
%                 xHere_tdt = t_bins(i).x_bins_mat(row+lag,col);
%                 uHere_tdt = (LregHere_tdt-xHere_tdt)./LregHere_tdt;
%                 ERK_tdt_detrend_AphiFu = ERK_tdt-...
%                     (0.80824+0.1829.*(1-phiHere_tdt)+0.39449.*(1-phiHere_tdt).*uHere_tdt);

            % append to pairs-array
            if isfield(lag_mat,'lag') && ismember(lagHere,[lag_mat.lag])
                lag_mat(lagHere==[lag_mat.lag]).lag = lagHere;
                lag_mat(lagHere==[lag_mat.lag]).ERK_pairs = ...
                    [lag_mat(lagHere==[lag_mat.lag]).ERK_pairs;...
                    [ERK_t0,ERK_tdt]];
                lag_mat(lagHere==[lag_mat.lag]).ERK_pairs_detrend_meanERK = ...
                    [lag_mat(lagHere==[lag_mat.lag]).ERK_pairs_detrend_meanERK;...
                    [colHere_detrend_meanERK(row),colHere_detrend_meanERK(rowT)]];
                lag_mat(lagHere==[lag_mat.lag]).ERK_pairs_detrend_linearERK = ...
                    [lag_mat(lagHere==[lag_mat.lag]).ERK_pairs_detrend_linearERK;...
                    [colHere_detrend_linearERK(row),colHere_detrend_linearERK(rowT)]];
            else
                lag_mat(numel(lag_mat)+1).lag = lagHere;
                lag_mat(lagHere==[lag_mat.lag]).ERK_pairs = [ERK_t0,ERK_tdt];
                lag_mat(lagHere==[lag_mat.lag]).ERK_pairs_detrend_meanERK = ...
                    [colHere_detrend_meanERK(row),colHere_detrend_meanERK(rowT)];
                lag_mat(lagHere==[lag_mat.lag]).ERK_pairs_detrend_linearERK = ...
                    [colHere_detrend_linearERK(row),colHere_detrend_linearERK(rowT)];
            end



            
%             ERK_lag_here = [ERK_lag_here;[ERK_t0,ERK_tdt]];
%             ERK_lag_here_detrend_meanERK = [ERK_lag_here_detrend_meanERK;...
%                 [colHere_detrend_meanERK(row),colHere_detrend_meanERK(row+lag)]];
%             ERK_lag_here_detrend_linearERK = [ERK_lag_here_detrend_linearERK;...
%                 [colHere_detrend_linearERK(row),colHere_detrend_linearERK(row+lag)]];
%                 ERK_lag_here_detrend_AphiFu = [ERK_lag_here_detrend_AphiFu;...
%                     [ERK_t0_detrend_AphiFu,ERK_tdt_detrend_AphiFu]];
            end
        end
    end
end
% sort
[~,I] = sort([lag_mat.lag]);
lag_mat = lag_mat(I);
% calculate corr
for i = 1:numel(lag_mat)
    corr_mat = corr(lag_mat(i).ERK_pairs);
    lag_mat(i).corr = corr_mat(1,2);
    corr_mat_detrend_meanERK = corr(lag_mat(i).ERK_pairs_detrend_meanERK);
    lag_mat(i).corr_detrend_meanERK = corr_mat_detrend_meanERK(1,2);    
    corr_mat_detrend_linearERK = corr(lag_mat(i).ERK_pairs_detrend_linearERK);
    lag_mat(i).corr_detrend_linearERK = corr_mat_detrend_linearERK(1,2);
end

%step 3
% plot autocorr
xplot = [lag_mat.lag];
sizeplot = arrayfun(@(s)size(s.ERK_pairs,1),lag_mat);

% yplot = [lag_mat.corr];
yplot = [lag_mat.corr_detrend_meanERK];
% yplot = [lag_mat.corr_detrend_linearERK];

% yplot = [lag_mat.corr_detrend_AphiFu];


% fit spline
fit_obj = fit(xplot',yplot','smoothingspline','SmoothingParam',0.1);
yy = fit_obj(xplot);

% yy = smooth(xplot,yplot,9,'sgolay',3);
% yy = smooth(xplot,yplot,11,'lowess');

f = figure;
cm = lines(11);

colororder([cm(1,:);cm(3,:)])

% yyaxis left
% plot(xplot,yplot,'-o',color=cm(1,:),LineWidth=3,MarkerSize=10)
errorbar(xplot,yplot,sqrt((1-yplot.^2)./(sizeplot-2)),LineWidth=3,MarkerSize=10)
hold on;
plot(xplot,yy,'-',color=cm(2,:),LineWidth=3,MarkerSize=10)

% confidence
% pred_bound = predint(fit_obj,xplot,0.95,'function','on');
% p_bound = plot(xplot,pred_bound,'--','LineWidth',2, 'Color',[0.5,0.5,0.5]);

xlabel('Time (h)')
ylabel('Temporal Autocorrelation')
% ylim([0,inf])
xlim([0,65])
ylim([0,1])

% yyaxis right
% plot(xplot,sizeplot,'-o',color=cm(3,:),LineWidth=3,MarkerSize=10)
% ylabel('Number')

% title('no detrend')
%title('detrend by <ERK>')
% title('detrend by linearERK')
% title('detrend by spline-7ERK')

% title('detrend by A(\phi)f(u)')
config_plot(f)

%% Supp. Fig. 1F

% load "analysis_mat_threshold_8h.mat" from folder "SuppFig1"

%step 1
% oscillation of bin ktr against time

t_bins = [];
for numfish = unique([analysis_mat.fish])
    for numray = unique([analysis_mat.ray])
        % create cropped myScale_merged struct
        fishNray_data = analysis_mat([analysis_mat.fish]==numfish & [analysis_mat.ray]==numray);
        if numel(fishNray_data)<2
            continue
        end

        disp('------------');
        display(['fish' num2str(numfish) '_ray' num2str(numray)]);
        
        single_ray = [];
        single_ray.fish = numfish;
        single_ray.ray = numray;

        
        if ~isempty(fishNray_data)
            single_ray.hpa = [fishNray_data.hpa];
            single_ray.hpaTrue = [fishNray_data.hpaTrue];
            single_ray.L_reg = [fishNray_data.L_reg];
            single_ray.L_amp = [fishNray_data.L_amp];
            % time course matrix
            % binValue, average, sem
            t_mat_cell = arrayfun(@(x)flip(x.averageKTR),fishNray_data,...
                'UniformOutput', false);
            t_bins_mat = padcat(t_mat_cell{:});

            x_mat_cell = arrayfun(@(x)flip(x.binvalue),fishNray_data,...
                'UniformOutput', false);
            x_bins_mat = padcat(x_mat_cell{:});
            
            
            single_ray.t_bins_mat = t_bins_mat;
            single_ray.t_mat_cell = t_mat_cell;
            single_ray.x_bins_mat = x_bins_mat;
            single_ray.x_mat_cell = x_mat_cell;
            t_bins = [t_bins single_ray];
        end
    end
end

%step 2
% create shift matrix (and plot) - plot different lags together (AR) 
dx = 50;
dt = 8;

for lag = 0:4
%    lag

% lag = 4; % lag*dt
shift_mat = [];



for shift = 0:max(cell2mat(arrayfun(@(s)(max(s.x_bins_mat,[],2)-min(s.x_bins_mat,[],2))'./dx,t_bins,uni=0))) % calculate max space step
    ERK_shift_here = [];
    ERK_shift_here_detrend_meanERK = [];
    ERK_shift_here_detrend_linearERK = [];
%     ERK_shift_here_detrend_AphiFu = [];
    for i = 1:numel(t_bins)
        t_bins_mat_here = t_bins(i).t_bins_mat;
        x_bins_mat_here = t_bins(i).x_bins_mat;
        hpa_here = t_bins(i).hpa;

%         t_bins_mat_here = t_bins_mat_here(:,[1:6]);
        
        t_bins_detrend_meanERK = t_bins_mat_here-nanmean(t_bins_mat_here,2);
        t_bins_detrend_linearERK = nan(size(t_bins_mat_here));
        for rowERK = 1:size(t_bins_mat_here,1)
            yfit = t_bins_mat_here(rowERK,:);
            xfit = x_bins_mat_here(rowERK,:);
            datafit = remove_nan([xfit(:),yfit(:)]);
            if size(datafit,1)>2
                fit_obj = fit(datafit(:,1),datafit(:,2),fittype({'x','1'}));
%                 fit_obj = fit(datafit(:,1),datafit(:,2),'smoothingspline',SmoothingParam = 0.0000001);
%                 %%
%                 figure;
%                 plot(fit_obj,datafit(:,1),datafit(:,2))
%                 %%
                t_bins_detrend_linearERK(rowERK,~isnan(t_bins_mat_here(rowERK,:))) = ...
                    t_bins_mat_here(rowERK,~isnan(t_bins_mat_here(rowERK,:)))...
                    -fit_obj(datafit(:,1))';
            end
        end


        for row = 1:size(t_bins_mat_here,1)
            rowHere = t_bins_mat_here(row,:);
            rowHere_detrend_meanERK = t_bins_detrend_meanERK(row,:);
            rowHere_detrend_linearERK = t_bins_detrend_linearERK(row,:);

            if row+lag>size(t_bins_mat_here,1) % find different or the same time point for crosscorr
                continue
            end
            rowTHere = t_bins_mat_here(row+lag,:);
            rowTHere_detrend_meanERK = t_bins_detrend_meanERK(row+lag,:);
            rowTHere_detrend_linearERK = t_bins_detrend_linearERK(row+lag,:);


            for col = 1:numel(rowHere)
                ERK_x0 = rowHere(col);
                if isempty(ERK_x0)||isnan(ERK_x0) % skip empty
                    continue
                end
                if col+shift>numel(rowHere) % skip when timecourse not long enough
                    continue
                end
                ERK_xdx = rowTHere(col+shift);
                if isempty(ERK_xdx)||isnan(ERK_xdx) % skip empty
                    continue
                end


                % append to pairs-array
                ERK_shift_here = [ERK_shift_here;[ERK_x0,ERK_xdx]];
                ERK_shift_here_detrend_meanERK = [ERK_shift_here_detrend_meanERK;...
                    [rowHere_detrend_meanERK(col),rowTHere_detrend_meanERK(col+shift)]];
                ERK_shift_here_detrend_linearERK = [ERK_shift_here_detrend_linearERK;...
                    [rowHere_detrend_linearERK(col),rowTHere_detrend_linearERK(col+shift)]];
%                 ERK_shift_here_detrend_AphiFu = [ERK_shift_here_detrend_AphiFu;...
%                     [ERK_t0_detrend_AphiFu,ERK_tdt_detrend_AphiFu]];
            end
        end
    end
    shift_mat_here.shift = shift.*dx;
    shift_mat_here.ERK_pairs = ERK_shift_here;
    shift_mat_here.ERK_pairs_detrend_meanERK = ERK_shift_here_detrend_meanERK;
    shift_mat_here.ERK_pairs_detrend_linearERK = ERK_shift_here_detrend_linearERK;
%     shift_mat_here.ERK_pairs_detrend_AphiFu = ERK_lag_here_detrend_AphiFu;

    corr_mat = corr(ERK_shift_here);
    shift_mat_here.corr = corr_mat(1,2);
    corr_mat_detrend_meanERK = corr(ERK_shift_here_detrend_meanERK);
    shift_mat_here.corr_detrend_meanERK = corr_mat_detrend_meanERK(1,2);    
    corr_mat_detrend_linearERK = corr(ERK_shift_here_detrend_linearERK);
    shift_mat_here.corr_detrend_linearERK = corr_mat_detrend_linearERK(1,2);


%     corr_mat_detrend_AphiFu = corr(ERK_lag_here_detrend_AphiFu);

%     lag_mat_here.corr_detrend_AphiFu = corr_mat_detrend_AphiFu(1,2);

    shift_mat = [shift_mat,shift_mat_here];
end
%end

% plot autocorr
xplot = [shift_mat.shift];
sizeplot = arrayfun(@(s)size(s.ERK_pairs,1),shift_mat);

yplot = [shift_mat.corr];
% yplot = log([shift_mat.corr]);
% yplot = [shift_mat.corr_detrend_meanERK];
% yplot = [shift_mat.corr_detrend_linearERK];

% yplot = [shift_mat.corr_detrend_AphiFu];



% fit spline
datafit = remove_nan([xplot(:),yplot(:)]);
datafit = datafit(datafit(:,1)<1000,:);
% fit_obj = fit(datafit(:,1),datafit(:,2),'smoothingspline','SmoothingParam',0.1);
fit_obj = fit(datafit(:,1),datafit(:,2),{'x' '1'});

yy = fit_obj(datafit(:,1));

% yy = smooth(xplot,yplot,9,'sgolay',3);
% yy = smooth(xplot,yplot,11,'lowess');

f = figure;
cm = lines(11);

% colororder([cm(1,:);cm(3,:)])





% yyaxis left
% plot(xplot,yplot,'-o',color=cm(1,:),LineWidth=3,MarkerSize=10)
errorbar(xplot,yplot,sqrt((1-yplot.^2)./(sizeplot-2)),LineWidth=3,MarkerSize=10)

% Collect values for later plotting
if lag == 0
    xplot0 = xplot;
    yplot0 = yplot;
    sizeplot0 = sizeplot;
elseif lag == 1
    xplot8 = xplot;
    yplot8 = yplot;
    sizeplot8 = sizeplot;
elseif lag == 2
    xplot16 = xplot;
    yplot16 = yplot;
    sizeplot16 = sizeplot;
elseif lag == 3
    xplot24 = xplot;
    yplot24 = yplot;
    sizeplot24 = sizeplot;
end

hold on;
% plot(datafit(:,1),yy,'-',color=cm(2,:),LineWidth=3,MarkerSize=10)
hold off;

xlabel('Shift (\mum)')
ylabel('log(Autocorr)')
ylim([-0.5,1])
% xlim([0,65])

% yyaxis right
% plot(xplot,sizeplot,'-o',color=cm(3,:),LineWidth=3,MarkerSize=10)
% ylabel('Number')

% title(['detrend <ERK>, lag = ' num2str(lag.*dt) 'h'])
title(['no detrend, lag = ' num2str(lag.*dt) 'h'])
% title('detrend by <ERK>')
% title('detrend by linearERK')
% title('detrend by spline-7ERK')

% title('detrend by A(\phi)f(u)')

% legend({'data' 'fit linear'})
config_plot(f)
end

%step 3
% Plot spatial correlations together

tLevels = 5;
color_map=colormap(cbrewer2('Blues',tLevels));

f = figure;


a = errorbar(xplot0,yplot0,sqrt((1-yplot0.^2)./(sizeplot0-2)),'linewidth',2,'markersize',10,'color',color_map(2,:)); hold on;
b = errorbar(xplot8,yplot8,sqrt((1-yplot8.^2)./(sizeplot8-2)),'linewidth',2,'markersize',10,'color',color_map(3,:));
c = errorbar(xplot16,yplot16,sqrt((1-yplot16.^2)./(sizeplot16-2)),'linewidth',2,'markersize',10,'color',color_map(4,:));
d = errorbar(xplot24,yplot24,sqrt((1-yplot24.^2)./(sizeplot24-2)),'linewidth',2,'markersize',10,'color',color_map(5,:));

xlabel('Distance (\mum)')
ylabel('Spatial Cross-correlation')
ylim([-0.5,1])
xlim([0,1000]) % added by ashley 22Oct24

config_plot(f)

legend_name={strcat("No Lag")...
    strcat("8hr Lag"),...
        strcat("16hr Lag"),...
        strcat("24hr Lag")};
legend([a,b,c,d], legend_name,...
        'Location','southwest','color','none','box','off');



%% Supp. Fig. 1H
% Plot GEM (500um tip) v. ERK (500um Tip) color map Greens

%load "ERKtip" and "GemPlusTip" from "SuppFig1" folder

% set up color code data

colorCodeCollect = [];
colorCode = [];

for i = 1:size(ERKtip)
    if ERKtip(i,2) == 1 || ERKtip(i,2) == 2
        cl = 10;
    elseif ERKtip(i,2) == 3 || ERKtip(i,2) == 4
        cl = 5; 
    elseif ERKtip(i,2) == 15 || ERKtip(i,2) == 16
        cl = 4;
    elseif ERKtip(i,2) == 11 || ERKtip(i,2) == 12
        cl = 3;
    elseif ERKtip(i,2) == 5 || ERKtip(i,2) == 6
        cl = 2;
    elseif ERKtip(i,2) == 7 || ERKtip(i,2) == 8 || ERKtip(i,2) == 13 || ERKtip(i,2) == 14
        cl = 1;
    end
    
    colorCode(i,1) = cl;
    
    %colorCodeCollect = vertcat(colorCodeCollect,colorCode);
end

tLevels = 10;
color_map=colormap(cbrewer2('Greens',tLevels));

for q = 1:size(ERKtip,1)
    
    plot(ERKtip(q,1),GemPlusTip(q),'.','color',color_map(colorCode(q),:),'MarkerSize',15); hold on;

end


plot(ERKtip(:,1),GemPlusTip(:),'o','color',[.42 .42 .42],'MarkerSize',5); hold on;

ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.915735525189067;
[fitresult, gof] = fit(ERKtip(:,1), GemPlusTip, ft, opts );

ft_int = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts_int = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_int.Display = 'Off';
opts_int.StartPoint = [0.915735525189067 0.792207329559554];
[fitresult_int, gof_int] = fit(ERKtip(:,1), GemPlusTip, ft_int, opts_int);

% a = plot(fitresult);
% set(a,'linestyle','-')
b = plot(fitresult_int);
set(b,'linestyle','-','color',[0 0 0 0.5],'linewidth',2)



xlabel('Average ERK Activity (distal 500 \mum, A.U.)');
ylabel('Frac. Prolif. (GEM+, distal 500 \mum)');
set(gca,'FontSize',16);

c = colorbar;
c= colorbar('Ticks', 0:0.5:10,'TickLabels',{'', '0', '','1','','2.5','','3.5','','5','','','','','','','','','','10'});
caxis([0 10])
c.FontSize = 16;
c.Label.String = 'PD03 Concentration (uM)';
%set(c, 'YDir', 'reverse' );
xlim([0,0.5]);
ylim([0,0.5]);

legend_name={strcat("R^2=",num2str(round(gof_int.rsquare,1)))};
legend([b], legend_name,...
        'Location','northwest','color','none','box','off');
