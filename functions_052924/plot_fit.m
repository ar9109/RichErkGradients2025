function [p_fit] = plot_fit(fit_obj,fit_x,lgd,options)
arguments
    fit_obj
    fit_x
    lgd = [];
    options.binned = []; % binned average
    options.data = []; % un-binned data
    options.confidence = true % whether or not to show prediction interval
    options.visible = 'on'
    options.fit_color = 'k'
end
%PLOT_FIT Summary of this function goes here
%   Detailed explanation goes here
data = options.data;
binned = options.binned;
confidence = options.confidence;
fit_color = options.fit_color;

y = fit_obj(fit_x);

% plot
% f = figure('visible',options.visible);hold on;
hold on;
cm = lines(30);
% plot scatter
if ~isempty(data) && isempty(binned)
    p_data = plot(data.x,data.y,'o','Color',cm(1,:),'LineWidth',2);
elseif ~isempty(data) && ~isempty(binned)
    p_data = plot(data.x,data.y,'o','Color',cm(3,:),'LineWidth',2);
end
if ~isempty(binned)
    if isfield(binned,'sem')
        p_bin = errorbar(binned.x,binned.y,binned.sem,'-','color',cm(1,:),'LineWidth',2);
    else
        p_bin = plot(binned.x,binned.y,'o','color',cm(1,:),'LineWidth',2);
    end
end

% p_fit = plot(fit_x,y,'-','color',cm(2,:),'LineWidth',4);
p_fit = plot(fit_x,y,'-','color',fit_color,'LineWidth',4);

if confidence
    pred_bound = predint(fit_obj,fit_x,0.95,'function','on');
    p_bound = plot(fit_x,pred_bound,'--','LineWidth',2, 'Color',[0.5,0.5,0.5]);
end
leg.p = [p_fit];
leg.text = {};
% customize legend
if ~isempty(legend)
    if isstruct(lgd)
        leg.text = [leg.text, ['R^2 = ',num2str(round(lgd.rsquare,2))]];
    
    elseif isstring(lgd) || ischar(lgd)
        leg.text = [leg.text, [lgd,' = ',num2str(fit_obj.(lgd))]];
    
    else
        leg.text = [leg.text, ['Fit']];
    
    end
    
    if confidence
        leg.p = [leg.p, p_bound(1)];
        leg.text = [leg.text, '95% interval'];
    end
    if ~isempty(binned)
        leg.p = [leg.p, p_bin];
        leg.text = [leg.text, 'Binned average'];
    end
    
    legend(leg.p,...
        leg.text,...
        'Location','best','color','none','box','off')
end

hold off;
% config_plot(f);
end

