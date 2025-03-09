function [f] = config_plot(f,c)
arguments
    f
    c = [];
end
%CONFIG_PLOT Summary of this function goes here
%   Detailed explanation goes here

f.OuterPosition = [100,100,700,700];
set(f.CurrentAxes, 'fontsize', 24,'linewidth',3','Position' ,[0.2021    0.1618    0.6553    0.689]);


if ~isempty(c)
    f.OuterPosition = [100,100,900,700];
    set(f.CurrentAxes, 'fontsize', 24,'linewidth',3','Position' ,[0.1921    0.1618    0.5053    0.689]);
    c.LineWidth = 3;
    c.FontSize = 24;
end
box on;

% set plot background color
set(f,'color','w');

% set legend
legend_f = findobj(f.Children,'Type','Legend');
if ~isempty(legend_f)
    set(legend_f,'color', 'none','location','best','box','off','fontsize', 20);
end

