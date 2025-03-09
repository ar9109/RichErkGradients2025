function [f] = config_plot_flat(f)
arguments
    f

end
%CONFIG_PLOT Summary of this function goes here
%   Detailed explanation goes here

f.OuterPosition = [100,100,700,400];
set(f.CurrentAxes, 'fontsize', 24,'linewidth',3','Position' ,[0.1921    0.3118    0.6553    0.6]);


box off;

% set plot background color
set(f,'color','w');

% set legend
legend_f = findobj(f.Children,'Type','Legend');
if ~isempty(legend_f)
    set(legend_f,'color', 'none','location','best','box','off','fontsize', 20);
end

