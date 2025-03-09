function [f,c] = plot_colorcode(f,data_str,fieldname_x, fieldname_y,fieldname_c,options)
arguments
    f
    data_str
    fieldname_x
    fieldname_y
    fieldname_c
    options.color_res = 100;
    options.c_digit = -2;
    options.c_range = [];
    options.colorbar_label = fieldname_c;
    options.linewidth = 1.5;
    options.cm_func = @(x)flipud(cool(x));
    options.lineshape = '-';
    options.color_skip = 0;
end
color_res = options.color_res;
c_digit = options.c_digit;
c_range = options.c_range;
cm_func = options.cm_func;
color_skip = options.color_skip;

%PLOT_COLORCODE Summary of this function goes here
%   Detailed explanation goes here
if isstring(fieldname_c)||ischar(fieldname_c)
    color_code = [data_str.(fieldname_c)];
else
    color_code = arrayfun(fieldname_c,data_str);
end
if isstring(options.colorbar_label)||ischar(options.colorbar_label)
    colorbar_label = options.colorbar_label;
else
    colorbar_label = 'User defined fun';
end
% color_code = [data_str.(fieldname_c)];

% cm_matrix = cm_func(round(color_res));
% cm = colormap(cm_matrix); % config color
cm_matrix = cm_func(round(color_res+round(color_res.*color_skip)));
cm = colormap(cm_matrix((round(color_res.*color_skip)+1):end,:)); % config color


if ~isempty(c_range)&&~isinf(c_range(1))
    cmin = c_range(1);
else
    cmin = round(min(color_code),c_digit);
end
if ~isempty(c_range)&&~isinf(c_range(2))
    cmax = c_range(2);
else
    cmax = round(max(color_code),c_digit);
end
color_unit = (cmax-cmin)./(color_res-1); % config color

% hold on;
%%

% normby = data_str(1).(fieldname_y)(end);
for i = 1:numel(data_str)
        % create cropped myScale_merged struct
        s = data_str(i);
        if isstring(fieldname_x)||ischar(fieldname_x)
            xdata = s.(fieldname_x);
        else
            xdata = fieldname_x(s);
        end
        if isstring(fieldname_y)||ischar(fieldname_y)
            ydata = s.(fieldname_y);
        else % if fieldname is a function
            ydata = fieldname_y(s);
        end
        
        if numel(xdata)<=1
            continue
        end
        toberemoved = [xdata(:),ydata(:)];
        removed = remove_nan_inf(toberemoved);
        xdata = removed(:,1);
        ydata = removed(:,2);

        cl_idx = round((color_code(i)-cmin)./color_unit)+1;
        if cl_idx <= 0
            cl_idx = 1;
        elseif cl_idx > round(color_res)
            cl_idx = round(color_res);
        end
        cl = cm(cl_idx,:); % config color

        plot(xdata,ydata,options.lineshape,'color',[cl,0.6],'markersize',8,'LineWidth',options.linewidth);
        hold on;

end
%%
c = colorbar('Ticks',0:0.5:1,'TickLabels', round((linspace(cmin,cmax,3)),c_digit));
c.Label.String = colorbar_label;
hold off;


end

