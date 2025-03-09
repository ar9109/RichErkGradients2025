function [c,fit_obj,fit_obj_gof] = plot_Gx_by_Ex_noODR(analysis_mat,color_res,r)
%PLOT_GX_BY_EX Summary of this function goes here
% field_suffix = '';
field_suffix = '_u_timeave';
% field_suffix = '_u_20bins_timeave';

field_ktr = ['averageKTR' field_suffix];
field_gem = ['fractionGEM' field_suffix];

% field_gem = ['fractionGEM' '_u_timeresum'];


fit_obj=[];
fit_obj_gof=[];
phi = arrayfun(@(x)x.L_reg./x.L_amp,analysis_mat);
Lamp = arrayfun(@(x)x.L_amp,analysis_mat);
x_dim = isrow(analysis_mat(1).(field_ktr))+1; % cat in different dimension depending on whether data is row or column vector
y_dim = isrow(analysis_mat(1).(field_gem))+1;
bin_x = [];
bin_y = [];
% f = figure('visible','on');hold on;
%%%cm = colormap(cool(round(color_res)+1)); % config color
%%%color_unit = (1-min(phi))./color_res; % config color\
%%%color_unit = (1-0)./color_res; % config color

tLevels = 101;
color_map=colormap(cbrewer2('Reds',tLevels));
lampMinHere = 706;
lampMaxHere = 4053;



range = 0:0.1:1;
% r = 10;

bin_x_small =  [];
bin_y_small = [];
bin_x_med = [];
bin_y_med = [];
bin_x_large = [];
bin_y_large = [];

for i = numel(analysis_mat):-1:1
%     if phi(i)>range(r)&&phi(i)<=range(r+1)
    s = analysis_mat(i);
    sz = min([numel(s.(field_ktr)),numel(s.(field_gem))]);

    %%%cl = cm(round((min(phi(i),1)-min(phi))./color_unit)+1,:); % config color
    %%%cl = cm(round((min(phi(i),1)-0)./color_unit)+1,:); % config color
    cl = ((Lamp(i)-700)./3353).*100+1;

    lampHere = s.L_amp;

    x_plot = s.(field_ktr)(1:sz);
    y_plot = s.(field_gem)(1:sz);

    if lampHere < 1575
        bin_x_small = vertcat(bin_x_small,x_plot);
        bin_y_small = vertcat(bin_y_small,y_plot);
    elseif lampHere < 3089
        bin_x_med = vertcat(bin_x_med,x_plot);
        bin_y_med = vertcat(bin_y_med,y_plot);
    else
        bin_x_large = vertcat(bin_x_large,x_plot);
        bin_y_large = vertcat(bin_y_large,y_plot);
    end


%     x_plot = x_plot(s.binvalue>=s.lambda);
%     y_plot = y_plot(s.binvalue>=s.lambda);

    colorPlot = color_map((round(cl)),:);
    colorPlot = horzcat(colorPlot,0.5);

    plot(x_plot-0.8,y_plot,'.','color',colorPlot,'markersize',10);hold on;
    bin_x = cat(x_dim,bin_x,x_plot);
    bin_y = cat(y_dim,bin_y,y_plot);
%     end
end
% % % % %     c = colorbar('Ticks',0:0.5:1,'TickLabels', round((linspace(min(phi),1,3)),2));
% % % % %     c = colorbar('Ticks',0:0.5:1,'TickLabels', round((linspace(0,1,3)),2));
% % % % % 
% % % % %     c.Label.String = '\phi';
% % % % %     c.LineWidth = 3;
% % % % %     c.FontSize = 24;
% % % % % %     c.Ruler.TickLabelFormat = '%.1e';

color_map=colormap(cbrewer2('Reds',tLevels));
c = colorbar;
%set(g,'TickLabels',{'0','811','1621','2432','3242','4053'});%,'2740','3080','3420','3760','4100'})
set(c,'TickLabels',{'700','1380','2060','2740','3420','4100'});
c.Label.String = 'Length Amputated (um)';

% plot binned average
bin_x = bin_x(bin_y>=0);
bin_y = bin_y(bin_y>=0);
bin_x = bin_x(:);
bin_y = bin_y(:);

bin_x_small = bin_x_small(bin_y_small>=0);
bin_y_small = bin_y_small(bin_y_small>=0);
bin_x_small = bin_x_small(:);
bin_y_small = bin_y_small(:);

bin_x_med = bin_x_med(bin_y_med>=0);
bin_y_med = bin_y_med(bin_y_med>=0);
bin_x_med = bin_x_med(:);
bin_y_med = bin_y_med(:);

bin_x_large = bin_x_large(bin_y_large>=0);
bin_y_large = bin_y_large(bin_y_large>=0);
bin_x_large = bin_x_large(:);
bin_y_large = bin_y_large(:);

%bin_x-0.8; %ashley added
%[binned.y,binned.x,binned.sem]=bin_average(bin_x(bin_x>0.65&bin_x<1.8),bin_y(bin_x>0.65&bin_x<1.8),15);
%%%[binned.y,binned.x,binned.std,binned.n,binned.sem]=bin_average(bin_x(bin_x>0.65&bin_x<1.4),bin_y(bin_x>0.65&bin_x<1.4),15);
%%%%%%[binned.y,binned.x,binned.sem]=bin_average(bin_x(bin_x>0.65&bin_x<1.25),bin_y(bin_x>0.65&bin_x<1.25),15);
%%%[binned.y,binned.x,binned.sem]=bin_average(bin_x(bin_x>0.65&bin_x<1.2),bin_y(bin_x>0.65&bin_x<1.2),15);
 %[binned.y,binned.x,binned.sem]=bin_average(bin_x(bin_x>0.6&bin_x<1.4),bin_y(bin_x>0.6&bin_x<1.4),15);

[binned.y,binned.x,binned.sem]=bin_average(bin_x(bin_x>0.6&bin_x<1.25),bin_y(bin_x>0.6&bin_x<1.25),15,0.6,1.25);
[binned.ySmall,binned.xSmall,binned.semSmall]=bin_average(bin_x_small(bin_x_small>0.6&bin_x_small<1.25),bin_y_small(bin_x_small>0.6&bin_x_small<1.25),15,0.6,1.25);
[binned.yMed,binned.xMed,binned.semMed]=bin_average(bin_x_med(bin_x_med>0.6&bin_x_med<1.25),bin_y_med(bin_x_med>0.6&bin_x_med<1.25),15,0.6,1.25);
[binned.yLar,binned.xLar,binned.semLar]=bin_average(bin_x_large(bin_x_large>0.6&bin_x_large<1.25),bin_y_large(bin_x_large>0.6&bin_x_large<1.25),15,0.6,1.25);

% plot binned average without fitting
hold on;
cl0 = lines(30);
%errorbar(binned.x-0.8,binned.y,binned.sem,'-o','color',cl0(1,:),'linewidth',2,'CapSize',10)
%errorbar(binned.x-0.8,binned.y,binned.sem,'-o','color','k','linewidth',2,'CapSize',10)

errorbar(binned.xSmall-0.8,binned.ySmall,binned.semSmall,'-o','color',[1 0.5 0.5],'linewidth',3,'CapSize',10); hold on;
errorbar(binned.xMed-0.8,binned.yMed,binned.semMed,'-o','color',[1 0 0],'linewidth',3,'CapSize',10)
errorbar(binned.xLar-0.8,binned.yLar,binned.semLar,'-o','color',[0.5 0 0],'linewidth',3,'CapSize',10)
errorbar(binned.x-0.8,binned.y,binned.sem,'-o','color','k','linewidth',3,'CapSize',10); hold on;

% fit with offset

% fit_func_fit = @(a,b,c,d,x)(heaviside(c-x).*(b)+heaviside(x-c).*(a.*(x-c).^d+b));
% p0 = [1,0.01,0.8,1];
% lb = [0,0,0,0];
% ub = [Inf,1,1.2,Inf];

% fit no offset

% fit_func_fit = @(a,c,d,x)(heaviside(x-c).*(a.*(x-c).^d));
% p0 = [1,0.8,1];
% lb = [0,0,0];
% ub = [Inf,1.2,Inf];

fit_func_fit = @(a,d,x)(a.*(x).^d);
p0 = [1,1];
lb = [0,0];
ub = [Inf,Inf];

% fit and plot

% xdata = adjust_neg(binned.x-0.8);
% ydata = binned.y;

xdata = adjust_neg(remove_nan(bin_x(bin_x<=1.3))-0.8);
ydata = remove_nan(bin_y(bin_x<=1.3));


%%% [fit_obj,fit_obj_gof] = fit_lsqnonlin(xdata,ydata,fit_func_fit,p0,lb,ub,100); %%%
%%% plot_fit(fit_obj,linspace(0,0.6,100),'d') %%%

% % python odr
% myoutput = odr_scipy(xdata,ydata);
% myoutput.pprint();
% B = double(myoutput.beta);
% x = linspace(0,0.6,100);
% hold on;
% p_odr = plot(x,B(1)*x.^B(2),'k-',LineWidth=4);
% hold off;

% config

xlim([0,0.5])
ylim([0,0.8])

hold off;
xlabel('Binned Average ERK (AU)')
ylabel('Binned Fraction Cycling')
% f.OuterPosition = [100,100,650,600];
set(gca, 'fontsize', 20,'linewidth',3')%,'Position' ,[0.1269    0.1561    0.6053    0.7689]);
end

