function [p_odr,c,fit_obj,fit_obj_gof,B] = plot_Gx_by_Ex(analysis_mat,color_res,r)
%PLOT_GX_BY_EX Summary of this function goes here
% field_suffix = '';
% field_suffix = '_u_timeave';
% field_suffix = '_u_20bins_timeave';
field_suffix = '';
field_suffix_ktr = '_res';

field_ktr = ['averageKTR' field_suffix_ktr];
field_gem = ['fractionGEM' field_suffix];

% field_gem = ['fractionGEM' '_u_timeresum'];


fit_obj=[];
fit_obj_gof=[];
phi = arrayfun(@(x)x.L_reg./x.L_amp,analysis_mat);
x_dim = isrow(analysis_mat(1).(field_ktr))+1; % cat in different dimension depending on whether data is row or column vector
y_dim = isrow(analysis_mat(1).(field_gem))+1;
bin_x = [];
bin_y = [];
% f = figure('visible','on');hold on;
cm = colormap(cool(round(color_res)+1)); % config color
color_unit = (1-min(phi))./color_res; % config color\
color_unit = (1-0)./color_res; % config color

range = 0:0.1:1;
% r = 10;


for i = numel(analysis_mat):-1:1
%     if phi(i)>range(r)&&phi(i)<=range(r+1)
    s = analysis_mat(i);
    sz = min([numel(s.(field_ktr)),numel(s.(field_gem))]);

    cl = cm(round((min(phi(i),1)-min(phi))./color_unit)+1,:); % config color
    cl = cm(round((min(phi(i),1)-0)./color_unit)+1,:); % config color

    x_plot = s.(field_ktr)(1:sz);
    y_plot = s.(field_gem)(1:sz);

%     x_plot = x_plot(s.binvalue>=s.lambda);
%     y_plot = y_plot(s.binvalue>=s.lambda);

    plot(x_plot,y_plot,'.','color',cl,'markersize',10);hold on;
    bin_x = cat(x_dim,bin_x,x_plot);
    bin_y = cat(y_dim,bin_y,y_plot);
%     end
end
    c = colorbar('Ticks',0:0.5:1,'TickLabels', round((linspace(min(phi),1,3)),2));
    c = colorbar('Ticks',0:0.5:1,'TickLabels', round((linspace(0,1,3)),2));

    c.Label.String = '\phi';
    c.LineWidth = 3;
    c.FontSize = 24;
%     c.Ruler.TickLabelFormat = '%.1e';

% plot binned average
bin_x = bin_x(bin_y>=0);
bin_y = bin_y(bin_y>=0);
bin_x = bin_x(:);
bin_y = bin_y(:);
[binned.y,binned.x,binned.sem]=bin_average(bin_x(bin_x>0.65&bin_x<1.2),bin_y(bin_x>0.65&bin_x<1.2),15);
% [binned.y,binned.x,binned.sem]=bin_average(bin_x(bin_x>0.6&bin_x<1.4),bin_y(bin_x>0.6&bin_x<1.4),15);

% plot binned average without fitting
hold on;
cl0 = lines(30);
% errorbar(binned.x-0.8,binned.y,binned.sem,'-o','color',cl0(1,:),'linewidth',2,'CapSize',10)

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


% [fit_obj,fit_obj_gof] = fit_lsqnonlin(xdata,ydata,fit_func_fit,p0,lb,ub,100);
% plot_fit(fit_obj,linspace(0,0.6,100),'d')

% python odr
myoutput = odr_scipy(xdata,ydata);
myoutput.pprint();
B = double(myoutput.beta);
x = linspace(0,0.6,100);
hold on;
p_odr = plot(x,B(1)*x.^B(2),'k-',LineWidth=4);
hold off;

% config
% xlim([0.6,1.4])
% xlim([0,0.5])
% ylim([0,0.8])
% ylim([0,2])
hold off;
xlabel('averageERK(u)')
ylabel('%GEM+(u)')
f.OuterPosition = [100,100,650,600];
set(gca, 'fontsize', 20,'linewidth',3','Position' ,[0.1269    0.1561    0.6053    0.7689]);
end

