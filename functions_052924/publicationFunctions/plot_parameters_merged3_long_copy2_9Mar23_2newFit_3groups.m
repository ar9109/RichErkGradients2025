%function [f,fit_all] = plot_parameters_merged(fx,fy,para_mat,Xlabel,Ylabel,Title,ft)
function [f,fit_all,fit_l,fit_m] = plot_parameters_merged(fx,fy,para_mat,Xlabel,Ylabel,Title,ft)
%PLOT_PARAMETERS_MERGED Summary of this function goes here
%   Detailed explanation goes here
arguments
    fx
    fy
    para_mat
    Xlabel = 'x'
    Ylabel = 'y'
    Title = []
    ft = fittype({'x','1'})
    
end
L_amp = [para_mat.L_amp];

tLevels = 101;
color_map=colormap(cbrewer2('Reds',tLevels));
lampMinHere = 706;
lampMaxHere = 4053;

%f = figure('visible','on');
%%%cm = colormap(flipud(cool(round(max(L_amp))-round(min(L_amp))+1)));%AR added flipud
fit_l=[];
fit_m=[];
fit_s=[];
fit_all=[];
collectLamp_l = [];
collectLamp_m = [];
collectLamp_s = [];
collect_para_mat_l=[];
collect_para_mat_m=[];
collect_para_mat_s=[];
for i = 1:numel(para_mat)
    s = para_mat(i);

    if s.L_amp<1575
        %cl = 'b';
        fit_s=[fit_s;fx(s)' fy(s)'];
        num2gen = size(fx(s),2);
        Lamp_s_local = zeros(num2gen,1);
        Lamp_s_local(Lamp_s_local==0)=s.L_amp;
        collectLamp_s = [collectLamp_s;Lamp_s_local];
        collect_para_mat_s=[collect_para_mat_s;s];
    elseif s.L_amp<3089
        %cl = 'm';
        fit_m=[fit_m;fx(s)',fy(s)'];
        num2gen = size(fx(s),2);
        Lamp_m_local = zeros(num2gen,1);
        Lamp_m_local(Lamp_m_local==0)=s.L_amp;
        collectLamp_m = [collectLamp_m;Lamp_m_local];
        collect_para_mat_m=[collect_para_mat_m;s];
    else 
        %cl = 'm';
        fit_l=[fit_l;fx(s)',fy(s)'];
        num2gen = size(fx(s),2);
        Lamp_l_local = zeros(num2gen,1);
        Lamp_l_local(Lamp_l_local==0)=s.L_amp;
        collectLamp_l = [collectLamp_l;Lamp_l_local];
        collect_para_mat_l=[collect_para_mat_l;s];
    end
end

    idxsL=fit_l(:,2)>0;
    idxsM=fit_m(:,2)>0;
    idxsS=fit_s(:,2)>0;
    
    fit_l = fit_l(idxsL,:);
    fit_m = fit_m(idxsM,:);
    fit_s = fit_s(idxsS,:);
    
    collectLamp_l = collectLamp_l(idxsL,:);
    collectLamp_m = collectLamp_m(idxsM,:);
    collectLamp_s = collectLamp_s(idxsS,:);

    fit_all2 = vertcat(fit_l,fit_m, fit_s);

    for zz = 1:size(fit_l)
        colorCode = ((collectLamp_l(zz)-706)./3353).*100+1; %%% 699 will need to be edited 
        plot(fit_l(zz,1),fit_l(zz,2),'.','color',color_map((round(colorCode)),:),'MarkerSize',20);hold on;
    end

    for qqq = 1:size(fit_m)
        colorCode = ((collectLamp_m(qqq)-706)./3353).*100+1; %%% 699 will need to be edited 
        %colorCode
        plot(fit_m(qqq,1),fit_m(qqq,2),'.','color',color_map(round(colorCode),:),'MarkerSize',20);hold on;
    end

    for xxx = 1:size(fit_s)
        colorCode = ((collectLamp_s(xxx)-706)./3353).*100+1; %%% 699 will need to be edited 
        %colorCode
        plot(fit_s(xxx,1),fit_s(xxx,2),'.','color',color_map(round(colorCode),:),'MarkerSize',20);hold on;
    end

if ~isempty(ft)

    hold on

    fitxAll = fit_all2(:,1);
    fityAll = fit_all2(:,2);

    fit_l_x = fit_l(:,1);
    fit_l_y = fit_l(:,2);

    fit_m_x = fit_m(:,1);
    fit_m_y = fit_m(:,2);

    fit_s_x = fit_s(:,1);
    fit_s_y = fit_s(:,2);

    % Set up fittype and options. - Commented out 15Jan22
    %%ft = fittype( 'exp1' );
    %%opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %%opts.Display = 'Off';
    %%opts.StartPoint = [0.653167833208689 -0.0119972767220434];
    
    % Set up fittype and options.
    
    ft = fittype( 'a*x^b*exp(-c*x)', 'independent', 'x', 'dependent', 'y' );
    %1) ft = fittype( 'a*(exp(-b*x)-exp(-c*x))', 'independent', 'x', 'dependent', 'y' );
    %3) ft = fittype('a*exp(-(log(x)+b).^2/c)', 'independent', 'x', 'dependent', 'y' );
    %ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %opts = fitoptions( 'Method', 'SmoothingSpline' );
    %opts.SmoothingParam = 0.99;
    opts.Display = 'Off';
    opts.StartPoint = [0.915735525189067 0.792207329559554 0.959492426392903];

    % Fit model to data.
    [fitresult, gof] = fit(fitxAll, fityAll, ft, opts );
    [fitresult_l, gof_l] = fit(fit_l_x, fit_l_y, ft, opts );
    [fitresult_m, gof_m] = fit(fit_m_x, fit_m_y, ft, opts );
    [fitresult_s, gof_s] = fit(fit_s_x, fit_s_y, ft, opts );
    %[fitresult, gof] = fit_lsqnonlin(fitxAll, fityAll, ft, opts );
    
    %[fitresult, gof] = fit(fitxAll, fityAll, ft, opts);

    % Plot fit with data.
    fit_plot = plot(fitresult);
    set(fit_plot,'linewidth',2,'color','k');

    fit_plot_l = plot(fitresult_l);
    set(fit_plot_l,'linewidth',2,'color',[0.5 0 0]);

    fit_plot_m = plot(fitresult_m);
    set(fit_plot_m,'linewidth',2,'color',[1 0 0]);

    fit_plot_s = plot(fitresult_s);
    set(fit_plot_s,'linewidth',2,'color',[1 0.5 0.5]);
end

% axes spec
% xlim([0 2000]);
%ylim([0 0.5]);
xlabel(Xlabel,'Interpreter', 'none');
ylabel(Ylabel,'Interpreter', 'none');
title(Title,'Interpreter', 'none');
% % % c = colorbar('Ticks',0:0.5:1,'TickLabels', round((linspace(min(L_amp),max(L_amp),3)),-2));
c.Label.String = 'Length Amputated (\mum)';
c.LineWidth = 3;
c.FontSize = 20;

ax = gca;
ax.YAxis.Exponent = 0;

g = colorbar;
%set(g,'TickLabels',{'0','811','1621','2432','3242','4053'});%,'2740','3080','3420','3760','4100'})
set(g,'TickLabels',{'700','1380','2060','2740','3420','4100'});
g.Label.String = 'Length Amputated (um)';

f.OuterPosition = [1000,800,650,600];
set(gca, 'fontsize', 20,'linewidth',3);%,'Position' ,[0.1269    0.1561    0.6053    0.7689]);
%ylim([0 0.5])
if ~isempty(ft)
    %legend([fit_plot], {['R^2 = ' num2str(gof2.rsquare)]},'Location','best','color','none','box','off');
% % % % %     legend_name={strcat("R^2 = ",num2str(gof.rsquare)),...
% % % % %         strcat("R^2 = ", num2str(gof_mag.rsquare))};
legend_name={strcat("R^2 = ",num2str(gof.rsquare)),...
    strcat("R^2 = ",num2str(gof_l.rsquare)),...
    strcat("R^2 = ",num2str(gof_m.rsquare)),...
    strcat("R^2 = ",num2str(gof_s.rsquare))};
    %legend([fit_plot,fit_plot_mag], legend_name,...
        %'Location','best','color','none','box','off');
    legend([fit_plot,fit_plot_l,fit_plot_m,fit_plot_s], legend_name,...
        'Location','best','color','none','box','off');
    %legend([fit_plot,fit_plot_mag],["1","2"])
end

%hold off;
end

