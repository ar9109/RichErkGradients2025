function [ss,ub,a1,b1,a2,b2] = sumsqr_pred_lambda(analysis_mat,phi,lambda,E0,ET)
%SUMSQR_PRED_LAMBDA Summary of this function goes here
%   Detailed explanation goes here


% bin average for single ray-time
averageKTR_u_cell = arrayfun(@(s)bin_average(s.ccrot(s.trim_logical,1),...
    adjust_neg(s.ktr(s.trim_logical)-E0),...
    10,0,s.L_reg),...
    analysis_mat,'UniformOutput',false);

mean_ktr_bin10_offsetE0 = cellfun(@nanmean,averageKTR_u_cell);




%
%
%
%
% <Erk> vs phi coeff --- a1(1-phi)+b1

% bin average
data.x = phi;
data.y = mean_ktr_bin10_offsetE0;
[binned.y,binned.x,binned.sem] = bin_average(data.x,data.y,15);
fit_x = linspace(0,1,100);


% fit
xdata = binned.x(~isnan(binned.y)&~isnan(binned.sem));
ydata = binned.y(~isnan(binned.y)&~isnan(binned.sem));
ft = fittype({'-x','1'});
foKTR = fitoptions('Method','LinearLeastSquares',...
    'Lower',[0,0],...
    'Weights', [binned.sem(~isnan(binned.x)&~isnan(binned.sem)).^(-2)]);



[fit_obj1,fit_gof1] = fit(xdata(:),ydata(:),ft,foKTR);


a1 = fit_obj1.a;
b1 = fit_obj1.b - fit_obj1.a;
% a1 = 0.32872;
% b1 = adjust_neg(0.79655-E0);



%
%
%






% normalized Erk vs u coeff --- b2-a2u
% normalize, e = Erk(u), m = mean(Erk)
norm_averageKTR_u_cell = arrayfun(@(e,m)(e{:})./(m),averageKTR_u_cell,...
    mean_ktr_bin10_offsetE0,'UniformOutput',false);
norm_averageKTR_u = cell2mat(norm_averageKTR_u_cell);

% fit not-timeaveraged data with a line
% bin average
data.x = repmat([0.05:0.1:0.95]',[1,size(norm_averageKTR_u,2)]);
data.y = norm_averageKTR_u;
[binned.y,binned.x,binned.sem] = bin_average(data.x(:),data.y(:),10,0,1);

% fit to get coeff
xdata = data.x(~isnan(data.y));
ydata = data.y(~isnan(data.y));
ft = fittype({'-x','1'});
foKTR = fitoptions('Method','LinearLeastSquares',...
    'Lower',[0,0],...
    'Weights', []);
[fit_obj2,fit_gof2] = fit(xdata(:),ydata(:),ft,foKTR);

a2 = fit_obj2.a;
b2 = fit_obj2.b;


%
%
%
%
%
% calculate sum of square
lambda_pred = (b2-(ET-E0)./(a1.*(1-phi)+b1))./a2;
ss = sumsqr(adjust_neg(lambda_pred)-lambda);
% ss = sumsqr(lambda_pred(lambda_pred>=0)-lambda(lambda_pred>=0));
ub = b2*(a1.*(1-min(phi))+b1)+E0;

end

