function [fit_obj_out,fit_obj_gof,x] = fit_lsqnonlin(xdata,ydata,fit_func_fit,p0,lb,ub,seedno,weights,compact_output)
arguments
    xdata
    ydata
    fit_func_fit
    p0
    lb
    ub
    seedno = 100;
    weights = [];  % weights should be inverse of variance
    compact_output = false; % whether get one structure output
end
%FIT_LSQNONLIN Summary of this function goes here
%   Detailed explanation goes here
fit_obj_out = [];
fit_obj_gof = [];
% fit_func_lsqnonlin = @(p,x)fit_func_fit(p(1),p(2),p(3),p(4),x);
% residual_fun = @(p)(fit_func_lsqnonlin(p,xdata)-ydata);

if ~isempty(weights)
    residual_fun = @(p)(fit_func_lsqnonlin(p,xdata,fit_func_fit)-ydata).*weights.^(1/2);
else
    residual_fun = @(p)(fit_func_lsqnonlin(p,xdata,fit_func_fit)-ydata);
end


% tolerance for multistart
opts = optimoptions('lsqnonlin',...
    'MaxFunctionEvaluations',600,...
    'MaxIterations',600,...
    'FunctionTolerance',1e-7,...
    'OptimalityTolerance',1e-7,...
    'StepTolerance',1e-7);

% multistart
problem = createOptimProblem('lsqnonlin',...
'objective',residual_fun,...
'x0',p0,'ub',ub,'lb',lb,'options',opts);
ms = MultiStart('UseParallel',true);
[x,fval,eflag,output,manymins] = run(ms,problem,seedno);
%     ss = [manymins.X];

try
% use fit
% fit_func_fit = @(a,b,c,x)(heaviside(c-x).*(-a.*x+b)+heaviside(x-c).*(-a.*c+b));
% ft = fittype(fit_func_fit);
ft = fit_func_fit;
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',lb,...
               'Upper',ub,...
               'StartPoint',x,...
               'MaxIter', 1000,...
               'TolFun', 1e-8,...
               'TolX',1e-8, ...
               'DiffMaxChange',Inf,...
               'DiffMinChange',0);
if ~isempty(weights)
    fo.weights = weights;
end

[fit_obj,fit_obj_gof] = fit(xdata(:),ydata(:),ft,fo);

if compact_output
    fit_obj_out.fit_obj = fit_obj;
    fit_obj_out.gof = fit_obj_gof;
    fit_obj_out.lsqn_param = x;
else
    fit_obj_out = fit_obj;
end
catch
    disp("!!!!!!!!!!!!!!!!!!Error caught when using the fit function!!!!!!!!!!!!!!!!!!!")
    fit_obj_out.lsqn_param = x;
end

end

%%
function results = fit_func_lsqnonlin(p,x,fit_func_fit)
p_cell = num2cell(p);

results = fit_func_fit(p_cell{:},x);



end