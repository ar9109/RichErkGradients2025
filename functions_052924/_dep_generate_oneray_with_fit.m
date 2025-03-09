function [single_exp] = generate_oneray(single_exp,ktr,gem,ccrot,ampPointrot,endPointrot,dbins)
%GENERATE_ONERAY Summary of this function goes here
%   Detailed explanation goes here



% binning
%     dbins = 50;
%     bins = linspace(0,(s.endPointrot(1)-s.ampPointrot(1)),21);
bins = 0:dbins:(endPointrot(1)-ampPointrot(1));
%     bins = flip((endPointrot(1)-ampPointrot(1)):-dbins:0);

idxBin = discretize(ccrot(:,1),bins); % assigns each column to a particular bin
% Initialise
% KTR
varKTR_sc = idxBin;

averageKTR = nan(1,numel(bins)-1);
semKTR = nan(1,numel(bins)-1);
averagelogKTR = nan(1,numel(bins)-1);
semlogKTR = nan(1,numel(bins)-1);
varKTR = nan(1,numel(bins)-1);

% GEM
fractionGEM = nan(1,numel(bins)-1);
nucleiGEM = nan(1,numel(bins)-1);
nucleiALL = nan(1,numel(bins)-1);
pseudo_fractionGEM = nan(1,numel(bins)-1); % prevent zero variance
varGEM = nan(1,numel(bins)-1); % binomial variance

% binvalue
binvalue = nan(1,numel(bins)-1);

for bin = 1:(numel(bins)-1)
    % KTR
    if ~isempty(ktr)
        logktr = log(ktr);
        averageKTR(bin) = mean(ktr(idxBin==bin));
        semKTR(bin) = std(ktr(idxBin==bin))/sqrt(sum(idxBin==bin));
        averagelogKTR(bin) = mean(logktr(idxBin==bin));
        semlogKTR(bin) = std(logktr(idxBin==bin))/sqrt(sum(idxBin==bin));
        varKTR(bin) = var(ktr(idxBin==bin));
        % generate weight matrix from varKTR
        varKTR_sc(idxBin==bin) = varKTR(bin);

        % trim at amp points and end points
        trim_logical_KTR = ccrot(:,1)<=(endPointrot(1)-ampPointrot(1)) & ccrot(:,1)>=0 & ~isnan(idxBin) & varKTR_sc~=0; 
    end
    % GEM
    if ~isempty(gem)
        nucleiGEM(bin) = sum(idxBin==bin&gem>1.5);
        nucleiALL(bin) = sum(idxBin==bin);
        fractionGEM(bin) = nucleiGEM(bin)./nucleiALL(bin);
        if nucleiALL(bin)==1 && nucleiGEM(bin)==0
            pseudo_fractionGEM(bin) = (nucleiGEM(bin)+0.5)./nucleiALL(bin); % prevent zero variance
        elseif nucleiALL(bin)==1 && nucleiGEM(bin)==1
            pseudo_fractionGEM(bin) = (nucleiGEM(bin)-0.5)./nucleiALL(bin); % prevent zero variance
        % normal pseudo-count
        elseif fractionGEM(bin)==0
            pseudo_fractionGEM(bin) = (nucleiGEM(bin)+1)./nucleiALL(bin); % prevent zero variance
        elseif fractionGEM(bin)==1
            pseudo_fractionGEM(bin) = (nucleiGEM(bin)-1)./nucleiALL(bin); % prevent zero variance
        else
            pseudo_fractionGEM(bin) = fractionGEM(bin); 
        end
        varGEM(bin) = (pseudo_fractionGEM(bin).*(1-pseudo_fractionGEM(bin)))./nucleiALL(bin);
    end
    % binvalue
    binvalue(bin) = (bins(bin)+bins(bin+1))/2;
end

% trim at amp points and end points for GEM
trim_logical = ccrot(:,1)<=(endPointrot(1)-ampPointrot(1)) & ccrot(:,1)>=0; 
 
if false %numel(binvalue)>2 edited 022623
    
    if ~isempty(gem)
        %%%%% fit heaviside for GEM
        % seed
        b0GEM = max(fractionGEM);
        c0 = (endPointrot(1)-ampPointrot(1))/2;
        a0GEM = (b0GEM-fractionGEM(end))/c0;
       
        
        %%%%% 1 fit data with multistart lsqnonlin
        fit_func = @(p,x)(heaviside(p(3)-x).*(-p(1).*x+p(2))+heaviside(x-p(3)).*(-p(1).*p(3)+p(2)));
        xdata = binvalue(~isnan(fractionGEM));
        ydata = fractionGEM(~isnan(fractionGEM));
        residual_fun = @(p)(fit_func(p,xdata)-ydata);
        p0 = [a0GEM,b0GEM,c0];
        lb = [0,0,0];
        ub = [0.01,1,(endPointrot(1)-ampPointrot(1))];
        
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
        [x,fval,eflag,output,manymins] = run(ms,problem,100);
    %     ss = [manymins.X];
        
        % use fit
        fit_func_fit = @(a,b,c,x)(heaviside(c-x).*(-a.*x+b)+heaviside(x-c).*(-a.*c+b));
        ft = fittype(fit_func_fit);
        foGEM = fitoptions('Method','NonlinearLeastSquares',...
                       'Lower',lb,...
                       'Upper',ub,...
                       'StartPoint',x,...
                       'MaxIter', 1000,...
                       'TolFun', 1e-8,...
                       'TolX',1e-8, ...
                       'DiffMaxChange',Inf,...
                       'DiffMinChange',0);
        [fit_gem_heavi,fit_gem_heavi_gof] = fit(xdata(:),ydata(:),ft,foGEM);
    
    
    
        %%%%% 2 fit data with weighted multistart lsqnonlin by N
        residual_fun = @(p)(fit_func(p,xdata)-ydata).*nucleiALL(~isnan(fractionGEM)).^(1/2);
    
        % multistart
        problem = createOptimProblem('lsqnonlin',...
        'objective',residual_fun,...
        'x0',p0,'ub',ub,'lb',lb,'options',opts);
        ms = MultiStart('UseParallel',true);
        [x,fval,eflag,output,manymins] = run(ms,problem,100);
    %     ss = [manymins.X];
    
        % use fit
        fit_func_fit = @(a,b,c,x)(heaviside(c-x).*(-a.*x+b)+heaviside(x-c).*(-a.*c+b));
        ft = fittype(fit_func_fit);
        foGEM.startpoint = x;
        foGEM.weights = nucleiALL(~isnan(fractionGEM));
        [fit_gem_heavi_weightN,fit_gem_heavi_weightN_gof] = fit(xdata(:),ydata(:),ft,foGEM);
    
    
    
        %%%%% 3 fit data with weighted multistart lsqnonlin by binomial variance
        residual_fun = @(p)(fit_func(p,xdata)-ydata).*varGEM(~isnan(fractionGEM)).^(-1/2);
        % multistart
        problem = createOptimProblem('lsqnonlin',...
        'objective',residual_fun,...
        'x0',p0,'ub',ub,'lb',lb,'options',opts);
        ms = MultiStart('UseParallel',true);
        [x,fval,eflag,output,manymins] = run(ms,problem,100);
    %     ss = [manymins.X];
    
        % use fit
        fit_func_fit = @(a,b,c,x)(heaviside(c-x).*(-a.*x+b)+heaviside(x-c).*(-a.*c+b));
        ft = fittype(fit_func_fit);
        foGEM.startpoint = x;
        foGEM.weights = varGEM(~isnan(fractionGEM)).^(-1);
        [fit_gem_heavi_weightVar,fit_gem_heavi_weightVar_gof] = fit(xdata(:),ydata(:),ft,foGEM);
        
    end

    if ~isempty(ktr)
    
        %%%%% fit linear to KTR
        % seed
        b0KTR = max(averageKTR);
        c0 = (endPointrot(1)-ampPointrot(1))/2;
        a0KTR = (b0KTR-averageKTR(end))/c0;    
        
    
    
        %%%%% 1 fit KTR sc data with linear least square
        ft = fittype({'-x','1'});
        foKTR = fitoptions('Method','LinearLeastSquares',...
            'Lower',[-Inf,0],...
            'Weights', []);
    
        xdata = ccrot(trim_logical,1);
        ydata = ktr(trim_logical);
    
        [fit_ktr_linear,fit_ktr_linear_gof] = fit(xdata(:),ydata(:),ft,foKTR);
    
        %%%%% 2 fit KTR binned data with weighted linear least square
        ft = fittype({'-x','1'});
        foKTR = fitoptions('Method','LinearLeastSquares',...
            'Lower',[-Inf,0],...
            'Weights', [semKTR(~isnan(averageKTR)&semKTR~=0).^(-2)]);
    
        xdata = binvalue(~isnan(averageKTR)&semKTR~=0)';
        ydata = averageKTR(~isnan(averageKTR)&semKTR~=0)';
    
        [fit_ktr_linear_binned,fit_ktr_linear_binned_gof] = fit(xdata(:),ydata(:),ft,foKTR);

    end




    %%%%%
    % create analysis_mat single_experiment
    if ~isempty(ktr)
        single_exp.ktr = ktr;
        single_exp.logktr = logktr; 
        single_exp.trim_logical_KTR = trim_logical_KTR;
        single_exp.fit_ktr_linear = fit_ktr_linear;
        single_exp.fit_ktr_linear_gof = fit_ktr_linear_gof;
        single_exp.fit_ktr_linear_binned = fit_ktr_linear_binned;
        single_exp.fit_ktr_linear_binned_gof = fit_ktr_linear_binned_gof;
        single_exp.averageKTR = averageKTR;
        single_exp.semKTR = semKTR;
        single_exp.averagelogKTR = averagelogKTR;
        single_exp.semlogKTR = semlogKTR;
        single_exp.varKTR = varKTR;
    else
        single_exp.ktr = [];
        single_exp.logktr = []; 
        single_exp.trim_logical_KTR = [];
        single_exp.fit_ktr_linear = [];
        single_exp.fit_ktr_linear_gof = [];
        single_exp.fit_ktr_linear_binned = [];
        single_exp.fit_ktr_linear_binned_gof = [];
        single_exp.averageKTR = [];
        single_exp.semKTR = [];
        single_exp.averagelogKTR = [];
        single_exp.semlogKTR = [];
        single_exp.varKTR = [];
    end
    if ~isempty(gem)
        single_exp.gem = gem;
        single_exp.fit_gem_heavi = fit_gem_heavi;
        single_exp.fit_gem_heavi_gof = fit_gem_heavi_gof;
        single_exp.fit_gem_heavi_weightN = fit_gem_heavi_weightN;
        single_exp.fit_gem_heavi_weightN_gof = fit_gem_heavi_weightN_gof;
        single_exp.fit_gem_heavi_weightVar = fit_gem_heavi_weightVar;
        single_exp.fit_gem_heavi_weightVar_gof = fit_gem_heavi_weightVar_gof;
        single_exp.nucleiGEM = nucleiGEM;
        single_exp.nucleiALL = nucleiALL;
        single_exp.fractionGEM = fractionGEM;
        single_exp.varGEM = varGEM;
    else
        single_exp.gem = [];
        single_exp.fit_gem_heavi = [];
        single_exp.fit_gem_heavi_gof = [];
        single_exp.fit_gem_heavi_weightN = [];
        single_exp.fit_gem_heavi_weightN_gof = [];
        single_exp.fit_gem_heavi_weightVar = [];
        single_exp.fit_gem_heavi_weightVar_gof = [];
        single_exp.nucleiGEM = [];
        single_exp.nucleiALL = [];
        single_exp.fractionGEM = [];
        single_exp.varGEM = [];
    end


    single_exp.ccrot = ccrot;
    single_exp.ampPointrot = ampPointrot;
    single_exp.endPointrot = endPointrot;
    single_exp.L_reg = endPointrot(1)-ampPointrot(1);


    single_exp.trim_logical = trim_logical;
    single_exp.binvalue = binvalue;

else % in case there's only one/two bin; whe ray is too short

    if ~isempty(ktr)
        %%%%% fit linear to KTR
        % seed
        b0KTR = max(averageKTR);
        c0 = (endPointrot(1)-ampPointrot(1))/2;
        a0KTR = (b0KTR-averageKTR(end))/c0;    
        
    
    
        %%%%% fit KTR sc data with linear least square
        ft = fittype({'-x','1'});
        foKTR = fitoptions('Method','LinearLeastSquares',...
            'Lower',[0,0],...
            'Weights', []);
    
        xdata = ccrot(trim_logical,1);
        ydata = ktr(trim_logical);
    
%         [fit_ktr_linear,fit_ktr_linear_gof] = fit(xdata,ydata,ft,foKTR);
        fit_ktr_linear = [];
        fit_ktr_linear_gof = [];
    
        single_exp.ktr = ktr;
        single_exp.logktr = logktr;
        single_exp.trim_logical_KTR = [];
        single_exp.fit_ktr_linear = fit_ktr_linear;
        single_exp.fit_ktr_linear_gof = fit_ktr_linear_gof;
        single_exp.fit_ktr_linear_binned = [];
        single_exp.fit_ktr_linear_binned_gof = [];
        single_exp.averageKTR = averageKTR;
        single_exp.semKTR = semKTR;
        single_exp.averagelogKTR = averagelogKTR;
        single_exp.semlogKTR = semlogKTR;
        single_exp.varKTR = varKTR;
    else
        single_exp.ktr = [];
        single_exp.logktr = [];
        single_exp.trim_logical_KTR = [];
        single_exp.fit_ktr_linear = [];
        single_exp.fit_ktr_linear_gof = [];
        single_exp.fit_ktr_linear_binned = [];
        single_exp.fit_ktr_linear_binned_gof = [];
        single_exp.averageKTR = [];
        single_exp.semKTR = [];
        single_exp.averagelogKTR = [];
        single_exp.semlogKTR = [];
        single_exp.varKTR = [];
    end

    if ~isempty(gem)
        single_exp.gem = gem;
        single_exp.fit_gem_heavi = [];
        single_exp.fit_gem_heavi_gof = [];
        single_exp.fit_gem_heavi_weightN = [];
        single_exp.fit_gem_heavi_weightN_gof = [];
        single_exp.fit_gem_heavi_weightVar = [];
        single_exp.fit_gem_heavi_weightVar_gof = [];
        single_exp.nucleiGEM = nucleiGEM;
        single_exp.nucleiALL = nucleiALL;
        single_exp.fractionGEM = fractionGEM;
        single_exp.varGEM = varGEM;
    else
        single_exp.gem = [];
        single_exp.fit_gem_heavi = [];
        single_exp.fit_gem_heavi_gof = [];
        single_exp.fit_gem_heavi_weightN = [];
        single_exp.fit_gem_heavi_weightN_gof = [];
        single_exp.fit_gem_heavi_weightVar = [];
        single_exp.fit_gem_heavi_weightVar_gof = [];
        single_exp.nucleiGEM = [];
        single_exp.nucleiALL = [];
        single_exp.fractionGEM = [];
        single_exp.varGEM = [];
    end
    single_exp.ccrot = ccrot;
    single_exp.ampPointrot = ampPointrot;
    single_exp.endPointrot = endPointrot;
    single_exp.L_reg = endPointrot(1)-ampPointrot(1);


    single_exp.trim_logical = trim_logical;

    single_exp.binvalue = binvalue;


end

