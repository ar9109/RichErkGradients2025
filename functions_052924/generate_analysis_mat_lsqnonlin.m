function [analysis_mat_gem] = generate_analysis_mat_lsqnonlin(st_dir,pxsize,dbins,fish_modifier)
%GENERATE_ANALYSIS_MAT Summary of this function goes here
%   fish_modifier is used to distinguish fish from different expriments but
%   have the same fish number
arguments
    st_dir
    pxsize
    dbins = 50
    fish_modifier = 0
end

analysis_mat_gem = [];
for i =1:numel(st_dir)
    disp('------------');
    display([st_dir(i).folder filesep st_dir(i).name]);
    
    % create analysis_mat struct for each single experiments
    single_exp = [];
    load([st_dir(i).folder filesep st_dir(i).name]);
    
    % Load and remove NaN/Inf; invert gem to get nuc/cyt if needed
    cc=myScale.TGMM_hypo_eq_ch2.centers;
    ktr=myScale.TGMM_hypo_eq_ch2.ktr;
    %idxs=FUCCI>1; %All the FUCCI values that are not infinity
    idxs=~isnan(ktr) & ~isinf(ktr); %All the FUCCI values that are not NaN
    %idxs=~isinf(idxs);
    cc=cc(idxs,:);%
    ktr=ktr(idxs);%
    logktr = log(ktr);
    if isfield(myScale.TGMM_hypo_eq_ch2,'GEM')
        gem=myScale.TGMM_hypo_eq_ch2.GEM(idxs);
    elseif isfield(myScale.TGMM_hypo_eq_ch2,'greenCytNuc')
        gem=1./myScale.TGMM_hypo_eq_ch2.greenCytNuc(idxs);
    end

    %BoneRemoved = myScale.TGMM_hypo_eq_ch3.BoneRemoved;
    
    ampPointLabel = 'amputationPoint';
    endPointLabel = 'endPoint';
    
    ampPoint = myScale.TGMM_hypo_eq_ch2.(ampPointLabel);
    endPoint = myScale.TGMM_hypo_eq_ch2.(endPointLabel);
    % Rotate centers and amp/endPoints
    [ccrot,ampPointrot,endPointrot] = process_coords(cc,ampPoint,endPoint,pxsize);
    
    
    
    % binning
%     dbins = 50;
%     bins = linspace(0,(s.endPointrot(1)-s.ampPointrot(1)),21);
%     bins = 0:dbins:(endPointrot(1)-ampPointrot(1));
    bins = flip((endPointrot(1)-ampPointrot(1)):-dbins:0);
    
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
        averageKTR(bin) = mean(ktr(idxBin==bin));
        semKTR(bin) = std(ktr(idxBin==bin))/sqrt(sum(idxBin==bin));
        averagelogKTR(bin) = mean(logktr(idxBin==bin));
        semlogKTR(bin) = std(logktr(idxBin==bin))/sqrt(sum(idxBin==bin));
        varKTR(bin) = var(ktr(idxBin==bin));
        % generate weight matrix from varKTR
        varKTR_sc(idxBin==bin) = varKTR(bin);
        % GEM
        nucleiGEM(bin) = sum(idxBin==bin&gem>1.5);
        nucleiALL(bin) = sum(idxBin==bin);
        fractionGEM(bin) = nucleiGEM(bin)./nucleiALL(bin);
        if fractionGEM(bin)==0
            pseudo_fractionGEM(bin) = (nucleiGEM(bin)+1)./nucleiALL(bin); % prevent zero variance
        elseif fractionGEM(bin)==1
            pseudo_fractionGEM(bin) = (nucleiGEM(bin)-1)./nucleiALL(bin); % prevent zero variance
        else
            pseudo_fractionGEM(bin) = fractionGEM(bin);
        end
        varGEM(bin) = (pseudo_fractionGEM(bin).*(1-pseudo_fractionGEM(bin)))./nucleiALL(bin);
        % binvalue
        binvalue(bin) = (bins(bin)+bins(bin+1))/2;
    end

    % trim at amp points and end points for GEM
    trim_logical_GEM = ccrot(:,1)<=(endPointrot(1)-ampPointrot(1)) & ccrot(:,1)>=0; 
    % trim at amp points and end points
    trim_logical_KTR = ccrot(:,1)<=(endPointrot(1)-ampPointrot(1)) & ccrot(:,1)>=0 & ~isnan(idxBin) & varKTR_sc~=0; 



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
    ub = [0.01,Inf,(endPointrot(1)-ampPointrot(1))];

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
    [x_onerun,sse,res,~,~,~,J] = lsqnonlin(residual_fun,x,lb,ub,opts);


    fit_gem_heavi.a = x_onerun(1);
    fit_gem_heavi.b = x_onerun(2);
    fit_gem_heavi.c = x_onerun(3);
    fit_gem_heavi.res = res;
    fit_gem_heavi.sse = sse;
    fit_gem_heavi.J = J;


    %%%%% 2 fit data with weighted multistart lsqnonlin by N
    residual_fun = @(p)(fit_func(p,xdata)-ydata).*nucleiALL(~isnan(fractionGEM)).^(1/2);

    % multistart
    problem = createOptimProblem('lsqnonlin',...
    'objective',residual_fun,...
    'x0',p0,'ub',ub,'lb',lb,'options',opts);
    ms = MultiStart('UseParallel',true);
    [x,fval,eflag,output,manymins] = run(ms,problem,100);
%     ss = [manymins.X];
    [x_onerun,sse,res,~,~,~,J] = lsqnonlin(residual_fun,x,lb,ub,opts);
    
    
    fit_gem_heavi_weightN.a = x_onerun(1);
    fit_gem_heavi_weightN.b = x_onerun(2);
    fit_gem_heavi_weightN.c = x_onerun(3);
    fit_gem_heavi_weightN.res = res;
    fit_gem_heavi_weightN.sse = sse;
    fit_gem_heavi_weightN.J = J;



    %%%%% 3 fit data with weighted multistart lsqnonlin by binomial variance
    residual_fun = @(p)(fit_func(p,xdata)-ydata).*varGEM(~isnan(fractionGEM)).^(-1/2);

    % multistart
    problem = createOptimProblem('lsqnonlin',...
    'objective',residual_fun,...
    'x0',p0,'ub',ub,'lb',lb,'options',opts);
    ms = MultiStart('UseParallel',true);
    [x,fval,eflag,output,manymins] = run(ms,problem,100);
%     ss = [manymins.X];
    [x_onerun,sse,res,~,~,~,J] = lsqnonlin(residual_fun,x,lb,ub,opts);
    
    
    fit_gem_heavi_weightVar.a = x_onerun(1);
    fit_gem_heavi_weightVar.b = x_onerun(2);
    fit_gem_heavi_weightVar.c = x_onerun(3);
    fit_gem_heavi_weightVar.res = res;
    fit_gem_heavi_weightVar.sse = sse;
    fit_gem_heavi_weightVar.J = J;
    



    %%%%% fit linear to KTR
    % seed
    b0KTR = max(averageKTR);
    c0 = (endPointrot(1)-ampPointrot(1))/2;
    a0KTR = (b0KTR-averageKTR(end))/c0;    
    


    %%%%% fit KTR data with weighted linear least square
    ft = fittype({'-x','1'});
    foKTR = fitoptions('Method','LinearLeastSquares',...
        'Lower',[0,0],...
        'Weights', varKTR_sc(trim_logical_KTR).^(-1));

    xdata = ccrot(trim_logical_KTR,1);
    ydata = ktr(trim_logical_KTR);

    [fit_ktr_linear,fit_ktr_linear_gof] = fit(xdata,ydata,ft,foKTR);






    %%%%%
    % create analysis_mat single_experiment
    single_exp.name = st_dir(i).name;
    single_exp.fish = myScale.fish+fish_modifier;
    single_exp.ray = myScale.ray;
    single_exp.hpp = myScale.hpp;
    single_exp.hppTrue = myScale.hppTrue;
    single_exp.ktr = ktr;
    single_exp.gem = gem;
    single_exp.ccrot = ccrot;
    single_exp.ampPointrot = ampPointrot;
    single_exp.endPointrot = endPointrot;

    single_exp.trim_logical_GEM = trim_logical_GEM;
    single_exp.trim_logical_KTR = trim_logical_KTR;

    single_exp.fit_gem_heavi = fit_gem_heavi;
    single_exp.fit_gem_heavi_weightN = fit_gem_heavi_weightN;
    single_exp.fit_gem_heavi_weightVar = fit_gem_heavi_weightVar;
    single_exp.nucleiGEM = nucleiGEM;
    single_exp.nucleiALL = nucleiALL;
    single_exp.fractionGEM = fractionGEM;
    single_exp.varGEM = varGEM;
    
    single_exp.fit_ktr_linear = fit_ktr_linear;
    single_exp.fit_ktr_linear_gof = fit_ktr_linear_gof;
    single_exp.averageKTR = averageKTR;
    single_exp.semKTR = semKTR;
    single_exp.averagelogKTR = averagelogKTR;
    single_exp.semlogKTR = semlogKTR;
    single_exp.varKTR = varKTR;

    single_exp.binvalue = binvalue;

    analysis_mat_gem = [analysis_mat_gem single_exp];
end
end

