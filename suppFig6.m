%% Supplemental Figure 6

%% Supp. Fig. 6A
%load "analysis_mat_dbin150.mat" from "SuppFig6" folder

%step 1 - create analysis matrix that only contains 120 hpa (pre) and 144 hpa (24 hr post) *

analysis_mat = analysis_mat_dbin150;


analysis_mat_trim = [];
for i = 1:size(analysis_mat,2)
    if analysis_mat(i).hpa == 120 || analysis_mat(i).hpa == 132
        analysis_mat_trim = vertcat(analysis_mat_trim,analysis_mat(i));
    else
        continue
    end
    
end

%step 2 - create treatment mat *
treatmentMat = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1];

%step 3 - Plot difference at each position *

analysis_mat_combined = analysis_mat_trim;

%fishSet = [1,2,3,4,5,6,7,8,11,12,13,14,15,16,18];
preVals = linspace(1,15,8);
postVals = linspace(2,16,8);

%h4 = figure;
diffTreat = [];
diffCont = [];
xTreat = [];
xCont = [];

for i = 1:size(preVals,2)

    q=preVals(i);
    j=postVals(i);
    treatment = treatmentMat(q);

    KTRi = analysis_mat_combined(q).ktr;
    KTRf = analysis_mat_combined(j).ktr;
    centersi = analysis_mat_combined(q).ccrot;
    centersf = analysis_mat_combined(j).ccrot;
    centersXi = analysis_mat_combined(q).ccrot(:,1);
    centersXf = analysis_mat_combined(j).ccrot(:,1);
    
    maxXi = max(centersXi);
    minXi = min(centersXi);
    flipCentersXi = (centersXi-maxXi)*-1;
    flipMaxXi = max(flipCentersXi);
    flipMinXi = min(flipCentersXi);
    normCentersXi = (flipCentersXi-flipMinXi)/(flipMaxXi-flipMinXi);

    maxXf = max(centersXf);
    minXf = min(centersXf);
    flipCentersXf = (centersXf-maxXf)*-1;
    flipMaxXf = max(flipCentersXf);
    flipMinXf = min(flipCentersXf);
    normCentersXf = (flipCentersXf-flipMinXf)/(flipMaxXf-flipMinXf);
    
    dbins=5;
    
    [yMeani,xMeani,ySEMi,yNi] = bin_average(normCentersXi,KTRi,dbins,min(normCentersXi),max(normCentersXi));

    %set up plot label
    fish = analysis_mat_combined(q).fish;
    fish2 = analysis_mat_combined(j).fish;
    ray = analysis_mat_combined(q).ray;
    ray2 = analysis_mat_combined(j).ray;
    hpa = analysis_mat_combined(q).hpa;
    hpa2 = analysis_mat_combined(j).hpa;

    [yMeanf,xMeanf,ySEMf,yNf] = bin_average(normCentersXf,KTRf,dbins,min(normCentersXf),max(normCentersXf));

    diffHere = yMeanf-yMeani;

    if treatment == 1
        diffTreat = vertcat(diffTreat,diffHere');
        xTreat = vertcat(xTreat,xMeani');
    elseif treatment == 0
        diffCont = vertcat(diffCont,diffHere');
        xCont = vertcat(xCont,xMeani');
    end
    
end

%Prep data for box plot
plotData = horzcat(diffCont(:,1),diffTreat(:,1),diffCont(:,2),diffTreat(:,2),diffCont(:,3),diffTreat(:,3),diffCont(:,4),diffTreat(:,4),diffCont(:,5),diffTreat(:,5));


xMat = zeros(size(plotData,1),1);



h4 = figure;
for i = 1:size(diffCont,2)*2
    if i == 1 || i == 3 || i    == 5 || i == 7 || i == 9
        cl = [0 0 0 0.5];
    elseif i == 2 || i == 4 || i == 6 || i == 8 || i == 10
        cl = [1 0 0 0.5];
    end
    boxchart(xMat+i,plotData(:,i),'BoxFaceColor',cl,'markercolor',cl); hold on;
end

xticks([1 2 3 4 5 6 7 8 9 10])
xticklabels({'Amp. Site','Amp. Site','1/4 Regen.','1/4 Regen.','Middle','Middle','3/4 Regen.','3/4 Regen.','Tip','Tip'})

ylabel('Difference in Avg. Erk Activity (Pre-Post)')
xlabel('Position')
set(gca,'fontsize',16)

%% Supp. Figure 6B
% load "collectDataNew.mat" from "SuppFig6" folder

% load excel data

% Set up the Import Options and import the data - import updated Lreg measurement data (only measure converted ray)
opts = spreadsheetImportOptions("NumVariables", 17);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A3:Q14";

% Specify column names and types
opts.VariableNames = ["VarName1", "fish", "ray", "pre", "post", "Lamp", "Lamp_um", "convDay", "Lreg_conv", "phi_t0", "time_t0", "Lreg_t1", "phi_t1", "time_t1", "Lreg_t2", "phi_t2", "time_t2"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
LampMeasurements19Mar26 = readtable("/Users/ashleyrich/Documents/GitHub/RichErkGradients2025/SuppFig6/LampMeasurements_19Mar26_updateLreg.xlsx", opts, "UseExcel", false);

% Convert to output type
LampData = table2array(LampMeasurements19Mar26);

% Clear temporary variables
clear opts

% Fish list
fishList = [1 3 4 6 7 8 9 10 11 12 16 163];

% Use this to predict 27Mar26 - w/ adjustment - Ashley 9 April
% Plot actual velocity as function of space v. predicted velo as func of
% space
% Adjust prediction for discrepancy between actual Lreg and predicted Lreg
% - use single factor for all rays
% normalize x values by lamp

close all

%set up color map
%colorsHere = colormap(parula(size(fishList,2)));
tLevels = 80;
colorsHere = flipud(viridis(tLevels)); % Create the flipped matrix
colormap(colorsHere);                  % Apply it to the FIGURE
allPhis = LampData(:,10);
allFish = LampData(:,2);
[allPhiSort,sortIdx] = sort(allPhis);
minPhi = min(allPhis);
maxPhi = max(allPhis);



for qq = 1:size(fishList,2)
    fishNum = fishList(qq);

    for cc = 1:size(LampData,1)
        if LampData(cc,2) == fishNum
        LregFish0 = LampData(cc,9);
        LampFish0 = LampData(cc,7);
        time0 = LampData(cc,11);
        time1 = LampData(cc,14);
        phi0 = LregFish0./LampFish0;
        else
            continue
        end
    end

    baseFish = ['fish' num2str(fishNum)];

    veloHere = collectDataNew.(baseFish).veloHere;
    veloHere = horzcat(veloHere,0);
    centroidLoc = collectDataNew.(baseFish).centroidLoc_toPlot;
    centroidLoc = horzcat(centroidLoc,0);
    
        
    % parameters
    a1 = 0.05088; % dLdt/L vs total fractionGEM
    % a2 = 7.06473694; % G(u) vs E(u)
    a2 = 7.19766154; % G(u) vs E(u)
    % alpha = 2.4097; % G(u) vs E(u)
    alpha = 2.39807; % G(u) vs E(u)
    A0 = 0.35998; % A(phi)
    af = 1.2; % f(u)
    bf = 0.4; % f(u)
        
    phi = LregFish0/LampFish0;
    u=0:0.01:1;

    %get prediction v from equations
    for kk=1:length(phi)
        u1=u*phi(kk);
        v=a1*a2./(alpha+1)./af.*(A0.*(1-phi(kk))).^alpha.*phi.*((bf+af.*u).^(alpha+1)-bf^(alpha+1));
        vSD=a1*a2./(alpha+1)./af.*(A0.*(1-phi(kk))).^alpha.*phi.*((bf+af.*1).^(alpha+1)-bf^(alpha+1));
        v = v .* LampFish0; 
        vSD = vSD .* LampFish0;

        %alvin
        Cf = 1/af*1/(alpha+1)*((af+bf)^(alpha+1)-bf^(alpha+1)); % integral of f(u)^alpha
        gamma = a1*a2.*A0^alpha.*Cf; % lumped parameter
        % alvin equations
        dphidt = gamma.*(1-phi).^alpha.*phi;
        dphidt = a1*a2.*A0^alpha.*(1-phi).^alpha.*phi.*1/af*1/(alpha+1)*((af+bf)^(alpha+1)-bf^(alpha+1));
        dphidt = dphidt.*LampFish0;
        %v is velocity/lamp
        %v = v.*LampFish0;
        u1 = u1.*LampFish0;
        %v = dphidt;

        vTipPre = dphidt;

        adjFacNew = veloHere(1)./vTipPre;

        for z = 1:size(allPhiSort)
            if phi0 == allPhiSort(z)
                colorNum = z;
            else
                continue
            end
        end

        colorCode = round(phi*100);

        
        plot(u1./LampFish0,(v./LregFish0).*adjFacNew,'Color',colorsHere(colorCode,:),'LineWidth',1.5); hold on    
    end
    
    %overlay actual data on prediction data
    plot(centroidLoc./LampFish0,(veloHere./LregFish0),'o--','Color',colorsHere(colorCode,:),'LineWidth',2); hold on;

end

c = colorbar()
c.Label.String = 'Fraction Regenerated';
set( c, 'YDir', 'reverse' );

%title('Normalize X by Lamp')
ylabel('Velocity (\mum/hr) / L_{reg}')
xlabel('Position (\mum)')
set(gca,'fontsize',16)
ylim([-0.005,0.055])
xlim([0 0.8])

%saveas(gcf,['/Users/ashleyrich/Documents/DiTaliaLab/Experiments_inProg/2Mar26_mEOSconversions/centriods/plots_veloPred_veloAct_9Apr26_adjustRayFactor_normVlreg_yLim_updateColors.jpg'])

%% Supp. Fig. 6C
% Plot change in integrated density - use this one

close all
figure;

xVals = [ 2 2 2 2 2 2 1 1 1 1 1 1];
yVals = [33300.276 187217.631 412525.253 128668.503 122534.436 200587.695 677157.025 472295.684 488427.916 361966.942 338343.435 675032.14];

[h,p,ci,stats] = ttest2(yVals(7:12),yVals(1:6));

% yValsCyclo = [33300.276 187217.631 412525.253 128668.503 122534.436 200587.695];
% yValsWater = [677157.025 472295.684 488427.916 361966.942 338343.435 675032.14];
% 
% [h2,p2,ci2,stats2] = ttest2(yVals(7:12),yVals(1:6));

%bar(xVals,yVals)
boxchart(xVals,yVals); hold on
plot(xVals,yVals,'.','MarkerSize',20,'color','k')

xticks([1 2])
xticklabels({'+H2O','+Cyclo'})
%xtickangle(45)

xstat = [1,2];
ystat = [690000,690000];
plot(xstat,ystat,'-','color','k','linewidth',1)
txt = ['p = ' num2str(round(p,1,"significant"))];

text(1.4,670000,txt,'FontSize',16)

ylabel('Change in Int. Density')

set(gca,'fontsize',16)

%yline(200000,'--')
%xline(6.5,'--')

%% Supp. Fig. 6D

% build lagrangian matrix
% set initial conditions

time_dependent = 0;

t0 = 72;
phi0 = 0.115;
% t0 = 96;
% phi0 = 0.24;
% t0 = 24;
% phi0 = 0.025;

xsize =100;



% initialize
lagr_mat = [];
matHere = [];
for i = 0:xsize
    xi0 = i.*phi0./xsize;
    
    [t,y]=ode45(@(t,y)growthmodel_lagrangian(t,y,time_dependent),[t0:12:24*28],[phi0,xi0]); % y: 1st column is phi, 2nd is xi


    A0 = 0.35998; % A(phi)

    % f(u) 
    if time_dependent %from linear model Erk ~ 1 + phi_ + u:phi_ + phi_:dpa + u:phi_:dpa, averaged data
    %     af = (0.16677+t./24.*0.034583)./A0; % f(u)
    %     bf = (0.30788+t./24.*(-0.043484))./A0; % f(u)
        af = (0.17355+t./24.*0.032759)./A0; % f(u)
        bf = (0.32718+t./24.*(-0.023222))./A0; % f(u)
    else
        % af = 0.68059; % f(u)
        % bf = 0.66059; % f(u)
        af = 1.2; % f(u)
        bf = 0.4; % f(u)
    end


    ERKi = A0.*(1-y(:,1)).*(af.*(y(:,2)./y(:,1))+bf);
    
    
    % dydt = cell2mat(arrayfun(@(x,y)growthmodel_lagrangian(0,[x;y],time_dependent)',y(:,1),y(:,2),uni=0));
    dydt = cell2mat(arrayfun(@(tt,x,y)growthmodel_lagrangian(tt,[x;y],time_dependent)',t,y(:,1),y(:,2),uni=0));

    vi = dydt(:,2);
    
    
    % append
    matHere.t = t;
    matHere.xi = y(:,2);
    matHere.ERKi = ERKi;
    matHere.dERKidt =gradient(ERKi,t);
    matHere.ERKnormi = ERKi./ERKi(1);
    matHere.logERKi = log(ERKi);
    matHere.dlogERKidt = gradient(log(ERKi),t);
%     matHere.dlogERKidt2 = matHere.dERKidt./ERKi;
    matHere.vi = vi;
    
    matHere.i = i;
    matHere.ui0 = i/xsize;
    
    lagr_mat = [lagr_mat,matHere];

end
[lagr_mat.ui] = expand_cell(arrayfun(@(s)s.xi./lagr_mat(end).xi,lagr_mat,uni=0));

% compute dilusion via -rho*dv/dx
t_array = lagr_mat(1).t;
i_array = [lagr_mat.i];
x_t_by_i = cat(2,lagr_mat.xi);
v_t_by_i = cat(2,lagr_mat.vi);
div_t_by_i = nan(size(x_t_by_i));
for ti = 1:size(div_t_by_i,1)
    div_t_by_i(ti,:) = gradient(v_t_by_i(ti,:),x_t_by_i(ti,:));
end
[lagr_mat.div_i] = expand_cell(num2cell(div_t_by_i,1)); % dv/dx
[lagr_mat.dil_i] = expand_cell(arrayfun(@(s)-s.div_i.*s.ERKi,lagr_mat,uni=0)); % -rho*dv/dx



% plot predicted J(source) from tau(ligand decay time)
% ***builds: lagr_mat_bytime

tau = 24*4; % hours

% predict J 
[lagr_mat.J_i] = expand_cell(arrayfun(@(s)s.dERKidt+s.div_i.*s.ERKi+s.ERKi.*1/tau,lagr_mat,uni=0)); % drho/dt + rho*dv/dx +1/tau.*rho

% make matrix by time

lagr_mat_bytime = struct('t',num2cell(lagr_mat(1).t));
fieldnames_ = fieldnames(lagr_mat);
for f_i = 2:numel(fieldnames_) % skip t
    fieldHere = fieldnames_{f_i};
    if numel(lagr_mat(1).(fieldHere))==1
        [lagr_mat_bytime.(fieldHere)] = expand_cell(repmat({[lagr_mat.(fieldHere)]},size(lagr_mat_bytime)));
    else
        [lagr_mat_bytime.(fieldHere)] = expand_cell(num2cell(cat(2,lagr_mat.(fieldHere)),2));
    end
end
 
[lagr_mat_bytime.xi_tip] = expand_cell(arrayfun(@(s)s.xi(end)-s.xi,lagr_mat_bytime,uni=0));

% Generate Supplemental Figure 6D
% Predict Lfinal as a function of integral of J
[lagr_mat_bytime.int_J] = expand_cell(arrayfun(@(s)trapz(fliplr(s.xi_tip),fliplr(s.J_i)),lagr_mat_bytime,uni=0));
[lagr_mat_bytime.int_J_02] = expand_cell(arrayfun(@(s)trapz(fliplr(s.xi_tip(s.xi_tip<=0.2)),fliplr(s.J_i(s.xi_tip<=0.2))),lagr_mat_bytime,uni=0));


% plot
f = figure;

xname = 'int_J';
% xname = 'int_J_02';

i=0;
xdata = [];
ydata = [];
for Lamp = 1000:1000:5000
    i = i+1;
    xdata(i) = [lagr_mat_bytime(1).(xname)].*Lamp;
    ydata(i) = lagr_mat_bytime(end).xi(end).*Lamp;
end

plot(xdata,ydata,'-o',LineWidth=2)

config_plot(f);

xlim([0,7])
ylim([0,5000])

ylabel('L(28dpa)')
ylabel('Length_{Reg} (28 dpa)')
%xlabel('Integral of J')
xlabel('Integral of Production Souce')