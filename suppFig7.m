%% Supp. Fig. 7A
% Plot GFP domain length / Lamp v. Lamp - here
chToResize = [1]; %reference channel to resize 

%paths.objFolder is located in "SuppFig7" folder under the name
%"objects_old_22Apr24_new_copy_no0"
%Please adjust the paths.masterFolder and paths.objFolder lines (13 & 14)
%to where you store the data


paths=[];

paths.masterFolder='/Users/ashleyrich/Documents/GitHub/RichErkGradients2025/SuppFig7/'; %folder where data is stored
paths.objFolder= [paths.masterFolder 'objects_old_22Apr24_new_copy_no0/']; %objects folder

fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='72';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa*' '.mat']);

opts=[];
%opts.resizeFactor = 0.25; %resizing factor

optsSave=[];
optsSave.compress='lzw';
optsSave.overwrite=true;

h = figure;

lampCollect = [];
lGFPcollect = [];

for i = 1:size(st_dir)
    display([paths.objFolder filesep st_dir(i).name]);
    load([paths.objFolder filesep st_dir(i).name]);
    fish = myScale.fish;
    ray = myScale.ray;
    hpa = myScale.hpp;
    lGFP = myScale.GFPlength;
    lamp = myScale.lamp;

    if hpa > 70
        if fish < 5
            cl = 'k';
            plot(lamp,lGFP,'.','MarkerSize',15,'color',cl); hold on;
            lampCollect = vertcat(lampCollect,lamp);
            lGFPcollect = vertcat(lGFPcollect,lGFP);
%         elseif fish <20
%             cl = 'r';
        elseif fish > 20
            cl = 'k';
            plot(lamp,lGFP,'.','MarkerSize',15,'color',cl); hold on;
            lampCollect = vertcat(lampCollect,lamp);
            lGFPcollect = vertcat(lGFPcollect,lGFP);
        end

        %plot(lamp,lGFP./lamp,'.','MarkerSize',15,'color',cl); hold on;
    else
        continue
    end
end

% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit( lampCollect, lGFPcollect, ft );
% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
%h = 
f = plot(fitresult);
set(f,'linewidth',2,'color','k')
%legend(f, 'R^{2}', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );

legend_name={strcat("R^2 = ",num2str(round(gof.rsquare,1)))};

legend(f, legend_name,...
        'Location','best','color','none','box','off');

xlabel('Length Amputated (\mum)');
ylabel({'Length Fgf20a Expression Domain (\mum)'})
%ylabel('Length Fgf20a Expression Domain (\mum) / Length Amputated (\mum)');
set(gca,'fontsize',16);

%% Supp. Fig. 7B

%load "dataLOF2.mat" and "dataWT2.mat" from "SuppFig7" folder

close all
figure;

wtPlot = dataWT(:,7)./dataWT(:,5);
lofPlot = dataLOF(:,7)./dataLOF(:,5);

x_wt = linspace(0.75,1.25,size(wtPlot,1));
x_lof = linspace(1.75,2.25,size(lofPlot,1));

stdWT = std(wtPlot);
varWT = stdWT^2;
stdLOF = std(lofPlot);
varLOF = stdLOF^2;

x1 = 1;
x_wt_box = repelem(x1,size(wtPlot,1));

x2 = 2;
x_lof_box = repelem(x2,size(lofPlot,1));


a = boxchart(x_wt_box,wtPlot); hold on;
plot(x_wt,wtPlot,'.','MarkerSize',20,'Color','k');
a.BoxFaceColor = 'k';
a.BoxFaceAlpha = 0;
a.MarkerStyle = 'none';

b = boxchart(x_lof_box,lofPlot); hold on;
plot(x_lof,lofPlot,'.','MarkerSize',20,'Color','r');
b.BoxFaceColor = 'r';
b.BoxFaceAlpha = 0;
b.MarkerStyle = 'none';

[h,p,ci,stats] = ttest2(wtPlot,lofPlot)

txt1 = strcat('s^2 = ',num2str(round(varWT,1,"significant")));
text(0.85,1.25,txt1,'fontsize',16)

txt2 = strcat('s^2 = ',num2str(round(varLOF,1,"significant")));
text(1.90,1.75,txt2,'fontsize',16)

ylim([0,2]);

xticks([1 2])
xticklabels({'Wildtype','\it{Longfin}'})

%yticks([0 1.0 2.0 3.0 4.0])

ylabel({'Fraction Regen.'})

set(gca,'fontsize',16)

%% Supp. Fig. 7C

%load "dataLOF2.mat" and "dataWT2.mat" from "SuppFig7" folder
close all
figure;

wtPlot = dataWT(:,7)./dataWT(:,4);
lofPlot = dataLOF(:,7)./dataLOF(:,4);

x_wt = linspace(0.75,1.25,size(wtPlot,1));
x_lof = linspace(1.75,2.25,size(lofPlot,1));

x1 = 1;
x_wt_box = repelem(x1,size(wtPlot,1));

x2 = 2;
x_lof_box = repelem(x2,size(lofPlot,1));


a = boxchart(x_wt_box,wtPlot); hold on;
plot(x_wt,wtPlot,'.','MarkerSize',20,'Color','k');
a.BoxFaceColor = 'k';
a.BoxFaceAlpha = 0;
a.MarkerStyle = 'none';

b = boxchart(x_lof_box,lofPlot); hold on;
plot(x_lof,lofPlot,'.','MarkerSize',20,'Color','r');
b.BoxFaceColor = 'r';
b.BoxFaceAlpha = 0;
b.MarkerStyle = 'none';

[h,p,ci,stats] = ttest2(wtPlot,lofPlot);

txt1 = strcat('p = ',num2str(round(p,1,"significant")));
text(1.35,3.5,txt1,'fontsize',16)
plot([1,2],[3.2,3.2],'-','color','k','linewidth',1);

ylim([0,4]);

xticks([1 2])
xticklabels({'Wildtype','\it{Longfin}'})

yticks([0 1.0 2.0 3.0 4.0])

ylabel({'Length Regen. /','Fgf20a Expressor (\mum)'})

set(gca,'fontsize',16)

saveas(gca,'/Users/ashleyrich/Documents/DiTaliaLab/Experiments_inProg/16Jan26_fgf20a_krt19_merge_redo/21feb26_plots/lreg_GFPplus_LOF_and_WT_stats.jpg')


%"lreg_GFPplus_LOF_and_WT.jpg"

%% Supp. Fig. 7D

%load "fgf20a_timecourse.mat" from "SuppFig7"

close all;

data = exp2DataFull;

fig = figure;

fishSet = [1 2 3 4 11 12 13 14];
raySet = [2 3 6 7];

allData = {};


for p= 1:size(fishSet,2)
    fishHere = fishSet(p);
    for t = 1:size(raySet,2)
        rayHere = raySet(t);
        timeCollect = [];
        countCollect = [];  
        for q = 1:size(exp2DataFull,1)
            

            if data(q,2) == fishHere && data(q,3) == rayHere
                timeHere = data(q,4);
                countHere = data(q,5);
                timeCollect = horzcat(timeCollect,timeHere);
                countCollect = horzcat(countCollect,countHere);
            else 
                continue
            end


        end
        if ~isempty(countCollect)
            oneRayData = {fishHere,rayHere,timeCollect,countCollect};
            allData = vertcat(allData,oneRayData);
            oneRayData = {};
        else
            continue
        end
    end
end

% Plot GFP expression in WT and LOF over time
close all;
f = figure;
allCounts = [];
WTcounts = [];
LOFcounts = [];
for i = 1:size(allData)
    if cell2mat(allData(i,1)) < 10
        cl = [0 0 0];
    else
        cl = [1 0 0];
    end
    cl = horzcat(cl,0.3);
    plot(cell2mat(allData(i,3))./24,cell2mat(allData(i,4))./max(cell2mat(allData(i,4))),'.-','color',cl,'markersize',10,'linewidth',1); hold on;

    %collect data for averaging
    tempTime = cell2mat(allData(i,3))./24;
    tempCount = cell2mat(allData(i,4));
    if size(tempTime,2) == 7
        adjCount = tempCount;
    elseif size(tempTime,2) == 6
        adjCount = horzcat(tempCount,NaN);
    elseif size(tempTime,2) == 5
        adjCount = horzcat(tempCount,NaN,NaN);
    end
        allCounts = vertcat(allCounts,adjCount);
end

    %separate collected values into WT and LOF
    for qq = 1:size(allData)
        if cell2mat(allData(qq,1)) < 10
            WTcounts = vertcat(WTcounts,allCounts(qq,:)./max(allCounts(qq,:)));
        else
            LOFcounts = vertcat(LOFcounts,allCounts(qq,:)./max(allCounts(qq,:)));
        end
    end

    timeVals = [1 2 3 5 7 10 14];
    errorbar(timeVals,nanmean(WTcounts,1),nanstd(WTcounts,0,1)./sqrt(size(WTcounts,1)),'.-','markersize',15,'linewidth',3,'color',[0 0 0])
    errorbar(timeVals,nanmean(LOFcounts,1),nanstd(LOFcounts,0,1)./sqrt(size(LOFcounts,1)),'.-','markersize',15,'linewidth',3,'color',[1 0 0])
ylabel('Normalized #fgf20a Expressors')
xlabel('Time (days post amp.)')
yticks([0:0.2:1])
% ylim([0,1])
% xlim([0,15])
set(gca,'fontsize',16)

%% Supp. Fig. 7E

%load "timeCollect_WTnoAvg.mat", "timeCollect_LOFnoAvg.mat",
%"meanAmpCollect_WTnoAvg.mat", and "meanAmpCollect_LOFnoAvg.mat" from
%"SuppFig7" folder

% Compare decay rate in WT and LOF fish
close all;
f = figure;
plot(timeWTnoAvg./24,meanWTnoAvg,'.','markersize',15,'color',[0 0 0 0.5]); hold on;
plot(timeLOF./24,meanLOF,'.','markersize',15,'color',[1 0 0 0.5]);

ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.StartPoint = [0.6 0.1 0.75];
opts.Upper = [1 1 1];

% Fit model to data.
[fitresultWT, gofWT] = fit(timeWTnoAvg./24, meanWTnoAvg, ft, opts );
[fitresultLOF, gofLOF] = fit(timeLOF./24, meanLOF, ft, opts );

xWT = [48:10:348];

h = plot(fitresultWT);
set(h,'color','k','linewidth',3);
h2 = plot(fitresultLOF);
set(h2,'color','r','linewidth',3);

ylabel('Average ERK Activity @ Amp. Site (AU)')
xlabel('Time (dpa)')
ylim([0.6,1.4])
%yticks([0:0.1:0.5])
set(gca,'fontsize',16)

legend_name={
    strcat("WT 1/b=",num2str(1/fitresultWT.b)),...
        strcat("LOF 1/b=",num2str(1/fitresultLOF.b))};
legend([h,h2], legend_name,...
        'Location','northeast','color','none','box','off');


%% Supp. Fig. 7F

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
% ***builds lagr_mat

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


% fit 3d
% *** generates fit_obj



xx = [lagr_mat_bytime.xi_tip];
% xx = [lagr_mat_bytime.ui];

tt = cell2mat(arrayfun(@(s)repmat(s.t,size(s.xi_tip)), lagr_mat_bytime, uni=0)');

JJ = [lagr_mat_bytime.J_i];
divv = [lagr_mat_bytime.div_i];
vv = [lagr_mat_bytime.vi];
dERKK = [lagr_mat_bytime.dERKidt];
ERKK = [lagr_mat_bytime.ERKi];

fit_obj = fit([xx(:),tt(:)],JJ(:),"cubicinterp",ExtrapolationMethod="linear"); % only work with 2023b
% fit_obj = fit([xx(:),tt(:)],JJ(:),"cubicinterp",ExtrapolationMethod="nearest"); % only work with 2023b

% solve for a spectrum of tau and Lfin
J_cubic = fit_obj;

% build lagrangian matrix
% set initial conditions

t0 = 72;
phi0 = 0.115;
% t0 = 96;
% phi0 = 0.24;
% t0 = 24;
% phi0 = 0.025;


xsize = 500;

% initial condition for rho (ERK gradient)

aHere = (lagr_mat_bytime(1).ERKi(end)-lagr_mat_bytime(1).ERKi(1))./lagr_mat_bytime(1).xi(end); % (E2-E1)/Length
bHere = lagr_mat_bytime(1).ERKi(1); % E1



% construct ode solver

i_array = 0:xsize;
x0 = i_array.*phi0./xsize;
rho0 = x0.*aHere+bHere;

N = numel(i_array);

% solve ode
% wt
tic
tau = 3.7*24;
[t_wt,y_wt]=ode45(@(t,y)growthmodel_lagrangian_characteristics(t,y,J_cubic,tau),[t0:12:336],[rho0 x0]); % y: 1st column is phi, 2nd is xi, 3rd is rhoi
phi_wt = y_wt(:,end);
phi_fin_wt = phi_wt(end);
toc

tic
phi_fin = [];
for tau = 1:0.1:10
    tau
    [t_here,y_here]=ode45(@(t,y)growthmodel_lagrangian_characteristics(t,y,J_cubic,tau.*24),[t0:6:28*24],[rho0 x0]); % y: 1st column is phi, 2nd is xi, 3rd is rhoi
    phi_here = y_here(:,end);
    phi_fin_here = phi_here(end);
    phi_fin = [phi_fin phi_fin_here];

end
toc

% plot
f = figure;
plot(1:0.1:10, phi_fin./phi_fin_wt, 'LineWidth',2); hold on;
y1 = phi_fin(1:0.1:10==3.7)./phi_fin_wt;
line([3.7 3.7 0],[0 y1 y1],'Color','k','LineWidth',1.5)
y2 = phi_fin(1:0.1:10==6)./phi_fin_wt;
line([6 6 0],[0 y2 y2],'Color','r','LineWidth',1.5)
xlabel('\tau (day)')
ylabel('Lfin/Lfin_{wt}')
xticks(0:2:10)
config_plot(f)














































