%% preparation
%%%%
% this pipeline is compatible with MATLAB2021b
clear;clc;close all;
% addpath(genpath('D:\BIO\PhD\ditalia\zebrafish\github\ditalia-zebrafish\code\functions\'))
addpath(genpath('./functions/'))

%% Initialization
paths =[];
paths.expName = '19Aug22_tph1b_H2A10b_tph1b_ERK30';
% paths.expName = '23Sep22_tph1b_H2A10b_tph1b_ERK40';
% paths.expName = '24Nov22_MEK_inh_tph1b_H2A10b_tph1b_ERK30';
% paths.expName = '12Jan23_MEK_inh_tph1b_H2A10b_tph1b_ERK40';
paths.expName = '08Feb23_MEK_inh_tph1b_H2A10b_ERK40_GEM2';
paths.expFolder = ['D:\BIO\PhD\ditalia\zebrafish\github\ditalia-zebrafish\data\fibroblasts\' filesep paths.expName];
% paths.expFolder = 'E:\19Aug22_tph1b_H2A10b_tph1b_ERK30\';

paths.objFolder = [paths.expFolder filesep 'objects\'];
paths.matFolder = [paths.expFolder filesep 'analysis/'];
% paths.lampFolder = 'E:\19Aug22_tph1bRFPtph1b_ERK_dissecting/';
paths.lampFolder = 'D:\BIO\PhD\ditalia\zebrafish\github\ditalia-zebrafish\data\length_measurement';
% paths.lampFolder = paths.expFolder;
mkdir(paths.matFolder);

% Piexel size
dbins = 50;

fish_modifier= 0;

osteoblast = 0;

L_amp = [];
%% Calculate L_amp
excel_all = dir([paths.lampFolder,filesep,'*xlsx*']);


L_amp_cell = {};
%
excel1 = readtable([paths.lampFolder filesep excel_all(1).name]);
L_amp_tbl = excel1(1:24,[1:2 5]);
% L_amp_tbl.Fish = cellfun(@(x)str2double(x(end)), L_amp_tbl.fish);
% L_amp_tbl.Ray = cellfun(@(x)str2double(x(end)), L_amp_tbl.ray);
L_amp_tbl.Fish = L_amp_tbl.fish;
L_amp_tbl.Ray = L_amp_tbl.ray;
L_amp_tbl.AmoutnAmputated = cellfun(@(x)x.*96000./281.33, num2cell(L_amp_tbl.AmountCut));
L_amp_cell = [L_amp_cell,{L_amp_tbl}];
%
% excel2 = readtable([paths.lampFolder filesep excel_all(2).name]);
% L_amp_tbl = excel2(2:15,[1:2 10]);
% L_amp_tbl.Fish = cellfun(@(x)str2double(x(end)), L_amp_tbl.Fish);
% L_amp_tbl.Ray = cellfun(@(x)str2double(x(end)), L_amp_tbl.Ray);
% L_amp_tbl.AmoutnAmputated = cellfun(@(x)str2double(x).*96000./281.33, L_amp_tbl.AmoutnAmputated);
% L_amp_cell = [L_amp_cell,{L_amp_tbl}];


L_amp = L_amp_cell{1};
%% Generate analysis matrix
if osteoblast
    tgmmFolder = 'TGMM_hypo_eq_ch2';
else
    tgmmFolder = 'TGMM_equalized_ch2';
end



fish=''; % put the number of the fish/scale/hpp to process or '' to choose them all
scale=''; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*ray' num2str(scale) '*' num2str(hpp) 'hpa' '.mat']);

% generate analysis mat
analysis_mat = generate_analysis_mat(st_dir,tgmmFolder,dbins,...
    fish_modifier,L_amp, true);


[~,idx] = unique([[analysis_mat.fish];[analysis_mat.ray];[analysis_mat.hpa]]','rows');
analysis_mat = analysis_mat(idx);
%% calculate growth rate
analysis_mat_backup = analysis_mat;
analysis_mat = [];
for fish_ray = unique([analysis_mat_backup.fish;analysis_mat_backup.ray]','rows','stable')'
        disp('------------');
        display(['fish' num2str(fish_ray(1)) '_ray' num2str(fish_ray(2))]);
        
        single_ray = [];
        single_ray.fish = fish_ray(1);
        single_ray.ray = fish_ray(2);

        fishNray_data = analysis_mat_backup( ...
            [analysis_mat_backup.fish]==fish_ray(1) & ...
            [analysis_mat_backup.ray]==fish_ray(2));

        dLdt = num2cell(...
            find_derivative(...
            [fishNray_data.hpaTrue], [fishNray_data.L_reg], [fishNray_data.hpaTrue])...
            );

        [fishNray_data.dLdt] = dLdt{:};


        analysis_mat = [analysis_mat,fishNray_data];
end
%% save
% Save analysis_mat
save([paths.matFolder filesep 'analysis_mat.mat'],'analysis_mat');
disp([paths.matFolder filesep 'analysis_mat.mat' '  Saved']);

%% load data
load([paths.matFolder filesep 'analysis_mat.mat'],'analysis_mat')


%%
%%
%%
%% plot ktr fit
paths.matFolder = [paths.expFolder '/analysis/'];
paths.ktrFolder = [paths.expFolder '/fitting/ktr/'];

mkdir(paths.ktrFolder);

% load analysis matrix
load([paths.matFolder filesep 'analysis_mat.mat']);

plot_fit_ktr(analysis_mat,paths.ktrFolder)


%% plot ktr fit bin
paths.matFolder = [paths.expFolder '/analysis/'];
paths.ktrFolder = [paths.expFolder '/fitting/ktr_binned/'];

mkdir(paths.ktrFolder);

% load analysis matrix
load([paths.matFolder filesep 'analysis_mat.mat']);

plot_fit_ktr(analysis_mat,paths.ktrFolder,'_binned')







%% plot fit gem
paths.matFolder = [paths.expFolder '/analysis/'];
paths.gemFolder = [paths.expFolder '/fitting/gem/'];

mkdir(paths.gemFolder);

% load analysis matrix
load([paths.matFolder filesep 'analysis_mat.mat']);


plot_fit_gem(analysis_mat,paths.gemFolder)


%% plot fit gem weightN
paths.matFolder = [paths.expFolder '/analysis/'];
paths.gemFolder = [paths.expFolder '/fitting/gem_weightN/'];

mkdir(paths.gemFolder);

% load analysis matrix
load([paths.matFolder filesep 'analysis_mat.mat']);


plot_fit_gem(analysis_mat,paths.gemFolder,'_weightN')
%% plot fit gem weightVar
paths.matFolder = [paths.expFolder '/analysis/'];
paths.gemFolder = [paths.expFolder '/fitting/gem_weightVar/'];

mkdir(paths.gemFolder);

% load analysis matrix
load([paths.matFolder filesep 'analysis_mat.mat']);


plot_fit_gem(analysis_mat,paths.gemFolder,'_weightVar')