function [analysis_mat] = generate_analysis_mat(st_dir,tgmmFolder,dbins,fish_modifier,L_amp, gem_thres, qc)
%GENERATE_ANALYSIS_MAT Summary of this function goes here
%   fish_modifier is used to distinguish fish from different expriments but
%   have the same fish number
arguments
    st_dir
    tgmmFolder
    dbins = 50
    fish_modifier = 0
    L_amp = []
    gem_thres = 1.5;
    qc = false;
end

analysis_mat = [];
for i =1:numel(st_dir)
    disp('------------');
    display([st_dir(i).folder filesep st_dir(i).name]);
    
    % create analysis_mat struct for each single experiments
    load([st_dir(i).folder filesep st_dir(i).name]);  %load myScale

    if qc
        [myScale,QC_logical] = segmentation_qualitycontrol(myScale);
    
        % QC KTR, edited 020823, this part now included in segmentation_qualitycontrol      
%         QC_logical = myScale.(tgmmFolder).CoV_nuccyt<0.4 & myScale.(tgmmFolder).median_nuccyt>=10;
%         for j = fieldnames(myScale.(tgmmFolder))'
%             if ismember(numel(QC_logical),size(myScale.(tgmmFolder).(j{:})))
%                 myScale.(tgmmFolder).(j{:}) = myScale.(tgmmFolder).(j{:})(QC_logical,:);
%             end
%         end
    end

    % Load and remove NaN/Inf; invert gem to get nuc/cyt if needed
%     ampPointLabel = 'amputationPoint';
%     endPointLabel = 'endPoint';
%     ampPoint = myScale.(tgmmFolder).(ampPointLabel);
%     endPoint = myScale.(tgmmFolder).(endPointLabel);
    % Rotate centers and amp/endPoints
    ampPointrot = myScale.ampPoint.*myScale.pxsize;
    endPointrot = myScale.endPoint.*myScale.pxsize;
    ccrot = myScale.(tgmmFolder).ccrot;
    ccrot(:,1) = endPointrot - myScale.(tgmmFolder).ccrot(:,1).*myScale.pxsize;
%     [ccrot,ampPointrot,endPointrot] = process_coords(cc,ampPoint,endPoint,pxsize);

    % ktr
    if isfield(myScale.(tgmmFolder),'ktr')
        ktr=myScale.(tgmmFolder).ktr_nonoverlap;
        %idxs=FUCCI>1; %All the FUCCI values that are not infinity
        idxs=~isnan(ktr) & ~isinf(ktr); %All the FUCCI values that are not NaN
        %idxs=~isinf(idxs);
        ccrot=ccrot(idxs,:);%
        ktr=ktr(idxs);%
        ktrAve = myScale.(tgmmFolder).ktrAve;
    else
        idxs = ones(size(ccrot));
        ktr = [];
        ktrAve = [];
    end
    
    % gem
    if isfield(myScale.(tgmmFolder),'gem')
        gem=myScale.(tgmmFolder).gem(idxs);
    elseif isfield(myScale.(tgmmFolder),'GEM')
        gem=myScale.(tgmmFolder).GEM(idxs);
    elseif isfield(myScale.(tgmmFolder),'greenCytNuc')
        gem=1./myScale.(tgmmFolder).greenCytNuc(idxs);
    else
        gem = [];
    end

    %BoneRemoved = myScale.TGMM_hypo_eq_ch3.BoneRemoved;
    

    
    single_exp = [];
    single_exp.name = st_dir(i).name;
    single_exp.name = ['fish',num2str(myScale.fish+fish_modifier),...
        '_ray',num2str(myScale.ray),'_hpa',num2str(myScale.hpp)];
    single_exp.fish = myScale.fish+fish_modifier;
    single_exp.ray = myScale.ray;
    single_exp.hpa = myScale.hpp;
    single_exp.hpaTrue = myScale.hppTrue;
    if ~isempty(L_amp)
        single_exp.L_amp = L_amp.Lamp(L_amp.fish==myScale.fish & L_amp.ray==myScale.ray);
%         single_exp.L_amp = L_amp.AmoutnAmputated(L_amp.Fish==myScale.fish & L_amp.Ray==myScale.ray); % old way of storing
    else
        single_exp.L_amp = [];
    end
    single_exp.ktrAve = ktrAve;
    single_exp.averageH2ANuc = myScale.(tgmmFolder).averageH2ANuc;
    single_exp.ratioH2ANucCyt_nonoverlap = myScale.(tgmmFolder).ratioH2ANucCyt_nonoverlap;

    single_exp = generate_oneray(single_exp,ktr,gem,ccrot,ampPointrot,endPointrot,dbins,gem_thres);

    % crosssection Area
    crossArea = [];
    volumeBinned = [];
    if isfield(myScale,'crossArea')
        crossArea = flip(myScale.crossArea); % zero at distal tip
        areaHere = crossArea(crossArea>0);
        xHere = (0:(numel(areaHere)-1)).*myScale.pxsize;

        bins = 0:dbins:(endPointrot(1)-ampPointrot(1));
        for bin = 1:(numel(bins)-1)
            volumeHere = sum(areaHere(xHere>=bins(bin)&xHere<bins(bin+1)),'all').*(myScale.pxsize)^3; %um^3
            volumeBinned = [volumeBinned, volumeHere];
        end
    end
    single_exp.crossArea = crossArea;
    single_exp.volumeBinned = volumeBinned;

    % crosssection Area filled
    crossArea_fill = [];
    volumeBinned_fill = [];
    if isfield(myScale,'crossArea_fill')
        crossArea_fill = flip(myScale.crossArea_fill); % zero at distal tip
        areaHere = crossArea_fill(crossArea_fill>0);
        xHere = (0:(numel(areaHere)-1)).*myScale.pxsize;

        bins = 0:dbins:(endPointrot(1)-ampPointrot(1));
        for bin = 1:(numel(bins)-1)
            volumeHere = sum(areaHere(xHere>=bins(bin)&xHere<bins(bin+1)),'all').*(myScale.pxsize)^3; %um^3
            volumeBinned_fill = [volumeBinned_fill, volumeHere];
        end
    end
    single_exp.crossArea_fill = crossArea_fill;
    single_exp.volumeBinned_fill = volumeBinned_fill;


    if numel(single_exp.nucleiALL)==numel(single_exp.volumeBinned)
        single_exp.densityBinned = single_exp.nucleiALL./single_exp.volumeBinned;
        single_exp.densityBinned_fill = single_exp.nucleiALL./single_exp.volumeBinned_fill;
    else
        single_exp.densityBinned  = [];
        single_exp.densityBinned_fill = [];
    end


    
    analysis_mat = [analysis_mat single_exp];
end
end

