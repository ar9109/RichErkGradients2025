function [analysis_mat] = generate_analysis_mat(st_dir,pxsize,dbins,fish_modifier,L_amp)
%GENERATE_ANALYSIS_MAT Summary of this function goes here
%   fish_modifier is used to distinguish fish from different expriments but
%   have the same fish number
arguments
    st_dir
    pxsize
    dbins = 50
    fish_modifier = 0
    L_amp = []
end

analysis_mat = [];
for i =1:numel(st_dir)
    disp('------------');
    display([st_dir(i).folder filesep st_dir(i).name]);
    
    % create analysis_mat struct for each single experiments
    
    load([st_dir(i).folder filesep st_dir(i).name]);
    
    % Load and remove NaN/Inf; invert gem to get nuc/cyt if needed
    cc=myScale.TGMM_hypo_eq_ch2.centers;
    ktr=myScale.TGMM_hypo_eq_ch2.ktr;
    %idxs=FUCCI>1; %All the FUCCI values that are not infinity
    idxs=~isnan(ktr) & ~isinf(ktr); %All the FUCCI values that are not NaN
    %idxs=~isinf(idxs);
    cc=cc(idxs,:);%
    ktr=ktr(idxs);%
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
    
    single_exp = [];
    single_exp.name = st_dir(i).name;
    single_exp.name = ['fish',num2str(myScale.fish+fish_modifier),...
        '_ray',num2str(myScale.ray),'_hpa',num2str(myScale.hpp)];
    single_exp.fish = myScale.fish+fish_modifier;
    single_exp.ray = myScale.ray;
    single_exp.hpa = myScale.hpp;
    single_exp.hpaTrue = myScale.hppTrue;
    if ~isempty(L_amp)
        single_exp.L_amp = L_amp.AmoutnAmputated(L_amp.Fish==myScale.fish & L_amp.Ray==myScale.ray);
    else
        single_exp.L_amp = [];
    end
    single_exp = generate_oneray(single_exp,ktr,gem,ccrot,ampPointrot,endPointrot,dbins);

    
    analysis_mat = [analysis_mat single_exp];
end
end

