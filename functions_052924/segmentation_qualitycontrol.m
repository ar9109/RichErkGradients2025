function [myScale_new,QC_logical] = segmentation_qualitycontrol(myScale,option)
arguments
    myScale;
    option.volume = 0;
    option.averageH2ANuc = 0;
    option.volume_high = 0.9;
%     option.ktrAve = 0.25;
    option.CoV_nuccyt = 0.4;
    option.median_nuccyt = 10;
end
%SEGMENTATION_QUALITYCONTROL Summary of this function goes here
%   Detailed explanation goes here
myScale_new = myScale;

field_names = fieldnames(myScale);
if sum(contains(field_names,'TGMM'))>1
    error('More than one TGMM field detected.')
else
    TGMM_field = field_names{contains(field_names,'TGMM')};
end

QC_logical = myScale.(TGMM_field).volume>=quantile(myScale.(TGMM_field).volume,option.volume) &...
    myScale.(TGMM_field).averageH2ANuc>=quantile(myScale.(TGMM_field).averageH2ANuc,option.averageH2ANuc) &...
    myScale.(TGMM_field).volume<=quantile(myScale.(TGMM_field).volume,option.volume_high) &...
    myScale.(TGMM_field).CoV_nuccyt<option.CoV_nuccyt &...
    myScale.(TGMM_field).median_nuccyt>=option.median_nuccyt;

% trim myScale
for i = fieldnames(myScale.(TGMM_field))'
    if ismember(numel(QC_logical),size(myScale.(TGMM_field).(i{:})))
        myScale_new.(TGMM_field).(i{:}) = myScale.(TGMM_field).(i{:})(QC_logical,:);
    end
end
end

